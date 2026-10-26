module InterlinkTerminalTests

# Interlink terminals from mounted hardware: terminal numbers, endpoints given as
# terminals, directed links, rule-based registration, the power a link draws, a
# controller's levels for a link-driven terminal, and shedding that closes links.

using Test
using Logging
using LinearAlgebra
using StaticArrays
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PLANET = make_no_gram_planet(:earth)
const RADIUS = PLANET.Rp_e + 700e3
const SPEED = sqrt(PLANET.μ / RADIUS)

# A 100 kg bus on a circular orbit, `y` metres along track from the reference
# point, with laser terminals of distinct powers and an optional battery.
function terminal_sat(index; y=0.0, terminals=1, n_terminal=0, capacity_wh=10.0, soc=0.8, efficiency=0.5)
    bus = Link(root=true, m=100.0)
    ic = CartesianInitialCondition(SVector(RADIUS, y, 0.0), SVector(0.0, SPEED, 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, initial_condition=ic, id=index,
        n_terminal=n_terminal)
    for k in 1:terminals
        add_actuator!(sc, bus, LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=50.0 + k, magnification=100.0,
            range_m=200e3, efficiency=efficiency))
    end
    capacity_wh === nothing || add_power!(sc, bus, Battery(capacity_wh=capacity_wh, initial_soc=soc, recovery_soc=0.9))
    return sc
end

# Holds levels for the mounted actuators of spacecraft 1.
mutable struct HoldLevels <: AbstractControlEffectorModel
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::HoldLevels) = :mounted
SpaceAGORA.actuator_levels(c::HoldLevels, i::Int) = i == 1 ? c.levels : nothing
SpaceAGORA.update_control!(::HoldLevels, u, p, t::Float64, i::Int64) = nothing

# A directed kind that records the endpoint order its hooks receive.
struct OneWayLink <: AbstractInterlink end
SpaceAGORA.link_range_m(::OneWayLink) = 1e9
SpaceAGORA.link_directed(::OneWayLink) = true
const SEEN = Tuple{Int, Int}[]
SpaceAGORA.link_available(::OneWayLink, a, b, ctx) = (push!(SEEN, (a.spacecraft_index, b.spacecraft_index)); true)

# A kind no laser terminal carries.
struct RadioLink <: AbstractInterlink end
SpaceAGORA.link_range_m(::RadioLink) = 1e9

function run_args(spacecraft, model; duration=60.0, control=ControlModel((), Float64[]), orientation=false)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, generate_plots=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, mission_time=duration,
            orientation_sim=orientation, data_rate=1.0),
        environment_model=make_no_gram_environment(planet=PLANET),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]), navigation_model=NavigationModel((), Float64[]),
        control_model=control, initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-8, dt_max_orbit=1.0),
        solver_config=SM.SolverConfig(solver_mode=:tsit5), interlink_model=model,
        scheduling_policy_model=SchedulingPolicyModel(SemimajorAxisRatePolicy()))
end

# Spacecraft 2 trails spacecraft 1, so a laser link pushes 1 forward and the
# default policy (targeting spacecraft 1) selects it.
function linked_pair(; kwargs...)
    spacecraft = [terminal_sat(1; kwargs...), terminal_sat(2; y=-1000.0, kwargs...)]
    model = InterlinkModel()
    first_terminal = only(mounted_link_terminals(spacecraft[1])).terminal
    register_link!(model, spacecraft, (1, first_terminal), (2, 1), LaserMomentumLink(power_w=40.0, magnification=100.0, range_m=200e3))
    return spacecraft, model
end

@testset "Interlink terminals" begin
    @testset "terminal numbers" begin
        sc = terminal_sat(1; terminals=2, n_terminal=2)
        @test terminal_count(sc) == 4
        mounted = mounted_link_terminals(sc)
        @test [entry.number for entry in mounted] == [3, 4] && [entry.position for entry in mounted] == [1, 2]
        @test terminal_number(sc, mounted[2].terminal) == 4
        @test_throws ArgumentError terminal_number(sc, LaserTerminal(boresight=[1.0, 0.0, 0.0], max_power_w=1.0))
        twin = terminal_sat(1; terminals=0)
        same = LaserTerminal(boresight=[1.0, 0.0, 0.0], max_power_w=5.0)
        add_actuator!(twin, twin.root, same); add_actuator!(twin, twin.root, same)
        @test_throws ArgumentError terminal_number(twin, same)                    # equal terminals: give the number
        @test terminal_count((n_terminal=3,)) == 3
        @test link_terminal(mounted[1].terminal) && !link_terminal(Thruster())
    end

    @testset "hardware endpoints" begin
        spacecraft = [terminal_sat(1; n_terminal=1), terminal_sat(2; n_terminal=1)]
        terminal = only(mounted_link_terminals(spacecraft[1])).terminal
        by_object = register_link!(InterlinkModel(), spacecraft, (1, terminal), (2, 2), LaserMomentumLink())
        by_number = register_link!(InterlinkModel(), spacecraft, (1, 2), (2, 2), LaserMomentumLink())
        @test by_object == by_number == ((1, 2), (2, 2))
        @test_throws ArgumentError register_link!(InterlinkModel(), spacecraft, (1, 3), (2, 1), LaserMomentumLink())
        @test_throws ArgumentError register_link!(InterlinkModel(), spacecraft, (1, 2), (2, 1), RadioLink())
        @test register_link!(InterlinkModel(), spacecraft, (1, 1), (2, 1), RadioLink()) == ((1, 1), (2, 1)) # virtual carries any
        @test terminal_supports(terminal, LaserMomentumLink()) && !terminal_supports(terminal, RadioLink())
        @test link_terminal_levels(LaserMomentumLink(power_w=7.0)) == (7.0, 7.0)
        @test link_terminal_levels(RadioLink()) == (0.0, 0.0)
    end

    @testset "directed links" begin
        spacecraft = [(n_terminal=1,), (n_terminal=1,)]
        model = InterlinkModel(forbidden_pairs=[(1, 2)])
        @test register_link!(model, spacecraft, (2, 1), (1, 1), OneWayLink()) == ((2, 1), (1, 1))  # transmitter first
        @test register_link!(InterlinkModel(), spacecraft, (2, 1), (1, 1), LaserMomentumLink()) == ((1, 1), (2, 1))
        state = (sc=[(pos=[RADIUS, 0.0, 0.0],), (pos=[RADIUS, 10.0, 0.0],)],)
        @test update_availability!(model, spacecraft, state) == [false]         # the forbidden pair, either way round
        open_model = InterlinkModel()
        register_link!(open_model, spacecraft, (2, 1), (1, 1), OneWayLink())
        empty!(SEEN)
        @test update_availability!(open_model, spacecraft, state) == [true]
        @test SEEN == [(2, 1)]                                                   # the hooks see the transmitter as a
        @test update_availability!(open_model, spacecraft, state; terminal_available=[true, false]) == [false]
    end

    @testset "rule-based registration" begin
        spacecraft = [(n_terminal=2,) for _ in 1:3]
        @test all_pairs(spacecraft) == [(1, 2), (1, 3), (2, 3)]
        model = InterlinkModel()
        keys = register_links!(model, spacecraft, all_pairs(spacecraft), LaserMomentumLink())
        @test keys == [((1, 1), (2, 1)), ((1, 2), (3, 1)), ((2, 2), (3, 2))]     # terminals in turn
        built = InterlinkModel()
        register_links!(built, spacecraft, [(1, 2), (2, 3)], (i, j) -> LaserMomentumLink(range_m=1000.0 * i))
        @test [link.range_m for link in built.graph.edge_properties] == [1000.0, 2000.0]
        chosen = InterlinkModel()
        @test register_links!(chosen, spacecraft, [(1, 3)], LaserMomentumLink(); terminals=(i, j) -> (2, 2)) ==
            [((1, 2), (3, 2))]
        # Wrapping: three links on a spacecraft with two terminals share its first again.
        ring = InterlinkModel()
        hub = [(n_terminal=2,), (n_terminal=1,), (n_terminal=1,), (n_terminal=1,)]
        @test first.(register_links!(ring, hub, [(1, 2), (1, 3), (1, 4)], LaserMomentumLink())) == [(1, 1), (1, 2), (1, 1)]
        lasers = [terminal_sat(1), terminal_sat(2)]
        @test_throws ArgumentError register_links!(InterlinkModel(), lasers, [(1, 2)], RadioLink())
        @test_throws ArgumentError register_links!(InterlinkModel(), spacecraft, [(1, 2)], (i, j) -> 1.0)
    end

    @testset "automatic registration skips virtual hardware endpoints" begin
        for (label, link) in (("mounted recoil", LaserBeamLink(power_w=5.0, magnification=1.0, range_m=1e6, recoil=:mounted)),
                              ("pointing", LaserBeamLink(power_w=5.0, magnification=1.0, range_m=1e6,
                                  pointing=TerminalPointing(field_of_regard_rad=pi / 2))))
            @testset "$label" begin
                spacecraft = [terminal_sat(i; terminals=2, n_terminal=1) for i in 1:4]
                model = InterlinkModel()
                keys = register_links!(model, spacecraft, [(1, 2), (1, 3), (1, 4)], link)
                @test keys == [((1, 2), (2, 2)), ((1, 3), (3, 2)), ((1, 2), (4, 2))]
                @test length(model.graph.edges) == 3
                @test all(pair -> all(endpoint -> endpoint[2] > 1, pair), keys)
                @test_throws ArgumentError register_links!(InterlinkModel(),
                    [(n_terminal=1,), (n_terminal=1,)], [(1, 2)], link)
                # Explicit endpoint selection keeps its existing validation.
                @test_throws ArgumentError register_links!(InterlinkModel(), spacecraft,
                    [(1, 2)], link; terminals=(i, j) -> (1, 1))
            end
        end
        @testset "Shannon radio budget" begin
            spacecraft = [terminal_sat(i; terminals=1, n_terminal=1) for i in 1:2]
            for sc in spacecraft
                add_actuator!(sc, sc.root, RadioTerminal(boresight=[0.0, 1.0, 0.0],
                                                       max_power_w=5.0, gain_dbi=6.0, efficiency=0.5))
            end
            link = DataLink(transmit_power_w=1.0,
                budget=ShannonBudget(bandwidth_hz=1e6, frequency_hz=2.2e9), range_m=1e6)
            model = InterlinkModel()
            # Skip both the virtual endpoint and the incompatible laser hardware.
            @test register_links!(model, spacecraft, [(2, 1)], link) == [((2, 3), (1, 3))]
            @test length(model.graph.edges) == 1
            @test_throws ArgumentError register_links!(InterlinkModel(),
                [(n_terminal=1,), (n_terminal=1,)], [(1, 2)], link)
        end
        # Links without hardware requirements retain virtual terminal numbering.
        virtual = [(n_terminal=1,), (n_terminal=1,)]
        @test register_links!(InterlinkModel(), virtual, [(1, 2)], LaserBeamLink(power_w=5.0, magnification=1.0, range_m=1e6)) == [((1, 1), (2, 1))]
    end

    @testset "automatic registration respects directed hardware roles" begin
        spacecraft = [terminal_sat(i) for i in 1:2]
        for sc in spacecraft
            add_actuator!(sc, sc.root, PowerReceiver(normal=[0.0, -1.0, 0.0],
                                                   radius_m=0.25, efficiency=0.4))
        end
        link = PowerBeamLink(power_w=5.0, range_m=1e6)
        for pair in ((1, 2), (2, 1))
            @testset "direction $pair" begin
                model = InterlinkModel()
                @test register_links!(model, spacecraft, [pair], link) ==
                    [((pair[1], 1), (pair[2], 2))]
                @test length(model.graph.edges) == 1
            end
        end
        # A factory uses the same role-aware selection for each pair.
        @test register_links!(InterlinkModel(), spacecraft, [(1, 2), (2, 1)],
            (i, j) -> link) == [((1, 1), (2, 2)), ((2, 1), (1, 2))]
        @test_throws ArgumentError register_links!(InterlinkModel(), spacecraft,
            [(1, 2)], link; terminals=(i, j) -> (2, 1))
        lasers = [terminal_sat(i) for i in 1:2]
        @test_throws ArgumentError register_links!(InterlinkModel(), lasers, [(1, 2)], link)
        receivers = [terminal_sat(i; terminals=0) for i in 1:2]
        for sc in receivers
            add_actuator!(sc, sc.root, PowerReceiver(normal=[0.0, -1.0, 0.0],
                                                   radius_m=0.25, efficiency=0.4))
        end
        @test_throws ArgumentError register_links!(InterlinkModel(), receivers, [(1, 2)], link)
    end

    @testset "pairs by plane and slot" begin
        factory(plane, slot, index) = (bus = Link(root=true, m=100.0);
            SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, n_terminal=4))
        walker = walker_constellation(factory; total=12, planes=3, phasing=1, a=RADIUS, i=53)
        plane_of(index) = walker[index].constellation_slot.plane
        slot_of(index) = walker[index].constellation_slot.slot
        same = same_plane_pairs(walker)
        @test length(same) == 18 && all(((i, j),) -> i < j && plane_of(i) == plane_of(j), same)
        ring = ring_pairs(walker)
        @test length(ring) == 12 && all(index -> count(pair -> index in pair, ring) == 2, eachindex(walker))
        @test all(((i, j),) -> plane_of(i) == plane_of(j) && mod(slot_of(j) - slot_of(i), 4) in (1, 3), ring)
        grid = grid_pairs(walker)
        across = setdiff(grid, ring)
        @test length(grid) == 24 && length(across) == 12
        @test all(((i, j),) -> slot_of(i) == slot_of(j) && plane_of(i) != plane_of(j), across)
        @test all(index -> count(pair -> index in pair, grid) == 4, eachindex(walker))
        open_grid = grid_pairs(walker; wrap_planes=false)
        @test length(open_grid) == 20 && !any(((i, j),) -> Set((plane_of(i), plane_of(j))) == Set((1, 3)), setdiff(open_grid, ring))
        # Spacecraft without a slot take part in no pair; small planes.
        plain = [walker; SpacecraftModel()]
        @test same_plane_pairs(plain) == same && ring_pairs(plain) == ring && grid_pairs(plain) == grid
        two = build_constellation(factory, [OrbitalPlane(a=RADIUS, i=10, Ω=0, M=[0, 180]), OrbitalPlane(a=RADIUS, i=10, Ω=90, M=[0])])
        @test ring_pairs(two) == [(1, 2)] && same_plane_pairs(two) == [(1, 2)]
        @test grid_pairs(two) == [(1, 2), (1, 3)]
        @test isempty(ring_pairs([(n_terminal=1,), (n_terminal=1,)]))
        # Two spacecraft in one slot have no single neighbour: refused.
        twin = [walker[1], factory(1, 1, 13)]
        twin[2].constellation_slot = walker[1].constellation_slot
        for rule in (same_plane_pairs, ring_pairs, grid_pairs)
            @test_throws ArgumentError rule(twin)
        end
        # The +grid on four terminals each: every spacecraft's four links take its four terminals.
        model = InterlinkModel()
        keys = register_links!(model, walker, grid, LaserMomentumLink())
        @test length(keys) == 24
        @test all(index -> sort([e[2] for key in keys for e in key if e[1] == index]) == [1, 2, 3, 4], eachindex(walker))
    end

    @testset "a hardware link draws its power" begin
        spacecraft, model = linked_pair()
        solution = run_simulation(run_args(spacecraft, model); return_solution=true)
        ran = solution.prob.p.args.interlink_model
        @test ran.terminal_endpoints == [[(1, 1, 1)], [(1, 2, 1)]]
        @test all(sample -> sample.active == [((1, 1), (2, 1))], ran.history)
        # Each end draws 40 W at 50 % efficiency for 60 s; the force is the link's.
        for i in 1:2
            @test solution.u[1].sc[i].battery_wh - solution.u[end].sc[i].battery_wh ≈ 80.0 * 60.0 / 3600.0 rtol=1e-9
        end
        @test norm(solution.u[end].sc[1].interlink_dv) ≈ 100.0 * 40.0 / 299_792_458.0 / 100.0 * 60.0 rtol=1e-6
        result = run_simulation(run_args(linked_pair()...); return_results=true)
        @test all(==(40.0), result.table.sc1_actuator_level_1) && all(==(40.0), result.table.sc2_actuator_level_1)
        @test all(≈(80.0), result.table.sc1_actuator_power)
        @test all(≈(80.0), result.table.sc1_power_load_w)
    end

    @testset "a virtual-terminal link draws nothing" begin
        spacecraft = [terminal_sat(1; n_terminal=1), terminal_sat(2; y=-1000.0, n_terminal=1)]
        model = InterlinkModel()
        register_link!(model, spacecraft, (1, 1), (2, 1), LaserMomentumLink(power_w=40.0))
        solution = run_simulation(run_args(spacecraft, model); return_solution=true)
        @test isempty(solution.prob.p.args.interlink_model.terminal_endpoints[1])
        @test solution.u[end].sc[1].battery_wh == solution.u[1].sc[1].battery_wh
        @test norm(solution.u[end].sc[1].interlink_dv) > 0
    end

    @testset "a controller's levels for a link-driven terminal apply nothing" begin
        # Unlinked, the controller's 30 W pushes spacecraft 1 back (recoil) and draws 60 W.
        free = [terminal_sat(1), terminal_sat(2; y=-1000.0)]
        held = HoldLevels([30.0])
        control = ControlModel(control_effectors=(held,), control_rates=[1.0])
        unlinked = run_simulation(run_args(free, nothing; control, orientation=true); return_solution=true)
        @test unlinked.u[1].sc[1].battery_wh - unlinked.u[end].sc[1].battery_wh ≈ 60.0 * 60.0 / 3600.0 rtol=1e-9
        recoil = unlinked.u[end].sc[1].vel - unlinked.u[1].sc[1].vel
        @test norm(recoil) > 0
        # Linked, the link drives it: its 40 W, no recoil, and setup says so once.
        spacecraft, model = linked_pair()
        logger = Test.TestLogger(min_level=Logging.Info)
        linked = Logging.with_logger(logger) do
            run_simulation(run_args(spacecraft, model; control=ControlModel(control_effectors=(HoldLevels([30.0]),),
                control_rates=[1.0]), orientation=true); return_solution=true)
        end
        @test count(entry -> occursin("Interlinks drive their terminal hardware for the whole run (spacecraft 1 at mounted actuator positions 1; spacecraft 2", string(entry.message)), logger.logs) == 1
        @test linked.u[1].sc[1].battery_wh - linked.u[end].sc[1].battery_wh ≈ 80.0 * 60.0 / 3600.0 rtol=1e-9
        @test SE._mounted_control_effect(HoldLevels([30.0]), linked.u[end].sc[1], linked.prob.p, 1, 60.0)[1] ==
            zeros(3)                                                            # no recoil from the terminal
    end

    @testset "shedding closes the link" begin
        # A 0.1 Wh battery feeding 80 W empties in 4.5 s; tier 1 then sheds the
        # laser terminal, whose link closes until the battery recovers (it cannot:
        # there are no arrays).
        spacecraft, model = linked_pair(capacity_wh=0.5, soc=0.2)
        solution = run_simulation(run_args(spacecraft, model; duration=20.0); return_solution=true)
        history = solution.prob.p.args.interlink_model.history
        closed = findfirst(sample -> isempty(sample.active), history)
        @test closed !== nothing && 4.0 < history[closed].time < 6.0
        @test all(sample -> isempty(sample.active), history[closed:end])
        @test all(sample -> !isempty(sample.active), history[1:closed - 1])
        @test solution.prob.p.args.dynamics_model.spacecraft[1].root.batteries[1].shed_tier >= 1
        flags = SM.InterlinkTerminals.terminal_availability(
            solution.prob.p.args.interlink_model, solution.prob.p.args.dynamics_model.spacecraft)
        @test flags == [false, false]
    end

    @testset "a link beyond its terminal's limits is refused" begin
        spacecraft = [terminal_sat(1), terminal_sat(2)]
        model = InterlinkModel()
        register_link!(model, spacecraft, (1, 1), (2, 1), LaserMomentumLink(power_w=500.0))
        err = try
            run_simulation(run_args(spacecraft, model; duration=1.0))
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError && occursin("outside its limits", sprint(showerror, err))
    end

    @testset "checkpoints record the terminal hardware" begin
        spacecraft, model = linked_pair()
        args = run_args(spacecraft, model)
        laser = (module_path=(:SpaceAGORA, :SimulationModel, :Components), name=:LaserTerminal, parameters=())
        hardware = SE._interlink_checkpoint_configuration(args).hardware_terminals
        @test [[terminal.kind for terminal in terminals] for terminals in hardware] == [[laser], [laser]]
        # Every field, floats with their bits: the efficiency sets the link-driven draw.
        fields = only(hardware[1]).fields
        @test keys(fields) == fieldnames(LaserTerminal)
        @test fields.efficiency == (type=(module_path=(:Core,), name=:Float64, parameters=()), value=0.5,
            bits=reinterpret(UInt64, 0.5))
        virtual = [terminal_sat(1; terminals=0, n_terminal=1), terminal_sat(2; terminals=0, n_terminal=1)]
        plain = InterlinkModel()
        register_link!(plain, virtual, (1, 1), (2, 1), LaserMomentumLink())
        @test !haskey(SE._interlink_checkpoint_configuration(run_args(virtual, plain)), :hardware_terminals)
    end
end

end # module
