module InterlinkLaserPhysicsTests

# The terminal-aware link contract (the six-argument `link_available`,
# `interlink_wrench` and its torques, `link_applies_torque`, `terminal_pose`) and
# the laser beam link's opt-in physics: occultation, divergence, pointing and
# mounted recoil.

using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const IM = SM.InterLinkModels
const PLANET = make_no_gram_planet(:earth)
const RADIUS = PLANET.Rp_e + 700e3
const SPEED = sqrt(PLANET.μ / RADIUS)
const ZERO3 = SVector(0.0, 0.0, 0.0)

bits(x) = reinterpret(UInt64, collect(Float64, x))
sample(pos; index=1) = IM.StateSample(SVector{3, Float64}(pos), SVector(0.0, SPEED, 0.0), 100.0, nothing, nothing,
    nothing, index)

# A 100 kg bus on a circular orbit, `y` metres along track, with `terminals`
# laser terminals at `location` and `n_terminal` virtual terminals.
function bus_sat(index; y=0.0, terminals=0, n_terminal=1, location=ZERO3, ang_vel=ZERO3)
    bus = Link(root=true, m=100.0)
    ic = CartesianInitialCondition(SVector(RADIUS, y, 0.0), SVector(0.0, SPEED, 0.0); q=SVector(0.0, 0.0, 0.0, 1.0),
        ang_vel)
    sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, initial_condition=ic, id=index,
        n_terminal=n_terminal)
    for _ in 1:terminals
        add_actuator!(sc, bus, LaserTerminal(boresight=[0.0, 1.0, 0.0], location=location, max_power_w=50.0))
    end
    return sc
end

function run_args(spacecraft, model; policy=SemimajorAxisRatePolicy(), duration=60.0, orientation=false,
                  control=ControlModel((), Float64[]), settings=SimulationSettings(results=false, generate_plots=false))
    return SimulationConfiguration(
        simulation_settings=settings,
        mission_configuration=MissionConfiguration(mission_type=MissionTime, mission_time=duration,
            orientation_sim=orientation, data_rate=1.0),
        environment_model=make_no_gram_environment(planet=PLANET),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]), navigation_model=NavigationModel((), Float64[]),
        control_model=control, initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-8, dt_max_orbit=1.0),
        solver_config=SM.SolverConfig(solver_mode=:tsit5), interlink_model=model,
        scheduling_policy_model=SchedulingPolicyModel(policy))
end

# A kind written with the four-argument hooks only, as kinds were before the
# terminal-aware forms; and the same force given through `interlink_wrench` only.
struct FourArgumentLink <: AbstractInterlink
    range_m::Float64
end
SpaceAGORA.link_range_m(link::FourArgumentLink) = link.range_m
SpaceAGORA.link_available(link::FourArgumentLink, a, b, ctx) = norm(a.pos_ii - b.pos_ii) < link.range_m / 2
SpaceAGORA.link_force(link::FourArgumentLink, a, b, ctx) = (0.3 * (a.pos_ii - b.pos_ii) / 7.0, 0.3 * (b.pos_ii - a.pos_ii) / 7.0)

struct WrenchOnlyLink <: AbstractInterlink
    range_m::Float64
end
SpaceAGORA.link_range_m(link::WrenchOnlyLink) = link.range_m
SpaceAGORA.interlink_wrench(link::WrenchOnlyLink, terminal_a, terminal_b, a, b, ctx) =
    (0.3 * (a.pos_ii - b.pos_ii) / 7.0, ZERO3, 0.3 * (b.pos_ii - a.pos_ii) / 7.0, ZERO3)

# A kind that records the terminals its six-argument hooks receive and turns its
# two endpoints the opposite ways about body z, with no force.
const RECEIVED = Any[]
struct TurningLink <: AbstractInterlink
    torque_n_m::Float64
    declared::Bool
end
SpaceAGORA.link_range_m(::TurningLink) = 1e9
SpaceAGORA.link_applies_torque(link::TurningLink) = link.declared
SpaceAGORA.terminal_supports(::LaserTerminal, ::TurningLink) = true
SpaceAGORA.link_available(::TurningLink, terminal_a, terminal_b, a, b, ctx) = (push!(RECEIVED, (terminal_a, terminal_b)); true)
SpaceAGORA.interlink_wrench(link::TurningLink, terminal_a, terminal_b, a, b, ctx) =
    (ZERO3, SVector(0.0, 0.0, link.torque_n_m), ZERO3, SVector(0.0, 0.0, -link.torque_n_m))

# Selects every available link.
struct MixedStationInterrupted <: Exception end
struct EveryLinkPolicy <: AbstractSchedulingPolicy end
SpaceAGORA.candidate_score(::EveryLinkPolicy, candidate, ctx, t) = 1.0

# A beam model of one's own: a constant share of the near-field force.
struct ConstantBeam <: AbstractLaserBeam
    factor::Float64
end
SpaceAGORA.beam_force_factor(beam::ConstantBeam, link, terminal_a, terminal_b, distance_m) = beam.factor

# Holds levels for the mounted actuators of spacecraft 1.
mutable struct HoldLevels <: AbstractControlEffectorModel
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::HoldLevels) = :mounted
SpaceAGORA.actuator_levels(c::HoldLevels, i::Int) = i == 1 ? c.levels : nothing
SpaceAGORA.update_control!(::HoldLevels, u, p, t::Float64, i::Int64) = nothing

const NEAR_FIELD_N = 100.0 * 40.0 / 299_792_458.0   # magnification × power / c, for 40 W at 100
const rot = SM.QuaternionMath.rot

# Spacecraft 1 and, 1 km behind it along track, spacecraft 2, each with one laser
# terminal: boresights and locations as given (spacecraft 1's points back at 2,
# and 2's forward at 1, by default).
function beam_pair(link; boresights=([0.0, -1.0, 0.0], [0.0, 1.0, 0.0]), location=ZERO3, locations=(location, location))
    spacecraft = [bus_sat(1; n_terminal=0), bus_sat(2; y=-1000.0, n_terminal=0)]
    for (sc, boresight, offset) in zip(spacecraft, boresights, locations)
        add_actuator!(sc, sc.root, LaserTerminal(boresight=boresight, location=offset, max_power_w=50.0))
    end
    model = InterlinkModel()
    register_link!(model, spacecraft, (1, 1), (2, 1), link)
    return spacecraft, model
end

# The two ends of a `beam_pair` as the hooks see them, with attitudes `q1` and `q2`.
function beam_ends(spacecraft; q1=SVector(0.0, 0.0, 0.0, 1.0), q2=SVector(0.0, 0.0, 0.0, 1.0))
    ends = map(((1, 0.0, q1), (2, -1000.0, q2))) do (index, y, q)
        state = IM.StateSample(SVector(RADIUS, y, 0.0), SVector(0.0, SPEED, 0.0), 100.0, q, ZERO3, spacecraft[index], index)
        (state, only(mounted_link_terminals(spacecraft[index])).terminal)
    end
    return ends[1][1], ends[2][1], ends[1][2], ends[2][2]
end

# The gve example, in a module of its own.
module GveExample
include(joinpath(@__DIR__, "..", "..", "..", "examples", "gve_sma_interlinks.jl"))
end

function turning_pair(link; ang_vel=ZERO3)
    spacecraft = [bus_sat(1; terminals=1, n_terminal=0, ang_vel), bus_sat(2; y=-1000.0, ang_vel)]
    model = InterlinkModel()
    register_link!(model, spacecraft, (1, 1), (2, 1), link)
    return spacecraft, model
end

@testset "Interlink laser physics" begin
    @testset "the terminal-aware forms call the four-argument hooks" begin
        a, b = sample([RADIUS, 0.0, 0.0]), sample([RADIUS, -1500.0, 0.0]; index=2)
        far = sample([RADIUS, -9000.0, 0.0]; index=2)
        ctx = LinkContext(time_s=3.0, planet=PLANET)
        link = FourArgumentLink(10_000.0)
        for other in (b, far)
            @test link_available(link, nothing, nothing, a, other, ctx) === link_available(link, a, other, ctx)
        end
        @test link_available(link, nothing, nothing, a, b, ctx) && !link_available(link, nothing, nothing, a, far, ctx)
        force_a, torque_a, force_b, torque_b = interlink_wrench(link, nothing, nothing, a, b, ctx)
        expected = link_force(link, a, b, ctx)
        @test bits(force_a) == bits(expected[1]) && bits(force_b) == bits(expected[2])
        @test torque_a === ZERO3 && torque_b === ZERO3
        @test !link_applies_torque(link)
        # The built-in kinds, likewise.
        terminal = LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=50.0)
        for builtin in (LaserMomentumLink(power_w=40.0, magnification=100.0, range_m=200e3),
                        PowerBeamLink(power_w=40.0, range_m=200e3),
                        DataLink(transmit_power_w=1.0, budget=FixedRateBudget(1e6), range_m=200e3))
            for ends in ((nothing, nothing), (terminal, nothing))
                wrench = interlink_wrench(builtin, ends..., a, b, ctx)
                forces = link_force(builtin, a, b, ctx)
                @test bits(wrench[1]) == bits(forces[1]) && bits(wrench[3]) == bits(forces[2])
                @test iszero(wrench[2]) && iszero(wrench[4])
                @test link_available(builtin, ends..., a, b, ctx) === link_available(builtin, a, b, ctx)
            end
            @test !link_applies_torque(builtin)
        end
        # The power beam's occultation, through the six-argument form.
        behind = sample([-RADIUS, 0.0, 0.0]; index=2)
        @test !link_available(PowerBeamLink(power_w=40.0, range_m=3e7), nothing, nothing, a, behind, ctx)
    end

    @testset "the engine and the policies read interlink_wrench" begin
        spacecraft = [bus_sat(1), bus_sat(2; y=-1000.0)]
        models = map((FourArgumentLink(10_000.0), WrenchOnlyLink(10_000.0))) do link
            model = InterlinkModel()
            register_link!(model, spacecraft, (2, 1), (1, 1), link)
            model
        end
        args = run_args(spacecraft, first(models))
        state = SE.build_initial_conditions(args)
        mu = PLANET.μ
        scores = map(models) do model
            scheduler = SchedulingPolicyModel(SemimajorAxisRatePolicy())
            update_availability!(model, spacecraft, state)
            score_candidates!(model, scheduler, state, mu; spacecraft)
            copy(scheduler.matching.score)
        end
        @test bits(scores[1]) == bits(scores[2]) && all(>(0), scores[1])
        for model in models
            model.graph.active .= true
            force = interlink_force_on_spacecraft(model, state, 1; spacecraft)
            @test bits(force) == bits(0.3 * (state.sc[1].pos - state.sc[2].pos) / 7.0)
            @test !last(IM.interlink_wrench_on_spacecraft(model, state, 1; spacecraft))
        end
        # A run of each ends on the same bits, with the link selected.
        finals = map(models) do model
            solution = run_simulation(run_args(spacecraft, model); return_solution=true)
            @test any(sample -> !isempty(sample.active), solution.prob.p.args.interlink_model.history)
            solution.u[end]
        end
        @test bits(finals[1]) == bits(finals[2])
    end

    @testset "the hooks receive the terminal hardware" begin
        spacecraft, model = turning_pair(TurningLink(0.0, false))
        terminal = only(mounted_link_terminals(spacecraft[1])).terminal
        # Outside a run, from the spacecraft models; without them, none.
        empty!(RECEIVED)
        state = SE.build_initial_conditions(run_args(spacecraft, model))
        update_availability!(model, spacecraft, state)
        update_availability!(model, nothing, state)
        @test RECEIVED[1][1] === terminal && RECEIVED[1][2] === nothing
        @test RECEIVED[2] === (nothing, nothing)
        # In a run, as its setup recorded them.
        empty!(RECEIVED)
        solution = run_simulation(run_args(spacecraft, model; policy=EveryLinkPolicy(), duration=5.0); return_solution=true)
        @test !isempty(RECEIVED) && all(ends -> ends[1] === terminal && ends[2] === nothing, RECEIVED)
        @test solution.prob.p.args.interlink_model.edge_terminals == [(terminal, nothing)]
    end

    @testset "torques turn the endpoints by I⁻¹ τ" begin
        torque = 2.0e-3
        spacecraft, model = turning_pair(TurningLink(torque, true))
        duration = 60.0
        solution = run_simulation(run_args(spacecraft, model; policy=EveryLinkPolicy(), duration, orientation=true);
            return_solution=true)
        final = solution.u[end]
        # About a principal axis from rest, ω grows as τ t / I.
        for (index, sign) in ((1, 1.0), (2, -1.0))
            inertia = spacecraft[index].inertia_tensor
            @test final.sc[index].ω[3] ≈ sign * torque * duration / inertia[3, 3] rtol = 1e-9
            @test abs(final.sc[index].ω[1]) < 1e-15 && abs(final.sc[index].ω[2]) < 1e-15
        end
        # No force: the orbits are those of a run without the link (to the
        # integration's tolerance, since the attitude states differ).
        unlinked = run_simulation(run_args(spacecraft, InterlinkModel(); policy=EveryLinkPolicy(), duration,
            orientation=true); return_solution=true).u[end]
        @test final.sc[1].pos ≈ unlinked.sc[1].pos rtol = 1e-12
        # A link whose torques are zero leaves a tumbling spacecraft's attitude
        # equation as it is, bit for bit.
        tumble = SVector(0.01, -0.02, 0.03)
        still, still_model = turning_pair(TurningLink(0.0, true); ang_vel=tumble)
        resting = run_simulation(run_args(still, still_model; policy=EveryLinkPolicy(), duration, orientation=true);
            return_solution=true).u[end]
        free_pair, _ = turning_pair(TurningLink(0.0, true); ang_vel=tumble)
        free = run_simulation(run_args(free_pair, InterlinkModel(); policy=EveryLinkPolicy(), duration, orientation=true);
            return_solution=true).u[end]
        @test bits(resting) == bits(free)
        @test resting.sc[1].ω != tumble
    end

    @testset "a torque needs attitude" begin
        # Declared: refused at setup, before any solve.
        spacecraft, model = turning_pair(TurningLink(1e-3, true))
        declared = run_args(spacecraft, model; policy=EveryLinkPolicy(), duration=5.0)
        err = try
            SE._validate_interlinks!(declared, :tsit5)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("applies torques to its endpoints (link_applies_torque)", sprint(showerror, err))
        @test_throws ArgumentError run_simulation(declared)
        # Not declared: setup lets it through, and the run stops at its first
        # nonzero torque.
        spacecraft, model = turning_pair(TurningLink(1e-3, false))
        undeclared = run_args(spacecraft, model; policy=EveryLinkPolicy(), duration=5.0)
        @test SE._validate_interlinks!(undeclared, :tsit5) === nothing
        err = try
            run_simulation(undeclared)
            nothing
        catch e
            e
        end
        @test err !== nothing && occursin("applied a torque to spacecraft 1 in a run without attitude", sprint(showerror, err))
        # Not declared, with attitude: applied.
        spacecraft, model = turning_pair(TurningLink(1e-3, false))
        final = run_simulation(run_args(spacecraft, model; policy=EveryLinkPolicy(), duration=5.0, orientation=true);
            return_solution=true).u[end]
        @test final.sc[1].ω[3] > 0
    end

    @testset "terminal_pose" begin
        bus = Link(root=true, m=100.0)
        turn = SVector(0.0, 0.0, sin(pi / 4), cos(pi / 4))   # 90° about z
        boom = Link(root=false, m=5.0, r=MVector(0.0, 1.5, 0.2), q=MVector(turn...))
        sc = SpacecraftModel(root=bus, links=[bus, boom], inertia_tensor=bus.inertia, id=1, n_terminal=0,
            initial_condition=CartesianInitialCondition(SVector(RADIUS, 0.0, 0.0), SVector(0.0, SPEED, 0.0)))
        laser = LaserTerminal(boresight=[1.0, 0.0, 0.0], location=[0.25, 0.0, -0.1], max_power_w=10.0)
        add_actuator!(sc, boom, laser)
        receiver = PowerReceiver(normal=[0.0, 0.0, 1.0], location=[0.0, 0.3, 0.0], radius_m=0.5, efficiency=0.4)
        add_actuator!(sc, bus, receiver)
        ctx = SM.ActuatorMounting.actuator_context(boom)
        boresight, location = terminal_pose(sc, laser)
        @test boresight == ctx.link_to_body * laser.boresight
        @test location == ctx.link_to_body * laser.location + ctx.link_r_body
        @test boresight ≈ SVector(0.0, 1.0, 0.0) atol = 1e-15
        # A force at the terminal has the torque of the terminal's own recoil.
        recoil, recoil_torque = actuator_wrench(laser, 7.0, ctx)
        @test bits(cross(location, recoil)) == bits(recoil_torque)
        # A receiver points along its normal, from the root link.
        @test terminal_pose(sc, receiver) == (receiver.normal, receiver.location)
        # Equal terminals at different poses cannot be told apart; at one pose they can.
        add_actuator!(sc, bus, LaserTerminal(boresight=[1.0, 0.0, 0.0], location=[0.25, 0.0, -0.1], max_power_w=10.0))
        @test_throws ArgumentError terminal_pose(sc, laser)
        @test_throws ArgumentError terminal_pose(sc, LaserTerminal(boresight=[0.0, 0.0, 1.0], max_power_w=1.0))
        twin = bus_sat(2; terminals=2)
        @test terminal_pose(twin, mounted_link_terminals(twin)[1].terminal) == (SVector(0.0, 1.0, 0.0), ZERO3)
    end

    @testset "LaserBeamLink without options is LaserMomentumLink" begin
        ctx = LinkContext(time_s=0.0, planet=PLANET)
        plain = LaserMomentumLink(power_w=40.0, magnification=100.0, range_m=200e3)
        beam_link = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3)
        a = sample([RADIUS, 0.0, 0.0])
        for other in (sample([RADIUS, -1500.0, 20.0]; index=2), sample([RADIUS, 0.0, 0.0]; index=2),
                      sample([RADIUS, -3e5, 0.0]; index=2))
            wrench = interlink_wrench(beam_link, nothing, nothing, a, other, ctx)
            forces = link_force(plain, a, other, ctx)
            @test bits(wrench[1]) == bits(forces[1]) && bits(wrench[3]) == bits(forces[2])
            @test wrench[2] === ZERO3 && wrench[4] === ZERO3
            @test bits(link_force(beam_link, a, other, ctx)[1]) == bits(forces[1])
            @test link_available(beam_link, nothing, nothing, a, other, ctx) === link_available(plain, a, other, ctx)
        end
        @test link_terminal_levels(beam_link) == link_terminal_levels(plain)
        @test !link_applies_torque(beam_link) && !link_directed(beam_link)
        terminal = LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=50.0, magnification=80.0, range_m=1e5)
        @test LaserBeamLink(terminal) === LaserBeamLink(power_w=50.0, magnification=80.0, range_m=1e5)
        # A run with each ends on the same bits, with the same schedule.
        runs = map((plain, beam_link)) do link
            spacecraft = [bus_sat(1), bus_sat(2; y=-1000.0)]
            model = InterlinkModel()
            register_link!(model, spacecraft, (1, 1), (2, 1), link)
            solution = run_simulation(run_args(spacecraft, model; duration=120.0); return_solution=true)
            (solution.u[end], [sample.active for sample in solution.prob.p.args.interlink_model.history])
        end
        @test bits(runs[1][1]) == bits(runs[2][1]) && runs[1][2] == runs[2][2]
        @test any(!isempty, runs[1][2])
        # Option checks.
        @test_throws ArgumentError LaserBeamLink(power_w=40.0, magnification=100.0, range_m=1e5, recoil=:sideways)
        @test_throws ArgumentError LaserBeamLink(power_w=40.0, magnification=100.0, range_m=1e5, occultation=(; altitude=1.0))
        @test_throws ArgumentError LaserBeamLink(power_w=40.0, magnification=100.0, range_m=1e5,
            pointing=TerminalPointing(field_of_regard_rad=0.1, error_rad=1e-6))   # an error needs a GaussianBeam
        @test_throws ArgumentError TerminalPointing(field_of_regard_rad=4.0)
        @test_throws ArgumentError GaussianBeam(waist_radius_m=0.0, aperture_radius_m=0.1)
    end

    @testset "occultation" begin
        R = PLANET.Rp_e
        link(occultation) = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=3e7; occultation)
        ctx = LinkContext(planet=PLANET)
        pair(x) = (sample([x, -1e6, 0.0]), sample([x, 1e6, 0.0]; index=2))
        available(l, (a, b)) = link_available(l, nothing, nothing, a, b, ctx)
        through = (sample([R + 7e5, 0.0, 0.0]), sample([-(R + 7e5), 0.0, 0.0]; index=2))
        blocked = link(true)
        @test !available(blocked, through)
        @test !available(blocked, pair(R - 1.0))        # a segment that dips below the surface
        @test available(blocked, pair(R + 1.0))
        @test available(blocked, pair(R + 5e4))
        # The grazing margin raises the sphere.
        margin = link((; min_ray_altitude_m=1e5))
        @test !available(margin, pair(R + 5e4)) && available(margin, pair(R + 1.5e5))
        # The momentum link and the default beam link, documented to ignore the
        # planet, stay available through it; without a planet nothing blocks.
        @test link_available(LaserMomentumLink(power_w=40.0, magnification=100.0, range_m=3e7), through..., ctx)
        @test available(link(false), through)
        @test link_available(blocked, nothing, nothing, through..., LinkContext())
    end

    @testset "divergence" begin
        ctx = LinkContext(planet=PLANET)
        waist, aperture, wavelength = 0.05, 0.08, 1.064e-6
        beam = GaussianBeam(; wavelength_m=wavelength, waist_radius_m=waist, aperture_radius_m=aperture)
        rayleigh = pi * waist^2 / wavelength
        link = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=1e8; beam)
        a = sample([RADIUS, 0.0, 0.0])
        force_at(l, d) = interlink_wrench(l, nothing, nothing, a, sample([RADIUS, -d, 0.0]; index=2), ctx)[1]
        for d in (1e-3 * rayleigh, rayleigh, 1e3 * rayleigh)
            spot = waist * sqrt(1 + (d / rayleigh)^2)
            capture = -expm1(-2 * aperture^2 / spot^2)
            force = force_at(link, d)
            @test norm(force) ≈ capture * NEAR_FIELD_N rtol = 1e-12
            @test force / norm(force) ≈ SVector(0.0, 1.0, 0.0)          # pushes spacecraft 1 away from 2
            @test beam_force_factor(beam, link, nothing, nothing, d) ≈ capture rtol = 1e-14
        end
        # Near the waist the force is the near-field force, less the aperture's
        # loss; far beyond the Rayleigh range it falls as 1/d².
        @test norm(force_at(link, 1e-3 * rayleigh)) ≈ (1 - exp(-2 * (aperture / waist)^2)) * NEAR_FIELD_N rtol = 1e-6
        @test norm(force_at(link, 2e3 * rayleigh)) / norm(force_at(link, 1e3 * rayleigh)) ≈ 0.25 rtol = 1e-5
        # A beam model of one's own.
        own = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=1e8, beam=ConstantBeam(0.25))
        @test norm(force_at(own, 5e3)) ≈ 0.25 * NEAR_FIELD_N rtol = 1e-14
    end

    @testset "pointing" begin
        ctx = LinkContext(planet=PLANET)
        narrow = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3,
            pointing=TerminalPointing(field_of_regard_rad=0.1))
        spacecraft, _ = beam_pair(narrow)
        a, b, terminal_a, terminal_b = beam_ends(spacecraft)
        @test link_available(narrow, terminal_a, terminal_b, a, b, ctx)              # facing each other
        @test_throws ArgumentError link_available(narrow, nothing, terminal_b, a, b, ctx)
        # Turned 90° about z, spacecraft 1 looks away: the link cannot close; a
        # field of regard just past 90° lets it.
        turn = SVector(0.0, 0.0, sin(pi / 4), cos(pi / 4))
        a_turned, b_same, _, _ = beam_ends(spacecraft; q1=turn)
        @test !link_available(narrow, terminal_a, terminal_b, a_turned, b_same, ctx)
        wide = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3,
            pointing=TerminalPointing(field_of_regard_rad=pi / 2 + 0.01))
        @test link_available(wide, terminal_a, terminal_b, a_turned, b_same, ctx)
        # The line of sight is read in the body frame: a boresight along the turned
        # line of sight closes the link only in that attitude.
        sight_body = rot(turn) * SVector(0.0, -1.0, 0.0)
        turned_pair, _ = beam_pair(narrow; boresights=(sight_body, [0.0, 1.0, 0.0]))
        p, q, tp, tq = beam_ends(turned_pair; q1=turn)
        @test link_available(narrow, tp, tq, p, q, ctx)
        p0, q0, _, _ = beam_ends(turned_pair)
        @test !link_available(narrow, tp, tq, p0, q0, ctx)
        # Each end must see the other: spacecraft 2 looking away closes it.
        away, _ = beam_pair(narrow; boresights=([0.0, -1.0, 0.0], [1.0, 0.0, 0.0]))
        @test !link_available(narrow, beam_ends(away)[3], beam_ends(away)[4], beam_ends(away)[1], beam_ends(away)[2], ctx)
        # A static pointing error scales the force by the far-field Gaussian loss.
        beam = GaussianBeam(waist_radius_m=0.05, aperture_radius_m=0.08)
        divergence = beam.wavelength_m / (pi * beam.waist_radius_m)
        errored = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3; beam,
            pointing=TerminalPointing(field_of_regard_rad=0.1, error_rad=0.5 * divergence))
        ideal = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3; beam,
            pointing=TerminalPointing(field_of_regard_rad=0.1))
        @test norm(interlink_wrench(errored, terminal_a, terminal_b, a, b, ctx)[1]) ≈
            exp(-0.5) * norm(interlink_wrench(ideal, terminal_a, terminal_b, a, b, ctx)[1]) rtol = 1e-14
        # It needs terminal hardware at both ends, and attitude.
        virtual = [bus_sat(1; terminals=1, n_terminal=0), bus_sat(2; y=-1000.0)]
        @test_throws ArgumentError register_link!(InterlinkModel(), virtual, (1, 1), (2, 1), narrow)
        err = try
            run_simulation(run_args(beam_pair(narrow)...; duration=5.0)); nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("orientation_sim", sprint(showerror, err))
        # In a run with attitude, facing terminals link and turned-away ones do not.
        actives(link; kwargs...) = begin
            solution = run_simulation(run_args(beam_pair(link; kwargs...)...; duration=5.0, orientation=true);
                return_solution=true)
            [!isempty(sample.active) for sample in solution.prob.p.args.interlink_model.history]
        end
        @test all(actives(narrow))
        @test !any(actives(narrow; boresights=([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])))
    end

    @testset "mounted recoil" begin
        ctx = LinkContext(planet=PLANET)
        location = SVector(0.3, 0.0, 0.2)
        mounted = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3, recoil=:mounted)
        central = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3)
        @test link_applies_torque(mounted)
        spacecraft, _ = beam_pair(mounted; location)
        a, b, terminal_a, terminal_b = beam_ends(spacecraft)
        force_a, torque_a, force_b, torque_b = interlink_wrench(mounted, terminal_a, terminal_b, a, b, ctx)
        central_wrench = interlink_wrench(central, terminal_a, terminal_b, a, b, ctx)
        @test bits(force_a) == bits(central_wrench[1]) && bits(force_b) == bits(central_wrench[3])
        # The force at the terminal's location: its lever arm's torque, which is the
        # terminal's own recoil for the same force (its beam points at the other end).
        _, arm = terminal_pose(spacecraft[1], terminal_a)
        @test bits(torque_a) == bits(cross(arm, force_a))
        _, arm_b = terminal_pose(spacecraft[2], terminal_b)
        @test bits(torque_b) == bits(cross(arm_b, force_b))        # each end's own force at its own terminal
        @test !iszero(torque_b)
        recoil_force, recoil_torque = actuator_wrench(terminal_a, 100.0 * 40.0, SM.ActuatorMounting.actuator_context(spacecraft[1].root))
        @test recoil_force ≈ force_a rtol = 1e-14
        @test recoil_torque ≈ torque_a rtol = 1e-14
        # In attitude, the force is turned into the body frame first.
        turn = SVector(0.0, 0.0, sin(pi / 8), cos(pi / 8))
        a_turned, _, _, _ = beam_ends(spacecraft; q1=turn)
        turned_wrench = interlink_wrench(mounted, terminal_a, terminal_b, a_turned, b, ctx)
        @test turned_wrench[2] ≈ cross(arm, rot(turn) * turned_wrench[1]) rtol = 1e-14
        # It needs terminal hardware at both ends.
        virtual = [bus_sat(1; terminals=1, n_terminal=0), bus_sat(2; y=-1000.0)]
        @test_throws ArgumentError register_link!(InterlinkModel(), virtual, (1, 1), (2, 1), mounted)
        # A two-spacecraft run: the link's impulses cancel, and each body spins up
        # by I⁻¹ τ (trapezoidal over 10 s, the line of sight turning slowly).
        duration = 10.0
        solution = run_simulation(run_args(beam_pair(mounted; location)...; duration, orientation=true); return_solution=true)
        first_state, final = solution.u[1], solution.u[end]
        @test all(!isempty(sample.active) for sample in solution.prob.p.args.interlink_model.history)
        @test norm(100.0 * final.sc[1].interlink_dv + 100.0 * final.sc[2].interlink_dv) <= 1e-12 * norm(final.sc[1].interlink_dv)
        torque_at(state, index) = begin
            run_sc = solution.prob.p.args.dynamics_model.spacecraft
            ends = [IM.StateSample(SVector{3}(state.sc[i].pos), SVector{3}(state.sc[i].vel), state.sc[i].mass,
                SVector{4}(state.sc[i].q), SVector{3}(state.sc[i].ω), run_sc[i], i) for i in 1:2]
            terminals = [only(mounted_link_terminals(run_sc[i])).terminal for i in 1:2]
            interlink_wrench(mounted, terminals..., ends..., ctx)[2index]
        end
        for index in 1:2
            inertia = spacecraft[index].inertia_tensor
            expected = inertia \ ((torque_at(first_state, index) + torque_at(final, index)) / 2 * duration)
            @test final.sc[index].ω ≈ expected rtol = 1e-4
            @test norm(final.sc[index].ω) > 0
        end
        # D-I2-1: while the link drives the terminal, a controller's level for it
        # applies no recoil.
        held = run_simulation(run_args(beam_pair(mounted; location)...; duration=5.0, orientation=true,
            control=ControlModel(control_effectors=(HoldLevels([30.0]),), control_rates=[1.0])); return_solution=true)
        @test SE._mounted_control_effect(HoldLevels([30.0]), held.u[end].sc[1], held.prob.p, 1, 5.0)[1] == zeros(3)
    end

    @testset "mounted beam uses the terminal-to-terminal ray" begin
        mounted = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3, recoil=:mounted)
        spacecraft, _ = beam_pair(mounted; locations=(SVector(1.0, 0.0, 0.0), SVector(-1.0, 0.0, 0.0)))
        for angle in (0.0, pi / 4)
            turn = SVector(0.0, 0.0, sin(angle / 2), cos(angle / 2))
            a, b, ta, tb = beam_ends(spacecraft; q1=turn)
            arm_a = terminal_pose(a.spacecraft, ta)[2]
            arm_b = terminal_pose(b.spacecraft, tb)[2]
            offset_a = transpose(rot(a.q_ib)) * arm_a
            offset_b = transpose(rot(b.q_ib)) * arm_b
            ray = (a.pos_ii - b.pos_ii) + (offset_a - offset_b)
            fa, torque_a, fb, torque_b = interlink_wrench(mounted, ta, tb, a, b, LinkContext())
            @test fa ≈ NEAR_FIELD_N * normalize(ray) rtol = 1e-13
            @test fa + fb == ZERO3
            @test torque_a ≈ cross(arm_a, rot(a.q_ib) * fa) rtol = 1e-13
            @test torque_b ≈ cross(arm_b, rot(b.q_ib) * fb) rtol = 1e-13
            # About spacecraft b: orbital moment plus both body torques must
            # cancel for equal/opposite forces on the same physical ray.
            net_moment = cross(a.pos_ii - b.pos_ii, fa) +
                transpose(rot(a.q_ib)) * torque_a + transpose(rot(b.q_ib)) * torque_b
            @test norm(net_moment) <= 1e-12 * NEAR_FIELD_N
        end
        a, b, ta, tb = beam_ends(spacecraft)
        # A 2 m lateral offset over 1 km lies outside this 1 mrad cone,
        # though the spacecraft-center line lies exactly on both boresights.
        narrow = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3,
            pointing=TerminalPointing(field_of_regard_rad=0.001))
        @test !link_available(narrow, ta, tb, a, b, LinkContext())
        short = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=1000.001, recoil=:mounted)
        @test !link_available(short, ta, tb, a, b, LinkContext())
        beam = GaussianBeam(waist_radius_m=0.05, aperture_radius_m=0.08)
        diverging = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3, recoil=:mounted, beam=beam)
        fa = interlink_wrench(diverging, ta, tb, a, b, LinkContext())[1]
        distance = hypot(1000.0, 2.0)
        @test norm(fa) ≈ NEAR_FIELD_N * beam_capture_fraction(beam, beam.aperture_radius_m, distance) rtol = 1e-13
        # Mounted terminals clear a sphere that blocks the center-to-center ray.
        occulted = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3,
            recoil=:mounted, occultation=true)
        clear_pair, _ = beam_pair(occulted; location=SVector(2.0, 0.0, 0.0))
        ca, cb, cta, ctb = beam_ends(clear_pair)
        @test link_available(occulted, cta, ctb, ca, cb, LinkContext(planet=(Rp_e=RADIUS + 1.0,)))
    end

    @testset "checkpoints keep the options" begin
        mktempdir() do directory
            segments(resume) = SimulationSettings(results=false, generate_plots=false, checkpoint_enabled=true,
                checkpoint_interval_s=10.0, checkpoint_directory=directory, resume_from_checkpoint=resume)
            options = (beam=GaussianBeam(waist_radius_m=0.05, aperture_radius_m=0.08), occultation=true)
            pair(; changes...) = begin
                spacecraft = [bus_sat(1), bus_sat(2; y=-1000.0)]
                model = InterlinkModel()
                register_link!(model, spacecraft, (1, 1), (2, 1),
                    LaserBeamLink(; power_w=40.0, magnification=100.0, range_m=200e3, options..., changes...))
                (spacecraft, model)
            end
            partial = run_simulation(run_args(pair()...; duration=30.0, settings=segments(false)); return_solution=true)
            for change in ((beam=GaussianBeam(waist_radius_m=0.05, aperture_radius_m=0.09),),
                           (beam=ConstantBeam(0.5),), (occultation=(; min_ray_altitude_m=1e3),), (occultation=false,))
                @test_throws ArgumentError run_simulation(run_args(pair(; change...)...; duration=60.0, settings=segments(true)))
            end
            resumed = run_simulation(run_args(pair()...; duration=60.0, settings=segments(true)); return_solution=true)
            history = resumed.prob.p.args.interlink_model.history
            @test first(history).time == 0.0 && last(history).time == 60.0
            @test count(sample -> sample.time == 30.0, history) == 1
            prefix(h) = [(s.time, s.active, s.available, s.scores, s.terms) for s in h if s.time < 30.0]
            @test prefix(history) == prefix(partial.prob.p.args.interlink_model.history)
        end
    end

    @testset "checkpoints preserve terminal mounting geometry" begin
        function mounted_pair(; hinged=false, offset=ZERO3, attitude=SVector(0.0, 0.0, 0.0, 1.0),
                              axis=SVector(0.0, 0.0, 1.0), rate=1.0)
            link = LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3, recoil=:mounted)
            spacecraft, model = beam_pair(link)
            for (index, sc) in enumerate(spacecraft)
                mount = Link(r=index == 1 ? offset : ZERO3,
                    q=index == 1 ? attitude : SVector(0.0, 0.0, 0.0, 1.0))
                push!(sc.links, mount)
                terminal = pop!(sc.root.laser_terminals)
                add_actuator!(sc, mount, terminal)
                hinged && add_actuator!(sc, mount, PanelHinge(; axis, max_rate_rad_s=rate))
            end
            return spacecraft, model
        end
        for hinged in (false, true)
            mktempdir() do directory
                settings(resume) = SimulationSettings(results=false, generate_plots=false, checkpoint_enabled=true,
                    checkpoint_interval_s=1.0, checkpoint_directory=directory, resume_from_checkpoint=resume)
                control() = hinged ? ControlModel((HoldLevels([0.0, 0.3]),), [1.0]) : ControlModel((), Float64[])
                args(; resume=false, duration=3.0, changes...) = run_args(mounted_pair(; hinged, changes...)...;
                    duration, orientation=true, control=control(), settings=settings(resume))
                partial = run_simulation(args(); return_solution=true)
                if hinged
                    saved_sc = partial.prob.p.args.dynamics_model.spacecraft[1]
                    @test only(saved_sc.links[2].hinges).angle_rad ≈ 0.3
                end
                @test_throws ArgumentError run_simulation(args(; resume=true, duration=3.0, offset=SVector(1.0, 0.0, 0.0)))
                if hinged
                    @test_throws ArgumentError run_simulation(args(; resume=true, duration=3.0, axis=SVector(0.0, 1.0, 0.0)))
                    @test_throws ArgumentError run_simulation(args(; resume=true, duration=3.0, rate=0.5))
                else
                    @test_throws ArgumentError run_simulation(args(; resume=true, duration=3.0,
                        attitude=SVector(0.0, 0.0, sin(0.1), cos(0.1))))
                end
                resumed = run_simulation(args(; resume=true, duration=5.0); return_solution=true)
                history = resumed.prob.p.args.interlink_model.history
                @test first(history).time == 0.0
                @test last(history).time == 5.0
                @test count(sample -> sample.time == 3.0, history) == 1
                prefix(h) = [(s.time, s.active, s.available, s.scores, s.terms) for s in h if s.time < 3.0]
                @test prefix(history) == prefix(partial.prob.p.args.interlink_model.history)
                @test last(resumed.t) == 5.0
                if hinged
                    saved_sc = resumed.prob.p.args.dynamics_model.spacecraft[1]
                    @test only(saved_sc.links[2].hinges).angle_rad ≈ 0.3
                end
            end
        end
        @testset "mounted laser and station data resume together" begin
            function mixed_args(directory; resume=false, offset=SVector(1.0, 0.0, 0.0), mask=10.0)
                spacecraft, model = mounted_pair(; offset)
                radio() = RadioTerminal(boresight=[0.0, 0.0, 1.0], max_power_w=5.0, gain_dbi=0.0, efficiency=0.5)
                # Put the station under spacecraft 1 at the initial epoch.
                env = make_no_gram_environment(planet=PLANET)
                epoch = SM.EphemeridesModels.ephemerides_epoch_seconds(env.ephemerides_model, InitialTime())
                L0 = SM.EphemeridesModels.ephemerides_frame_lpi(env.ephemerides_model, PLANET, epoch)
                direction = normalize(L0 * SVector(RADIUS, 0.0, 0.0))
                st = GroundStation(name="mixed", lat_deg=rad2deg(asin(direction[3])),
                    lon_deg=rad2deg(atan(direction[2], direction[1])), min_elevation_deg=mask,
                    antennas=[radio()])
                add_ground_station!(model, st)
                sc = spacecraft[1]
                radio_device = radio()
                add_actuator!(sc, sc.links[2], radio_device)
                # The mounted laser remains terminal 1; resolve the radio by identity.
                register_link!(model, spacecraft, (1, radio_device), (-1, 1),
                    DataLink(transmit_power_w=1.0, budget=FixedRateBudget(1e3), range_m=3e6))
                add_data!(sc, sc.root, DataStore(capacity_bits=1e6, initial_bits=1e5))
                settings = SimulationSettings(results=false, generate_plots=false,
                    checkpoint_enabled=true, checkpoint_interval_s=1.0,
                    checkpoint_directory=directory, resume_from_checkpoint=resume)
                return run_args(spacecraft, model; policy=EveryLinkPolicy(), duration=3.0,
                    orientation=true, settings)
            end
            mktempdir() do reference_dir
                whole = run_simulation(mixed_args(reference_dir); return_solution=true)
                mktempdir() do split_dir
                    args = mixed_args(split_dir)
                    stop = SM.SimulationCallbacks.PresetTimeCallback([1.1],
                        i -> throw(MixedStationInterrupted()); save_positions=(false, false))
                    @test_throws MixedStationInterrupted run_simulation(args; extra_callbacks=(stop,))
                    saved = SE._load_checkpoint(args)
                    @test saved.t == 1.0
                    @test saved.u.gs[1].received_bits > 0.0
                    @test_throws ArgumentError run_simulation(mixed_args(split_dir; resume=true, mask=20.0))
                    @test_throws ArgumentError run_simulation(mixed_args(split_dir; resume=true,
                        offset=SVector(2.0, 0.0, 0.0)))
                    resumed = run_simulation(mixed_args(split_dir; resume=true); return_solution=true)
                    @test resumed.t == whole.t
                    @test all(isequal(collect(a), collect(b)) for (a,b) in zip(resumed.u, whole.u))
                    final = resumed.u[end]
                    @test final.gs[1].received_bits ≈ 3e3 rtol=1e-10
                    @test final.gs[1].received_bits == final.sc[1].data_sent_bits
                    @test norm(final.sc[1].ω) > 0.0
                    @test any(length(h.active) == 2 for h in resumed.prob.p.args.interlink_model.history)
                end
            end
        end

        root_args = run_args(beam_pair(LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3, recoil=:mounted))...)
        original = SE._interlink_checkpoint_configuration(root_args)
        root_args.dynamics_model.spacecraft[1].root.q .= SVector(0.0, 0.0, sin(0.2), cos(0.2))
        @test isequal(SE._interlink_checkpoint_configuration(root_args), original)
        root_args.dynamics_model.spacecraft[1].root.r[1] += 1.0
        @test !isequal(SE._interlink_checkpoint_configuration(root_args), original)
        panel_args = run_args(mounted_pair()...; control=ControlModel(
            (SolarPanelAngleOfAttackControlModel(controlled_panel_links=(2,)),), [1.0]))
        @test_throws ArgumentError SE._interlink_checkpoint_configuration(panel_args)
        # Legacy/default center-to-center beams keep the pre-geometry payload.
        default_args = run_args(beam_pair(LaserBeamLink(power_w=40.0, magnification=100.0, range_m=200e3))...)
        @test !hasproperty(SE._interlink_checkpoint_configuration(default_args), :terminal_mounts)
    end

    @testset "the gve example with a beam" begin
        example = GveExample
        beam = GaussianBeam(waist_radius_m=0.05, aperture_radius_m=0.08)
        plain_args = example.gve_sma_configuration(duration=120.0)
        beam_args = example.gve_sma_configuration(duration=120.0, beam=beam)
        @test all(link -> link isa LaserMomentumLink, plain_args.interlink_model.graph.edge_properties)
        @test all(link -> link isa LaserBeamLink{GaussianBeam} && link.occultation,
            beam_args.interlink_model.graph.edge_properties)
        gains = map((plain_args, beam_args)) do args
            final = run_simulation(args; return_solution=true).u[end]
            abs(final.sc[1].interlink_delta_sma)
        end
        @test 0 < gains[2] < gains[1]          # tens of kilometres apart, the beam delivers less
    end
end

end # module InterlinkLaserPhysicsTests
