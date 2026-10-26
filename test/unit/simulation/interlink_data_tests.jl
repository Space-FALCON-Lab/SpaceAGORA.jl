module InterlinkDataTests

# Data links: the fixed and Shannon budgets and their caps, the radio terminal,
# direction and hardware, occultation, the bits two spacecraft exchange, the
# radio's draw, contact windows, shedding at the sender, scheduling, the states
# and columns, and checkpoints.

using Test
using LinearAlgebra
using StaticArrays
using TOML
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PLANET = make_no_gram_planet(:earth)
const RADIUS = PLANET.Rp_e + 700e3
const SHANNON = ShannonBudget(bandwidth_hz=1e6, frequency_hz=2.2e9)

radio(; kwargs...) = RadioTerminal(; boresight=[0.0, 1.0, 0.0], max_power_w=5.0, gain_dbi=6.0, efficiency=0.5, kwargs...)

# A 100 kg bus `arc_m` along a circular orbit of radius `radius`, with a radio
# and a battery.
function radio_sat(index, arc_m; radius=RADIUS, capacity_wh=10.0, soc=0.5, battery=true, terminal=radio(), n_terminal=0)
    bus = Link(root=true, m=100.0)
    θ = arc_m / RADIUS
    speed = sqrt(PLANET.μ / radius)
    ic = CartesianInitialCondition(SVector(radius * cos(θ), radius * sin(θ), 0.0),
        SVector(-speed * sin(θ), speed * cos(θ), 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, initial_condition=ic, id=index,
        n_terminal=n_terminal)
    terminal === nothing || add_actuator!(sc, bus, terminal)
    battery && add_power!(sc, bus, Battery(capacity_wh=capacity_wh, initial_soc=soc))
    return sc
end

const LINK = DataLink(transmit_power_w=1.0, budget=SHANNON, range_m=2000e3)

function data_pair(; arc_m=100e3, link=LINK, sender=(;), receiver=(;))
    spacecraft = [radio_sat(1, 0.0; sender...), radio_sat(2, arc_m; receiver...)]
    model = InterlinkModel()
    register_link!(model, spacecraft, (1, 1), (2, 1), link)
    return spacecraft, model
end

function run_args(spacecraft, model; duration=60.0, policy=DataThroughputPolicy(),
                  settings=SimulationSettings(results=false, generate_plots=false))
    return SimulationConfiguration(
        simulation_settings=settings,
        mission_configuration=MissionConfiguration(mission_type=MissionTime, mission_time=duration,
            orientation_sim=false, data_rate=1.0),
        environment_model=make_no_gram_environment(planet=PLANET),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]), navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]), initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-8, dt_max_orbit=1.0),
        solver_config=SM.SolverConfig(solver_mode=:tsit5), interlink_model=model,
        scheduling_policy_model=policy)
end

sample(pos) = SM.InterLinkModels.link_endpoint((pos=pos,), nothing, 1)

# The closed-form Shannon rate of `budget` at `d`, with the default radios.
function shannon(budget, d; power_w=1.0, gain_dbi=6.0, temperature_k=500.0)
    λ = 299_792_458.0 / budget.frequency_hz
    received = power_w * 10^(2gain_dbi / 10) * (λ / (4π * d))^2
    noise = 1.380649e-23 * temperature_k * budget.bandwidth_hz
    return min(budget.implementation_efficiency * budget.bandwidth_hz * log2(1 + received / noise), budget.max_rate_bps)
end

chord(arc_m) = 2RADIUS * sin(arc_m / RADIUS / 2)
integral(x, y) = sum((x[k + 1] - x[k]) * (y[k] + y[k + 1]) / 2 for k in 1:length(x) - 1)

@testset "Data links" begin
    @testset "budgets" begin
        a = sample([RADIUS, 0.0, 0.0])
        fixed = DataLink(transmit_power_w=0.0, budget=FixedRateBudget(rate_bps=2e6), range_m=1e7)
        for d in (1e3, 1e5, 5e6)
            @test link_data_rate_bps(fixed, nothing, nothing, a, sample([RADIUS, d, 0.0]), LinkContext()) == 2e6
        end
        tx, rx = radio(), radio()
        for d in (10e3, 100e3, 1000e3, 5000e3)
            b = sample([RADIUS, d, 0.0])
            @test link_data_rate_bps(LINK, tx, rx, a, b, LinkContext()) ≈ shannon(SHANNON, d) rtol=1e-12
        end
        @test 5.6e6 < shannon(SHANNON, 10e3) < 5.8e6 && 1.6e5 < shannon(SHANNON, 1000e3) < 1.8e5
        # The receiver's noise temperature and both gains enter.
        b = sample([RADIUS, 100e3, 0.0])
        @test link_data_rate_bps(LINK, tx, radio(noise_temperature_k=1000.0), a, b, LinkContext()) ≈
            shannon(SHANNON, 100e3; temperature_k=1000.0) rtol=1e-12
        @test link_data_rate_bps(LINK, radio(gain_dbi=0.0), radio(gain_dbi=12.0), a, b, LinkContext()) ≈
            shannon(SHANNON, 100e3) rtol=1e-12
        efficient = ShannonBudget(bandwidth_hz=2e6, frequency_hz=8.4e9, implementation_efficiency=0.8)
        @test link_data_rate_bps(DataLink(transmit_power_w=1.0, budget=efficient, range_m=1e7), tx, rx, a, b, LinkContext()) ≈
            shannon(efficient, 100e3) rtol=1e-12
        # The caps: the budget's, the transmitter's and the receiver's.
        capped = DataLink(transmit_power_w=1.0, budget=ShannonBudget(bandwidth_hz=1e6, frequency_hz=2.2e9, max_rate_bps=1e6),
            range_m=1e7)
        @test link_data_rate_bps(capped, tx, rx, a, b, LinkContext()) == 1e6
        @test link_data_rate_bps(LINK, radio(max_rate_bps=2e5), rx, a, b, LinkContext()) == 2e5
        @test link_data_rate_bps(LINK, tx, radio(max_rate_bps=3e5), a, b, LinkContext()) == 3e5
        @test link_data_rate_bps(fixed, radio(max_rate_bps=1e5), nothing, a, b, LinkContext()) == 1e5
        # A Shannon budget reads both radios; other kinds carry no data.
        @test_throws ArgumentError link_data_rate_bps(LINK, nothing, rx, a, b, LinkContext())
        @test link_data_rate_bps(LaserMomentumLink(), tx, rx, a, b, LinkContext()) == 0.0
        @test link_directed(LINK) && link_terminal_levels(LINK) == (1.0, 0.0)
        @test link_force(LINK, a, b, LinkContext()) == (zeros(3), zeros(3))
        @test budget_rate_bps(FixedRateBudget(3.0), LINK, nothing, nothing, 1.0) == 3.0
        for bad in ((rate_bps=0.0,), (rate_bps=Inf,))
            @test_throws ArgumentError FixedRateBudget(; bad...)
        end
        for bad in ((bandwidth_hz=0.0,), (frequency_hz=-1.0,), (implementation_efficiency=1.5,), (max_rate_bps=0.0,))
            @test_throws ArgumentError ShannonBudget(; merge((bandwidth_hz=1e6, frequency_hz=1e9), bad)...)
        end
        for bad in ((transmit_power_w=-1.0,), (range_m=0.0,))
            @test_throws ArgumentError DataLink(; merge((transmit_power_w=1.0, budget=SHANNON, range_m=1.0), bad)...)
        end
        for bad in ((max_power_w=0.0,), (gain_dbi=Inf,), (efficiency=0.0,), (noise_temperature_k=0.0,), (max_rate_bps=0.0,),
                    (boresight=[0.0, 0.0, 0.0],))
            @test_throws ArgumentError radio(; bad...)
        end
    end

    @testset "terminals and direction" begin
        spacecraft, model = data_pair()
        @test only(model.graph.edges) |> e -> map(v -> model.graph.vertices[v], e) == ((1, 1), (2, 1))
        # A radio transmits or receives, so the reverse link registers too.
        @test register_link!(model, spacecraft, (2, 1), (1, 1), LINK) == ((2, 1), (1, 1))
        terminal = SM.mounted_link_terminals(spacecraft[1])[1].terminal
        @test link_terminal(terminal) && terminal_supports(terminal, LINK, :transmitter) &&
            terminal_supports(terminal, LINK, :receiver)
        @test !terminal_supports(terminal, LaserMomentumLink()) && !terminal_supports(terminal, PowerBeamLink(power_w=1.0, range_m=1.0))
        @test terminal_shed_tier(terminal) == 1
        @test terminal_shed_tier(LaserTerminal(boresight=[1.0, 0.0, 0.0], max_power_w=1.0)) == 1
        @test terminal_shed_tier(PowerReceiver(normal=[1.0, 0.0, 0.0], radius_m=1.0, efficiency=0.5)) == 0
        # Laser terminals and power receivers carry no data link.
        laser = radio_sat(1, 0.0; terminal=LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=5.0))
        receiver = radio_sat(2, 100e3; terminal=PowerReceiver(normal=[0.0, -1.0, 0.0], radius_m=0.2, efficiency=0.5))
        @test_throws ArgumentError register_link!(InterlinkModel(), [laser, spacecraft[2]], (1, 1), (2, 1), LINK)
        @test_throws ArgumentError register_link!(InterlinkModel(), [spacecraft[1], receiver], (1, 1), (2, 1), LINK)
        # A Shannon budget needs radios at both ends; a fixed rate runs between virtual terminals.
        virtual = [radio_sat(1, 0.0; terminal=nothing, n_terminal=1), radio_sat(2, 100e3; terminal=nothing, n_terminal=1)]
        err = try
            register_link!(InterlinkModel(), [virtual[1], spacecraft[2]], (1, 1), (2, 1), LINK)
        catch caught
            caught
        end
        @test err isa ArgumentError && occursin("virtual", sprint(showerror, err))
        fixed = DataLink(transmit_power_w=0.0, budget=FixedRateBudget(1e6), range_m=1e6)
        @test register_link!(InterlinkModel(), virtual, (1, 1), (2, 1), fixed) == ((1, 1), (2, 1))
        # A transmit power above the radio's limit is refused at setup.
        strong, strong_model = data_pair(link=DataLink(transmit_power_w=6.0, budget=SHANNON, range_m=2000e3))
        @test_throws ArgumentError run_simulation(run_args(strong, strong_model; duration=1.0); return_results=true)
    end

    @testset "occultation" begin
        ctx = LinkContext(planet=PLANET)
        @test link_available(LINK, sample([RADIUS, 0.0, 0.0]), sample([RADIUS, 1e5, 0.0]), ctx)
        far = DataLink(transmit_power_w=1.0, budget=SHANNON, range_m=1e8)
        @test !link_available(far, sample([RADIUS, 0.0, 0.0]), sample([-RADIUS, 0.0, 0.0]), ctx)   # through the Earth
        @test link_available(far, sample([RADIUS, 0.0, 0.0]), sample([-RADIUS, 0.0, 0.0]), LinkContext())
        @test !link_available(LINK, sample([RADIUS, 0.0, 0.0]), sample([RADIUS, 3e6, 0.0]), ctx)   # out of range
    end

    @testset "two spacecraft exchange bits" begin
        spacecraft, model = data_pair()
        t = run_simulation(run_args(spacecraft, model); return_results=true).table
        rate = shannon(SHANNON, chord(100e3))
        @test all(≈(rate; rtol=1e-6), t.sc1_data_rate_tx_bps) && all(≈(rate; rtol=1e-6), t.sc2_data_rate_rx_bps)
        @test all(==(0.0), t.sc1_data_rate_rx_bps) && all(==(0.0), t.sc2_data_rate_tx_bps)
        @test t.sc2_data_received_bits[end] ≈ rate * 60.0 rtol=1e-6
        @test t.sc1_data_sent_bits ≈ t.sc2_data_received_bits rtol=1e-12
        @test all(==(0.0), t.sc1_data_received_bits) && all(==(0.0), t.sc2_data_sent_bits)
        @test t.sc2_data_received_bits[end] ≈ integral(t.time, t.sc2_data_rate_rx_bps) rtol=1e-9
        # The transmitting radio draws its power over its efficiency, the receiving one nothing.
        @test all(≈(2.0), t.sc1_power_load_w) && all(==(0.0), t.sc2_power_load_w)
        @test t.sc1_battery_wh[1] - t.sc1_battery_wh[end] ≈ 2.0 * 60.0 / 3600.0 rtol=1e-9
        @test t.sc1_battery_wh[end] - t.sc1_battery_wh[1] ≈ integral(t.time, t.sc1_power_net_w) / 3600.0 rtol=1e-8
        @test all(==(1.0), t.sc1_actuator_level_1) && all(==(0.0), t.sc2_actuator_level_1)
    end

    @testset "the counters integrate a varying rate" begin
        # The receiver 50 km higher drifts back: the separation grows, the rate falls.
        spacecraft, model = data_pair(arc_m=0.0, receiver=(radius=RADIUS + 50e3,))
        t = run_simulation(run_args(spacecraft, model; duration=600.0); return_results=true).table
        @test t.sc2_data_rate_rx_bps[end] < 0.9 * t.sc2_data_rate_rx_bps[1]
        @test t.sc2_data_received_bits[end] ≈ integral(t.time, t.sc2_data_rate_rx_bps) rtol=1e-6
        @test t.sc1_data_sent_bits[end] ≈ t.sc2_data_received_bits[end] rtol=1e-12
    end

    @testset "contact windows" begin
        key = ((1, 1), (2, 1))
        history = [InterLinkHistory(0.0, InterLinkKey[], InterLinkKey[]), InterLinkHistory(1.0, [key], [key]),
                   InterLinkHistory(2.0, InterLinkKey[], [key]), InterLinkHistory(3.0, InterLinkKey[], InterLinkKey[]),
                   InterLinkHistory(4.0, InterLinkKey[], [key])]
        @test contact_windows(history) == [(key=key, start_s=1.0, end_s=3.0, duration_s=2.0),
                                           (key=key, start_s=4.0, end_s=4.0, duration_s=0.0)]
        @test contact_windows(history; active=true) == [(key=key, start_s=1.0, end_s=2.0, duration_s=1.0)]
        @test contact_time_s(history, key) == 2.0 && contact_time_s(history, key; active=true) == 1.0
        @test contact_time_s(history, ((1, 1), (3, 1))) == 0.0
        @test isempty(contact_windows([InterLinkHistory(0.0, [key]), InterLinkHistory(1.0, [key])]))   # nothing available recorded
        @test isempty(contact_windows(InterLinkHistory[]))
        # In a run the receiver drifts out of a 60 km range: one window from the start.
        reach = DataLink(transmit_power_w=1.0, budget=SHANNON, range_m=60e3)
        spacecraft, model = data_pair(arc_m=0.0, link=reach, receiver=(radius=RADIUS + 50e3,))
        solution = run_simulation(run_args(spacecraft, model; duration=600.0); return_solution=true)
        history = solution.prob.p.args.interlink_model.history
        separation = [norm(u.sc[1].pos - u.sc[2].pos) for u in solution.u]
        crossing = solution.t[findfirst(>(60e3), separation)]
        window = only(contact_windows(history))
        @test window.key == key && window.start_s == 0.0
        @test crossing - 1.0 - 1e-9 <= window.end_s <= crossing + 1e-9
        @test contact_windows(history; active=true) == [window]                  # selected while it could close
        @test contact_time_s(history, key) ≈ window.end_s
        # A policy that never selects it: contact, but no active time. The data
        # counters it leaves at zero take part in the step control, so its
        # accepted steps, and the window's end, differ by less than a step.
        spacecraft, model = data_pair(arc_m=0.0, link=reach, receiver=(radius=RADIUS + 50e3,))
        idle = run_simulation(run_args(spacecraft, model; duration=600.0, policy=SemimajorAxisRatePolicy()); return_solution=true)
        idle_history = idle.prob.p.args.interlink_model.history
        @test contact_time_s(idle_history, key) ≈ window.end_s atol=1.0
        @test contact_time_s(idle_history, key; active=true) == 0.0
    end

    @testset "the history records the available links in a data run, or when asked" begin
        key = ((1, 1), (2, 1))
        for (record, recorded) in ((nothing, false), (true, true), (false, false))
            lasers = [radio_sat(i, (i - 1) * 5000.0; terminal=LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=100.0))
                      for i in 1:2]
            model = InterlinkModel(; record_available=record)
            register_link!(model, lasers, (1, 1), (2, 1), LaserMomentumLink(power_w=40.0))
            history = run_simulation(run_args(lasers, model; duration=2.0, policy=SemimajorAxisRatePolicy());
                return_solution=true).prob.p.args.interlink_model.history
            @test all(sample -> sample.available == (recorded ? [key] : InterLinkKey[]), history)
            @test isempty(contact_windows(history)) == !recorded
        end
        spacecraft, model = data_pair()
        @test model.record_available === nothing && SM.InterLinkModels.records_available(model)
        model.record_available = false
        history = run_simulation(run_args(spacecraft, model; duration=5.0); return_solution=true).prob.p.args.interlink_model.history
        @test all(sample -> isempty(sample.available), history) && isempty(contact_windows(history))
        @test only(contact_windows(history; active=true)).key == key                 # selected windows either way
    end

    @testset "shedding at the sender closes the link" begin
        # 0.001 Wh at 2 W lasts 1.8 s; tier 1 then sheds the transmitting radio.
        spacecraft, model = data_pair(sender=(capacity_wh=0.002, soc=0.5))
        solution = run_simulation(run_args(spacecraft, model; duration=10.0); return_solution=true)
        history = solution.prob.p.args.interlink_model.history
        closed = findfirst(entry -> isempty(entry.active), history)
        @test closed !== nothing && 1.0 < history[closed].time < 3.0
        @test all(entry -> isempty(entry.active) && isempty(entry.available), history[closed:end])
        received = solution.u[end].sc[2].data_received_bits
        @test 0 < received < shannon(SHANNON, chord(100e3)) * 3.0
        @test received ≈ solution.u[end].sc[1].data_sent_bits rtol=1e-12
        # A shed radio: level 0, its demand kept in tier 1, no draw, unavailable.
        spacecraft, model = data_pair()
        terminals = SM.InterlinkTerminals
        @test terminals.prepare_terminal_endpoints!(model, spacecraft)
        model.graph.active[1] = true
        @test terminals.link_terminal_power(model, 1, spacecraft[1]) == ((2.0, 0.0, 0.0, 0.0), 2.0)
        @test terminals.link_terminal_power(model, 2, spacecraft[2]) == ((0.0, 0.0, 0.0, 0.0), 0.0)
        SM.battery(spacecraft[1]).shed_tier = 1
        levels = terminals.link_driven_levels!(zeros(1), model, 1)
        @test levels == [1.0] && terminals.shed_link_driven_levels!(copy(levels), model, 1, spacecraft[1]) == [0.0]
        @test terminals.link_terminal_power(model, 1, spacecraft[1]) == ((2.0, 0.0, 0.0, 0.0), 0.0)
        @test terminals.terminal_availability(model, spacecraft) == [false, true]
    end

    @testset "scheduling" begin
        for (policy, selected) in ((SemimajorAxisRatePolicy(), false), (PowerDeliveryPolicy(), false), (DataThroughputPolicy(), true))
            spacecraft, model = data_pair()
            solution = run_simulation(run_args(spacecraft, model; duration=5.0, policy); return_solution=true)
            history = solution.prob.p.args.interlink_model.history
            @test all(entry -> entry.active == (selected ? [((1, 1), (2, 1))] : InterLinkKey[]), history)
            @test all(entry -> entry.available == [((1, 1), (2, 1))], history)
        end
        # Two-way: the two directions share the radios, so one of them is active at a time.
        spacecraft, model = data_pair()
        register_link!(model, spacecraft, (2, 1), (1, 1), LINK)
        solution = run_simulation(run_args(spacecraft, model; duration=5.0); return_solution=true)
        @test all(entry -> length(entry.active) == 1 && length(entry.available) == 2,
            solution.prob.p.args.interlink_model.history)
    end

    @testset "a fixed rate between virtual terminals" begin
        spacecraft = [radio_sat(1, 0.0; terminal=nothing, n_terminal=1, battery=false),
                      radio_sat(2, 100e3; terminal=nothing, n_terminal=1, battery=false)]
        model = InterlinkModel()
        register_link!(model, spacecraft, (1, 1), (2, 1), DataLink(transmit_power_w=0.0, budget=FixedRateBudget(1e6), range_m=1e6))
        t = run_simulation(run_args(spacecraft, model; duration=20.0); return_results=true).table
        @test t.sc2_data_received_bits ≈ 1e6 .* t.time rtol=1e-12 atol=1e-3
        @test !any(name -> occursin("power_", name), string.(propertynames(t)))
    end

    @testset "states and columns only with a data link" begin
        lasers = [radio_sat(i, (i - 1) * 5000.0; terminal=LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=100.0)) for i in 1:2]
        model = InterlinkModel()
        register_link!(model, lasers, (1, 1), (2, 1), LaserMomentumLink(power_w=40.0))
        solution = run_simulation(run_args(lasers, model; duration=1.0, policy=SemimajorAxisRatePolicy()); return_solution=true)
        @test !hasproperty(solution.u[1].sc[1], :data_sent_bits)
        names_laser = string.(propertynames(run_simulation(run_args(lasers, model; duration=1.0, policy=SemimajorAxisRatePolicy());
            return_results=true).table))
        @test !any(name -> occursin("data_", name), names_laser)
        @test !SM.has_data_links(model) && !SM.has_data_links(nothing)
        spacecraft, data_model = data_pair()
        @test SM.has_data_links(data_model)
    end

    @testset "checkpoints" begin
        mktempdir() do directory
            settings = SimulationSettings(results=false, generate_plots=false, checkpoint_enabled=true,
                checkpoint_directory=directory)
            spacecraft, model = data_pair()
            args = run_args(spacecraft, model; settings)
            state = (sc=[(pos=[RADIUS, 0.0, 0.0], interlink_dv=zeros(3), interlink_delta_sma=0.0, data_sent_bits=0.0,
                data_received_bits=0.0) for _ in 1:2],)
            SE._write_checkpoint!(args, 5.0, state, "tsit5"; is_active=[true, true])
            checkpoint = SE._load_checkpoint(args)
            @test isnothing(SE._restore_interlink_checkpoint_state!(args, checkpoint, trues(2)))
            changed = (
                "budget kind" => DataLink(transmit_power_w=1.0, budget=FixedRateBudget(1e6), range_m=2000e3),
                "bandwidth" => DataLink(transmit_power_w=1.0, budget=ShannonBudget(bandwidth_hz=2e6, frequency_hz=2.2e9),
                    range_m=2000e3),
                "frequency" => DataLink(transmit_power_w=1.0, budget=ShannonBudget(bandwidth_hz=1e6, frequency_hz=8.4e9),
                    range_m=2000e3),
                "transmit power" => DataLink(transmit_power_w=2.0, budget=SHANNON, range_m=2000e3),
            )
            for (name, link) in changed
                other, other_model = data_pair(; link)
                untouched = trues(2)
                @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(run_args(other, other_model; settings),
                    checkpoint, untouched)
                @test all(untouched)
            end
            more, more_model = data_pair()
            add_actuator!(more[2], more[2].root, radio())
            @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(run_args(more, more_model; settings),
                checkpoint, trues(2))
            identity = only(SE._interlink_checkpoint_configuration(args).hardware_terminals[1])
            @test identity.kind == (module_path=(:SpaceAGORA, :SimulationModel, :Components), name=:RadioTerminal,
                parameters=())
            @test keys(identity.fields) == fieldnames(RadioTerminal)
            @test identity.fields.gain_dbi.value == 6.0 && identity.fields.max_rate_bps.bits == reinterpret(UInt64, Inf)
            # A rebuilt radio with any field that sets the rate or the draw changed is refused.
            radios = (
                "transmitter gain" => (sender=(terminal=radio(gain_dbi=7.0),),),
                "transmitter efficiency" => (sender=(terminal=radio(efficiency=0.6),),),
                "transmitter limit" => (sender=(terminal=radio(max_power_w=4.0),),),
                "receiver noise temperature" => (receiver=(terminal=radio(noise_temperature_k=400.0),),),
                "receiver rate limit" => (receiver=(terminal=radio(max_rate_bps=1e6),),),
            )
            for (name, change) in radios
                other, other_model = data_pair(; change...)
                flags = falses(2)
                @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(run_args(other, other_model; settings),
                    checkpoint, flags)
                @test !any(flags)
            end
        end
        # A resumed run continues the counters.
        straight_pair = data_pair()
        straight = run_simulation(run_args(straight_pair...); return_results=true).table
        mktempdir() do directory
            segments(resume) = SimulationSettings(results=resume, results_directory=directory, generate_plots=false,
                checkpoint_enabled=true, checkpoint_interval_s=10.0, checkpoint_directory=directory,
                resume_from_checkpoint=resume)
            first_leg = run_simulation(run_args(data_pair()...; duration=30.0, settings=segments(false)); return_results=true).table
            resumed = run_simulation(run_args(data_pair()...; settings=segments(true)); return_results=true).table
            # It resumed at 30 s rather than starting again, and its table holds the
            # whole run: the first leg's rows, kept beside its checkpoint, then its own.
            @test TOML.parsefile(joinpath(directory, "simulation_results.manifest.toml"))["resumed_from_s"] == 30.0
            @test resumed.time[1] == 0.0
            @test isequal(resumed[resumed.time .< 30.0, :], first_leg[first_leg.time .< 30.0, :])
            @test resumed.sc2_data_received_bits[end] ≈ straight.sc2_data_received_bits[end] rtol=1e-9
            # The counter at the resume time holds what the first leg received.
            at_resume = findfirst(t -> t >= 30.0 - 1e-9, resumed.time)
            @test resumed.time[at_resume] ≈ 30.0
            @test resumed.sc2_data_received_bits[at_resume] ≈ shannon(SHANNON, chord(100e3)) * resumed.time[at_resume] rtol=1e-6
        end
    end
end

end # module
