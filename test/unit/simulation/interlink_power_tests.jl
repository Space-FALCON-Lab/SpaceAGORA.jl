module InterlinkPowerTests

# Power links: the Gaussian beam, the receiver, direction, occultation, the
# energy a beam moves between two batteries, curtailment, shedding at the
# sender, scheduling, columns and checkpoint identity.

using Test
using LinearAlgebra
using StaticArrays
using TOML
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PLANET = make_no_gram_planet(:earth)
const RADIUS = PLANET.Rp_e + 700e3
const SPEED = sqrt(PLANET.μ / RADIUS)

# A 100 kg bus `arc_m` along a circular orbit, with a laser transmitter or a
# power receiver and a battery.
function beam_sat(index, arc_m; sender, capacity_wh=10.0, soc=0.5, receiver_radius=0.25,
                  receiver_efficiency=0.4, max_power_w=Inf, battery=true, laser_efficiency=0.5)
    bus = Link(root=true, m=100.0)
    θ = arc_m / RADIUS
    ic = CartesianInitialCondition(SVector(RADIUS * cos(θ), RADIUS * sin(θ), 0.0),
        SVector(-SPEED * sin(θ), SPEED * cos(θ), 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, initial_condition=ic, id=index,
        n_terminal=0)
    if sender
        add_actuator!(sc, bus, LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=100.0, efficiency=laser_efficiency))
    else
        add_actuator!(sc, bus, PowerReceiver(normal=[0.0, -1.0, 0.0], radius_m=receiver_radius,
            efficiency=receiver_efficiency, max_power_w=max_power_w))
    end
    battery && add_power!(sc, bus, Battery(capacity_wh=capacity_wh, initial_soc=soc))
    return sc
end

function beam_pair(; arc_m=5000.0, link=PowerBeamLink(power_w=80.0, range_m=50e3), sender=(;), receiver=(;))
    spacecraft = [beam_sat(1, 0.0; sender=true, sender...), beam_sat(2, arc_m; sender=false, receiver...)]
    model = InterlinkModel()
    register_link!(model, spacecraft, (1, 1), (2, 1), link)
    return spacecraft, model
end

function run_args(spacecraft, model; duration=60.0, policy=PowerDeliveryPolicy(), settings=SimulationSettings(results=false, generate_plots=false))
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

@testset "Power links" begin
    @testset "the beam" begin
        link = PowerBeamLink(power_w=80.0, wavelength_m=1.064e-6, waist_radius_m=0.05, range_m=1e6)
        rayleigh = pi * 0.05^2 / 1.064e-6
        for d in (0.0, 1e3, rayleigh, 1e4, 1e5, 1e6)
            spot = 0.05 * sqrt(1 + (d / rayleigh)^2)
            @test beam_spot_radius_m(link, d) ≈ spot rtol=1e-14
            @test beam_capture_fraction(link, 0.25, d) ≈ 1 - exp(-2 * 0.25^2 / spot^2) rtol=1e-12
        end
        @test beam_spot_radius_m(link, rayleigh) ≈ 0.05 * sqrt(2)
        @test beam_capture_fraction(link, 0.25, 1e5) ≈ 0.2374 atol=1e-4
        @test beam_capture_fraction(link, 0.25, 1e6) < 0.003
        @test link_directed(link) && link_terminal_levels(link) == (80.0, 0.0)
        @test link_force(link, sample([RADIUS, 0.0, 0.0]), sample([RADIUS, 1e3, 0.0]), LinkContext()) == (zeros(3), zeros(3))
        for bad in ((power_w=-1.0,), (wavelength_m=0.0,), (waist_radius_m=Inf,), (range_m=0.0,))
            @test_throws ArgumentError PowerBeamLink(; merge((power_w=1.0, range_m=1.0), bad)...)
        end
        # The receiver delivers capture times efficiency, up to its limit.
        receiver = PowerReceiver(normal=[0.0, -1.0, 0.0], radius_m=0.25, efficiency=0.4)
        a, b = sample([RADIUS, 0.0, 0.0]), sample([RADIUS, 1e5, 0.0])
        @test link_received_power_w(link, receiver, a, b, LinkContext()) ≈ 80.0 * beam_capture_fraction(link, 0.25, 1e5) * 0.4
        capped = PowerReceiver(normal=[0.0, -1.0, 0.0], radius_m=0.25, efficiency=0.4, max_power_w=5.0)
        @test link_received_power_w(link, capped, a, b, LinkContext()) == 5.0
        @test link_received_power_w(LaserMomentumLink(), receiver, a, b, LinkContext()) == 0.0
        for bad in ((radius_m=0.0,), (efficiency=1.5,), (normal=[0.0, 0.0, 0.0],), (max_power_w=0.0,))
            @test_throws ArgumentError PowerReceiver(; merge((normal=[1.0, 0.0, 0.0], radius_m=1.0, efficiency=0.5), bad)...)
        end
    end

    @testset "direction" begin
        spacecraft, model = beam_pair()
        @test only(model.graph.edges) |> e -> map(v -> model.graph.vertices[v], e) == ((1, 1), (2, 1))
        link = PowerBeamLink(power_w=80.0, range_m=50e3)
        err = try
            register_link!(InterlinkModel(), spacecraft, (2, 1), (1, 1), link)        # receiver as transmitter
        catch caught
            caught
        end
        @test err isa ArgumentError && occursin("as its transmitter", sprint(showerror, err))
        receivers = [beam_sat(1, 0.0; sender=false), beam_sat(2, 5000.0; sender=false)]
        @test_throws ArgumentError register_link!(InterlinkModel(), receivers, (1, 1), (2, 1), link)
        senders = [beam_sat(1, 0.0; sender=true), beam_sat(2, 5000.0; sender=true)]
        @test_throws ArgumentError register_link!(InterlinkModel(), senders, (1, 1), (2, 1), link)   # laser as receiver
        @test register_link!(InterlinkModel(), senders, (2, 1), (1, 1), LaserMomentumLink()) == ((1, 1), (2, 1))
        @test_throws ArgumentError register_link!(InterlinkModel(), [senders[1], receivers[2]], (1, 1), (2, 1), LaserMomentumLink())
        @test terminal_supports(SM.mounted_link_terminals(receivers[1])[1].terminal, link, :receiver)
        @test !terminal_supports(SM.mounted_link_terminals(receivers[1])[1].terminal, link, :transmitter)
        @test link_terminal(SM.mounted_link_terminals(receivers[1])[1].terminal)
    end

    @testset "occultation" begin
        link = PowerBeamLink(power_w=1.0, range_m=1e8)
        ctx = LinkContext(planet=PLANET)
        near = (sample([RADIUS, 0.0, 0.0]), sample([RADIUS, 1e4, 0.0]))
        @test link_available(link, near..., ctx)
        opposite = (sample([RADIUS, 0.0, 0.0]), sample([-RADIUS, 0.0, 0.0]))
        @test !link_available(link, opposite..., ctx)                              # through the Earth
        @test link_available(link, opposite..., LinkContext())                     # no planet, no blocking
        x = 2PLANET.Rp_e
        tangent = (sample([-x, PLANET.Rp_e, 0.0]), sample([x, PLANET.Rp_e, 0.0]))
        @test !link_available(link, tangent..., ctx)                               # a grazing ray is blocked
        clear = (sample([-x, PLANET.Rp_e * (1 + 1e-6), 0.0]), sample([x, PLANET.Rp_e * (1 + 1e-6), 0.0]))
        @test link_available(link, clear..., ctx)
        @test !link_available(link, sample([RADIUS, 0.0, 0.0]), sample([RADIUS, 0.0, 0.0]), ctx)   # coincident
        @test !link_available(PowerBeamLink(power_w=1.0, range_m=10.0), near..., ctx)              # out of range
    end

    @testset "a beam moves energy between two batteries" begin
        spacecraft, model = beam_pair()
        result = run_simulation(run_args(spacecraft, model); return_results=true)
        t = result.table
        link = PowerBeamLink(power_w=80.0, range_m=50e3)
        capture = beam_capture_fraction(link, 0.25, 5000.0)
        @test t.sc1_battery_wh[1] - t.sc1_battery_wh[end] ≈ 80.0 / 0.5 * 60.0 / 3600.0 rtol=1e-9
        @test t.sc2_battery_wh[end] - t.sc2_battery_wh[1] ≈ 80.0 * capture * 0.4 * 60.0 / 3600.0 rtol=1e-6
        @test all(≈(80.0 * capture * 0.4; rtol=1e-6), t.sc2_power_received_w)
        @test t.sc2_power_generated_w == t.sc2_power_received_w
        @test all(==(80.0), t.sc1_power_beamed_w) && all(==(0.0), t.sc2_power_beamed_w)
        @test all(==(0.0), t.sc1_power_received_w)
        @test all(≈(160.0), t.sc1_power_load_w)
        # The battery integrates the net power the columns report (received power
        # less the loads), at both ends: the snapshot behind the battery's rate and
        # the one behind the columns both see the transmitter's state.
        integral(x, y) = sum((x[k + 1] - x[k]) * (y[k] + y[k + 1]) / 2 for k in 1:length(x) - 1)
        @test t.sc2_power_net_w ≈ t.sc2_power_received_w - t.sc2_power_load_w
        for (stored, net) in ((t.sc1_battery_wh, t.sc1_power_net_w), (t.sc2_battery_wh, t.sc2_power_net_w))
            @test stored[end] - stored[1] ≈ integral(t.time, net) / 3600.0 rtol=1e-8
        end
    end

    @testset "a full receiver curtails, and the sender still pays" begin
        spacecraft, model = beam_pair(receiver=(soc=1.0,))
        result = run_simulation(run_args(spacecraft, model); return_results=true)
        t = result.table
        @test all(≈(10.0; rtol=1e-12), t.sc2_battery_wh)
        @test t.sc1_battery_wh[1] - t.sc1_battery_wh[end] ≈ 80.0 / 0.5 * 60.0 / 3600.0 rtol=1e-9
        @test all(>(0), t.sc2_power_received_w)
    end

    @testset "shedding at the sender closes the beam" begin
        # 0.05 Wh at 160 W lasts about 1.1 s; tier 1 then sheds the transmitter.
        spacecraft, model = beam_pair(sender=(capacity_wh=0.1, soc=0.5))
        solution = run_simulation(run_args(spacecraft, model; duration=10.0); return_solution=true)
        history = solution.prob.p.args.interlink_model.history
        closed = findfirst(entry -> isempty(entry.active), history)
        @test closed !== nothing && 0.5 < history[closed].time < 2.0
        @test all(entry -> isempty(entry.active), history[closed:end])
        received = solution.u[end].sc[2].battery_wh - solution.u[1].sc[2].battery_wh
        link = PowerBeamLink(power_w=80.0, range_m=50e3)
        @test 0 < received < 80.0 * beam_capture_fraction(link, 0.25, 5000.0) * 0.4 * 2.0 / 3600.0
    end

    @testset "scheduling" begin
        spacecraft, model = beam_pair()
        gve = run_simulation(run_args(spacecraft, model; duration=10.0, policy=SemimajorAxisRatePolicy()); return_solution=true)
        @test all(entry -> isempty(entry.active), gve.prob.p.args.interlink_model.history)   # no force, no score
        spacecraft, model = beam_pair()
        power = run_simulation(run_args(spacecraft, model; duration=10.0); return_solution=true)
        @test all(entry -> entry.active == [((1, 1), (2, 1))], power.prob.p.args.interlink_model.history)
        # Out of range: available nowhere, selected never.
        spacecraft, model = beam_pair(arc_m=60e3)
        far = run_simulation(run_args(spacecraft, model; duration=10.0); return_solution=true)
        @test all(entry -> isempty(entry.active), far.prob.p.args.interlink_model.history)
    end

    @testset "columns only with a power beam" begin
        spacecraft = [beam_sat(1, 0.0; sender=true), beam_sat(2, 5000.0; sender=true)]
        model = InterlinkModel()
        register_link!(model, spacecraft, (1, 1), (2, 1), LaserMomentumLink(power_w=40.0))
        names_laser = string.(propertynames(run_simulation(run_args(spacecraft, model; duration=1.0, policy=SemimajorAxisRatePolicy());
            return_results=true).table))
        @test !any(name -> occursin("power_received", name) || occursin("power_beamed", name), names_laser)
        @test any(name -> occursin("power_generated", name), names_laser)
        # A receiver alone makes a spacecraft carry power parts.
        @test SM.has_power_parts(beam_sat(1, 0.0; sender=false, battery=false))
    end

    @testset "checkpoint identity" begin
        mktempdir() do directory
            settings = SimulationSettings(results=false, generate_plots=false, checkpoint_enabled=true,
                checkpoint_directory=directory)
            spacecraft, model = beam_pair()
            args = run_args(spacecraft, model; settings)
            state = (sc=[(pos=[RADIUS, 0.0, 0.0], interlink_dv=zeros(3), interlink_delta_sma=0.0) for _ in 1:2],)
            SE._write_checkpoint!(args, 5.0, state, "tsit5"; is_active=[true, true])
            checkpoint = SE._load_checkpoint(args)
            @test isnothing(SE._restore_interlink_checkpoint_state!(args, checkpoint, trues(2)))
            changed = (
                "beam power" => PowerBeamLink(power_w=81.0, range_m=50e3),
                "wavelength" => PowerBeamLink(power_w=80.0, wavelength_m=1.55e-6, range_m=50e3),
                "waist" => PowerBeamLink(power_w=80.0, waist_radius_m=0.06, range_m=50e3),
            )
            for (name, link) in changed
                other, other_model = beam_pair(; link)
                untouched = trues(2)
                @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(run_args(other, other_model; settings),
                    checkpoint, untouched)
                @test all(untouched)
            end
            # A receiver added to the second spacecraft changes its terminal hardware.
            more, more_model = beam_pair()
            add_actuator!(more[2], more[2].root, PowerReceiver(normal=[0.0, 1.0, 0.0], radius_m=0.1, efficiency=0.3))
            @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(run_args(more, more_model; settings),
                checkpoint, trues(2))
            identity = only(SE._interlink_checkpoint_configuration(args).hardware_terminals[2])
            @test identity.kind == (module_path=(:SpaceAGORA, :SimulationModel, :Components), name=:PowerReceiver,
                parameters=())
            @test keys(identity.fields) == fieldnames(PowerReceiver)
            @test identity.fields.radius_m.value == 0.25 && identity.fields.efficiency.value == 0.4
            @test identity.fields.max_power_w.value == Inf && identity.fields.max_power_w.bits == reinterpret(UInt64, Inf)
        end
    end

    @testset "checkpoints record each terminal's fields" begin
        # The receiver's aperture, efficiency and limit set the power it delivers,
        # and the transmitter's efficiency its draw. Through the checkpoint file, a
        # rebuilt terminal with any of them changed is refused before the run state
        # is restored, and the file still resumes the unchanged run.
        straight = run_simulation(run_args(beam_pair()...); return_results=true).table
        mktempdir() do directory
            segments(resume) = SimulationSettings(results=resume, results_directory=directory, generate_plots=false,
                checkpoint_enabled=true, checkpoint_interval_s=10.0, checkpoint_directory=directory,
                resume_from_checkpoint=resume)
            first_leg = run_simulation(run_args(beam_pair()...; duration=30.0, settings=segments(false)); return_results=true).table
            changed = (
                "receiver radius" => (receiver=(receiver_radius=0.3,),),
                "receiver efficiency" => (receiver=(receiver_efficiency=0.5,),),
                "receiver limit" => (receiver=(max_power_w=50.0,),),
                "transmitter efficiency" => (sender=(laser_efficiency=0.6,),),
            )
            for (name, change) in changed
                args = run_args(beam_pair(; change...)...; settings=segments(true))
                flags = falses(2)
                @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(args, SE._load_checkpoint(args), flags)
                @test !any(flags)                                                  # the activity flags not restored
                @test_throws ArgumentError run_simulation(args; return_results=true)
            end
            resumed = run_simulation(run_args(beam_pair()...; settings=segments(true)); return_results=true).table
            # It resumed at 30 s rather than starting again, and its table holds the
            # whole run: the first leg's rows, kept beside its checkpoint, then its own.
            @test TOML.parsefile(joinpath(directory, "simulation_results.manifest.toml"))["resumed_from_s"] == 30.0
            @test resumed.time[1] == 0.0 && resumed.time[end] > 30.0
            @test isequal(resumed[resumed.time .< 30.0, :], first_leg[first_leg.time .< 30.0, :])
            @test resumed.sc2_battery_wh[end] ≈ straight.sc2_battery_wh[end] rtol=1e-9
            @test resumed.sc1_battery_wh[end] ≈ straight.sc1_battery_wh[end] rtol=1e-9
        end
    end
end

end # module
