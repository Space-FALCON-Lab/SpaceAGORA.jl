module ThrusterGimbalTests
using Test, StaticArrays, LinearAlgebra, SpaceAGORA
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel

mutable struct HoldLevels <: AbstractControlEffectorModel
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::HoldLevels) = :mounted
SpaceAGORA.actuator_levels(c::HoldLevels, i::Int) = i == 1 ? c.levels : nothing
SpaceAGORA.update_control!(::HoldLevels, u, p, t::Float64, i::Int64) = nothing

function far_spacecraft(; mass=500.0)
    root = SM.Link(root=true, m=mass, dims=MVector(1.0, 1.0, 1.0))
    ic = SM.CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), SVector(0.0, 0.0, 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    return SM.SpacecraftModel(links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia, initial_condition=ic, id=1), root
end

function held_run(sc, levels; mission_time=5.0, rate=1.0)
    base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=mission_time,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
        keplerian=false, EI_km=1.0, verbose=false, results=false)
    args = SM.SimConfig._with_configuration(base;
        control_model=SM.ControlModel(control_effectors=(HoldLevels(collect(Float64, levels)),), control_rates=[rate]))
    return run_simulation(args; isolate_state=false, return_results=true, return_solution=true)
end

# Holds the thrust on and turns the gimbal to `angle` at the first tick at or after `at_s`.
mutable struct StepLevels <: AbstractControlEffectorModel
    angle::Float64
    at_s::Float64
    levels::Vector{Float64}
end
StepLevels(angle, at_s) = StepLevels(angle, at_s, [1.0, 0.0])
SpaceAGORA.actuation_mode(::StepLevels) = :mounted
SpaceAGORA.actuator_levels(c::StepLevels, i::Int) = i == 1 ? c.levels : nothing
SpaceAGORA.update_control!(c::StepLevels, u, p, t::Float64, i::Int64) = (c.levels[2] = t >= c.at_s ? c.angle : 0.0; nothing)

engine(; location=(0.0, 0.0, -1.0)) = Thruster(max_thrust=2.0, direction=MVector(1.0, 0.0, 0.0), location=MVector(location...))

@testset "Thruster gimbals" begin
    @testset "the record and mounting" begin
        sc, root = far_spacecraft()
        jet = add_actuator!(sc, root, engine())
        @test_throws ArgumentError ThrusterGimbal(jet; axis=SVector(0.0, 0.0, 0.0))
        @test_throws ArgumentError ThrusterGimbal(jet; axis=SVector(0.0, 1.0, 0.0), min_angle_rad=0.2, max_angle_rad=0.1)
        @test_throws ArgumentError ThrusterGimbal(jet; axis=SVector(0.0, 1.0, 0.0), max_angle_rad=0.1, angle_rad=0.2)
        loose = Thruster(max_thrust=1.0, direction=MVector(1.0, 0.0, 0.0))
        @test_throws ArgumentError add_actuator!(sc, root, ThrusterGimbal(loose; axis=SVector(0.0, 1.0, 0.0)))   # its thruster is not on the link
        g = add_actuator!(sc, root, ThrusterGimbal(jet; axis=SVector(0.0, 2.0, 0.0)))
        @test_throws ArgumentError add_actuator!(sc, root, g)                                                     # mounted twice
        @test g.axis == SVector(0.0, 1.0, 0.0) && jet.gimballed
        @test [e.actuator for e in mounted_actuators(sc)] == [jet, g] && SM.mounted_gimbal_count(sc) == 1
        @test actuator_limits(g) == (-Inf, Inf)
        plain = Thruster(max_thrust=1.0, direction=MVector(0.0, 1.0, 0.0))
        @test thrust_direction(plain) === SVector(0.0, 1.0, 0.0)                                                 # no gimbal: the direction itself
    end

    @testset "a gimbal turns the thrust about its axis, not the thrust point" begin
        θ = 0.1
        sc, root = far_spacecraft()
        jet = add_actuator!(sc, root, engine())
        g = add_actuator!(sc, root, ThrusterGimbal(jet; axis=SVector(0.0, 1.0, 0.0), angle_rad=θ))
        @test thrust_direction(jet) ≈ SVector(cos(θ), 0.0, -sin(θ)) atol=1e-15
        force, torque = actuator_wrench(jet, 2.0, actuator_context(root))
        @test force ≈ 2.0 * SVector(cos(θ), 0.0, -sin(θ)) atol=1e-15
        @test torque ≈ cross(SVector(0.0, 0.0, -1.0), force) atol=1e-15                                     # about the fixed nozzle
        B = allocation_matrix(sc)
        @test B[1:3, 1] ≈ SVector(cos(θ), 0.0, -sin(θ)) && all(iszero, B[:, 2])                             # the gimbal has no column of its own
        # A second axis, about z, turns after the first.
        φ = 0.3
        add_actuator!(sc, root, ThrusterGimbal(jet; axis=SVector(0.0, 0.0, 1.0), angle_rad=φ))
        @test thrust_direction(jet) ≈ SVector(cos(φ) * cos(θ), sin(φ) * cos(θ), -sin(θ)) atol=1e-15
        @test applied_actuator_levels(sc, [2.0, 0.0, 0.0]) == [2.0, θ, φ]
    end

    @testset "a gimbal moves at control ticks within its angle and rate limits" begin
        sc, root = far_spacecraft()
        jet = add_actuator!(sc, root, engine())
        slow = add_actuator!(sc, root, ThrusterGimbal(jet; axis=SVector(0.0, 1.0, 0.0), max_angle_rad=0.3, max_rate_rad_s=0.05))
        SM.apply_gimbal_levels!(sc, [1.0, 1.0], 0.0)
        @test slow.angle_rad == 0.0                                                                          # a rate-limited first tick holds
        SM.apply_gimbal_levels!(sc, [1.0, 1.0], 2.0)
        @test slow.angle_rad ≈ 0.1 && thrust_direction(jet) ≈ SVector(cos(0.1), 0.0, -sin(0.1))
        SM.apply_gimbal_levels!(sc, [1.0, 1.0], 100.0)
        @test slow.angle_rad == 0.3                                                                          # within its angle limit
        @test SM.gimbal_step(slow, -1.0, 101.0) ≈ 0.25
    end

    @testset "a run vectors the thrust from the gimbal's first tick" begin
        # 1 N along x from the centre of mass, the gimbal commanded to 0.2 rad: the
        # thrust is along x until the first control tick at 1 s, then turned.
        θ = 0.2
        sc, root = far_spacecraft(mass=500.0)
        jet = add_actuator!(sc, root, engine(location=(0.0, 0.0, 0.0)))
        g = add_actuator!(sc, root, ThrusterGimbal(jet; axis=SVector(0.0, 1.0, 0.0)))
        res = held_run(sc, [1.0, θ])
        v = res.solution.u[end].sc[1].vel
        @test v[1] ≈ (1.0 + 4.0 * cos(θ)) / 500.0 rtol=1e-9
        @test v[3] ≈ -4.0 * sin(θ) / 500.0 rtol=1e-9
        @test norm(res.solution.u[end].sc[1].ω) < 1e-15                                                      # no torque about the centre of mass
        @test res.table.sc1_actuator_level_2[end] ≈ θ                                                        # the gimbal's column is its angle
    end

    @testset "the visualization saves a gimballed thruster's exhaust direction" begin
        θ = 0.2
        sc, root = far_spacecraft()
        jet = add_actuator!(sc, root, engine(location=(0.0, 0.0, 0.0)))
        add_actuator!(sc, root, ThrusterGimbal(jet; axis=SVector(0.0, 1.0, 0.0)))
        base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=3.0,
            initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
            keplerian=false, EI_km=1.0, verbose=false, results=false)
        settings = SM.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false,
            save_visualization_scene=true)
        args = SM.SimConfig._with_configuration(base; simulation_settings=settings,
            control_model=SM.ControlModel(control_effectors=(HoldLevels([1.0, θ]),), control_rates=[1.0]))
        @test :thruster_direction in [f.name for f in SM.SimulationCallbacks.visualization_save_fields(args)]
        table = run_simulation(args; isolate_state=false, return_results=true).table
        exhaust(t) = [table[t, "sc1_thruster_direction_$(c)"] for c in 1:3]
        @test exhaust(1) ≈ [-1.0, 0.0, 0.0]                                                                  # nominal before the first tick
        @test exhaust(size(table, 1)) ≈ -collect(thrust_direction(jet)) && abs(exhaust(size(table, 1))[3]) > 0.1
        # Ticks every 5 s and a sample every 0.5 s, so a step spans many samples:
        # the samples of the step that ends at the turning tick (10 s) are saved
        # after it, and keep the nominal direction the step flew with.
        turned, troot = far_spacecraft()
        tjet = add_actuator!(turned, troot, engine(location=(0.0, 0.0, 0.0)))
        add_actuator!(turned, troot, ThrusterGimbal(tjet; axis=SVector(0.0, 1.0, 0.0)))
        base20 = make_example_config(planet=SM.Earth(), spacecraft=turned, mission_time=20.0,
            initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
            keplerian=false, EI_km=1.0, verbose=false, results=false)
        stepped = SM.SimConfig._with_configuration(base20; simulation_settings=settings,
            mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, mission_time=20.0, number_of_orbits=1,
                keplerian=false, orientation_sim=true, num_steps_to_save=1000, data_rate=0.5),
            control_model=SM.ControlModel(control_effectors=(StepLevels(θ, 10.0),), control_rates=[5.0]))
        steps = run_simulation(stepped; isolate_state=false, return_results=true).table
        rows(f) = [r for r in 1:size(steps, 1) if f(steps.time[r])]
        sexhaust(r) = [steps[r, "sc1_thruster_direction_$(c)"] for c in 1:3]
        @test length(rows(t -> 5.0 < t < 10.0)) >= 9
        @test all(r -> sexhaust(r) ≈ [-1.0, 0.0, 0.0], rows(t -> t < 10.0))                  # before the tick
        @test all(r -> sexhaust(r) ≈ [-cos(θ), 0.0, sin(θ)], rows(t -> t >= 10.0))            # from the tick on
        plain, proot = far_spacecraft()
        add_actuator!(plain, proot, engine())
        @test !SM.SimulationCallbacks._thruster_direction_field_enabled(SM.SimConfig._with_configuration(args;
            dynamics_model=SM.DynamicsModel([plain], ())))                                                    # no gimbal, no field
    end

    @testset "a resumed run continues a gimbal's turn from the checkpoint" begin
        # Commanded to 0.2 rad at 0.03 rad/s with ticks every 0.3 s: still turning
        # at the 3 s checkpoint. Compared with a run checkpointed at the same
        # interval, whose segments start where the resumed run's do.
        function turning_run(; mission_time, settings=(;))
            sc, root = far_spacecraft()
            jet = add_actuator!(sc, root, engine(location=(0.0, 0.0, 0.0)))
            g = add_actuator!(sc, root, ThrusterGimbal(jet; axis=SVector(0.0, 1.0, 0.0), max_rate_rad_s=0.03))
            base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=mission_time,
                initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
                dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
                keplerian=false, EI_km=1.0, verbose=false, results=false)
            st = base.simulation_settings
            fields = NamedTuple{fieldnames(typeof(st))}(map(n -> getfield(st, n), fieldnames(typeof(st))))
            args = SM.SimConfig._with_configuration(base; simulation_settings=SM.SimulationSettings(; merge(fields, settings)...),
                control_model=SM.ControlModel(control_effectors=(HoldLevels([1.0, 0.2]),), control_rates=[0.3]))
            run_simulation(args; isolate_state=false)
            return g, jet
        end
        checkpointed(dir; kw...) = (checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=dir, kw...)
        straight, sjet = turning_run(mission_time=6.0, settings=checkpointed(mktempdir()))
        dir = mktempdir()
        turning_run(mission_time=3.0, settings=checkpointed(dir))
        resumed, rjet = turning_run(mission_time=6.0, settings=checkpointed(dir; resume_from_checkpoint=true))
        @test 0.0 < straight.angle_rad < 0.2                                                                 # still turning at the end
        @test resumed.angle_rad == straight.angle_rad && resumed.last_update_s == straight.last_update_s
        @test thrust_direction(rjet) == thrust_direction(sjet)
    end

    @testset "a resumed run turns a gimbal as a run in one piece, ticking at the checkpoints" begin
        # Ticks every 0.5 s fall on the 1 s checkpoint boundaries and on the 3 s
        # resume time, which the run that ended there did not tick at.
        function turning_run(; mission_time, settings=(;))
            sc, root = far_spacecraft()
            jet = add_actuator!(sc, root, engine(location=(0.0, 0.0, 0.0)))
            g = add_actuator!(sc, root, ThrusterGimbal(jet; axis=SVector(0.0, 1.0, 0.0), max_rate_rad_s=0.03))
            base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=mission_time,
                initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
                dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
                keplerian=false, EI_km=1.0, verbose=false, results=false)
            st = base.simulation_settings
            fields = NamedTuple{fieldnames(typeof(st))}(map(n -> getfield(st, n), fieldnames(typeof(st))))
            args = SM.SimConfig._with_configuration(base; simulation_settings=SM.SimulationSettings(; merge(fields, settings)...),
                control_model=SM.ControlModel(control_effectors=(HoldLevels([1.0, 0.2]),), control_rates=[0.5]))
            run_simulation(args; isolate_state=false)
            return g, jet
        end
        checkpointed(dir; kw...) = (checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=dir, kw...)
        straight, sjet = turning_run(mission_time=6.0)                                                     # in one piece
        split, _ = turning_run(mission_time=6.0, settings=checkpointed(mktempdir()))
        dir = mktempdir()
        early, _ = turning_run(mission_time=3.0, settings=checkpointed(dir))
        resumed, rjet = turning_run(mission_time=6.0, settings=checkpointed(dir; resume_from_checkpoint=true))
        @test early.last_update_s == 2.5 && straight.last_update_s == 5.5
        @test 0.0 < straight.angle_rad < 0.2                                                                 # still turning at the end
        @test split.angle_rad == straight.angle_rad && split.last_update_s == straight.last_update_s
        @test resumed.angle_rad == straight.angle_rad && resumed.last_update_s == straight.last_update_s
        @test thrust_direction(rjet) == thrust_direction(sjet)
    end

    @testset "a failed gimbal holds its angle" begin
        sc, root = far_spacecraft()
        jet = add_actuator!(sc, root, engine())
        g = add_actuator!(sc, root, ThrusterGimbal(jet; axis=SVector(0.0, 1.0, 0.0), angle_rad=0.05))
        add_fault!(sc, ActuatorFault(g; kind=:failed, start_s=0.0))
        held_run(sc, [1.0, 0.4]; mission_time=3.0)
        @test g.angle_rad == 0.05 && thrust_direction(jet) ≈ SVector(cos(0.05), 0.0, -sin(0.05))
    end
end

end # module ThrusterGimbalTests
