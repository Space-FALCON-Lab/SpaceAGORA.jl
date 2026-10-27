module RHSKernelParityTests
using Test, SpaceAGORA, StaticArrays, ComponentArrays
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine

function configuration(orientation)
    planet = SM.make_no_gram_planet(:earth)
    spacecraft = map(1:4) do i
        root = SM.Link(root=true, m=400.0 + i, ref_area=4.0)
        initial = SM.CartesianInitialCondition(
            SVector(planet.Rp_e + 120_000.0 + 100.0 * i, 0.0, 0.0),
            SVector(0.0, 7_800.0 + i, 0.0);
            q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.01, -0.02, 0.03))
        SM.SpacecraftModel(; joints=SM.Joint[], links=[root], root,
            initial_condition=initial, inertia_tensor=root.inertia)
    end
    base = SpaceAGORA.TelemetryVerification.make_example_config(
        planet=planet, spacecraft=first(spacecraft), mission_time=10.0,
        initial_time=SM.InitialTime(year=2021, month=3, day=4, hour=5),
        dynamic_effectors=(SM.InverseSquaredGravityModel(), SM.AerodynamicCoefficientfM()),
        density_model=SM.ExponentialAtmosphereModel(planet),
        ephemerides_model=SM.SimpleEphemeridesModel(), EI_km=300.0,
        orientation_sim=orientation, results=false, verbose=false)
    return SM.SimConfig._with_configuration(base;
        dynamics_model=SM.DynamicsModel(spacecraft, base.dynamics_model.dynamic_effectors))
end

function parameters(args, mode)
    p = SM.ODEParams(n_sats=4, args=deepcopy(args))
    p.is_active .= true
    p.is_active[2] = false
    p.shared_buffers.rhs_env_config[] = withenv(
        "SPACEAGORA_RHS_EXECUTION_MODE" => mode,
        "SPACEAGORA_RHS_BATCH_PARALLEL" => (mode == "serial" ? "off" : "on"),
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off",
        "SPACEAGORA_PARALLEL_PROFILE" => "",
    ) do
        SE._snapshot_rhs_plan_env_config()
    end
    return p
end

const VARIANTS = (
    (:full, SE._spacecraft_dynamics_dispatch!),
    (:slow, SE.spacecraft_dynamics_slow!),
    (:implicit, SE.spacecraft_dynamics_implicit_atmosphere!),
    (:explicit, SE.spacecraft_dynamics_explicit_remainder!),
    (:fast_control, SE.spacecraft_dynamics_fast_control!),
)

@testset "serial and satellite kernels preserve each physics partition" begin
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())) do
        for orientation in (false, true)
            args = configuration(orientation)
            u = SE.build_initial_conditions(args)
            for (name, rhs!) in VARIANTS
                @testset "$name orientation=$orientation" begin
                    serial = parameters(args, "serial")
                    batch = parameters(args, "satellite")
                    @test !SE._rhs_batch_parallel_enabled(serial, 4)
                    if Threads.nthreads() > 1 && SE.Polyester.num_cores() > 1
                        @test SE._rhs_batch_parallel_enabled(batch, 4)
                    else
                        @test_skip "multiple threads required to exercise satellite batching"
                    end
                    du_serial = fill!(zero(u), 123.0)
                    du_batch = fill!(zero(u), 456.0)
                    rhs!(du_serial, u, serial, 0.0)
                    rhs!(du_batch, u, batch, 0.0)
                    @test all(isfinite, du_serial)
                    @test du_batch ≈ du_serial rtol=1e-13 atol=1e-13
                    @test all(iszero, du_serial.sc[2])
                    @test all(iszero, du_batch.sc[2])
                    @test serial.shared_buffers.current_time[] == batch.shared_buffers.current_time[] == 0.0
                    if name in (:implicit, :fast_control)
                        @test all(iszero, du_serial.sc[1].heat_loads)
                        orientation && @test all(iszero, du_serial.sc[1].q)
                    else
                        @test du_serial.sc[1].pos ≈ u.sc[1].vel
                        @test any(!iszero, du_serial.sc[1].vel)
                    end
                end
            end
        end
    end
end
end
