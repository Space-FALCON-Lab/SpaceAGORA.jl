module InitialForceOutputTests
# This module tests values and query counts, not timing. Lower optimization
# avoids pathological LLVM scheduling for the large inferred test fixtures;
# production benchmark processes keep their normal optimization settings.
Base.Experimental.@optlevel 1
using Test
using LinearAlgebra
using StaticArrays
using DiffEqCallbacks: SavedValues
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const CB = SM.SimulationCallbacks
const EM = SM.EnvironmentModels
const CH = SM.ControlHooks
const ZERO3 = SVector{3, Float64}(0.0, 0.0, 0.0)

# A wrapper with the public atmosphere/wrench contracts must not need to be
# named in a concrete-type allowlist to save correct initial forces.
struct WrappedAero <: SM.AbstractForceTorqueModel
    model::AerodynamicCoefficientConstant
end
SM.environment_requirements(m::WrappedAero) = SM.environment_requirements(m.model)
SM.solver_partition(m::WrappedAero) = SM.solver_partition(m.model)
SM.wrench(m::WrappedAero, x::SM.StateSample, env::SM.EnvironmentSample, t::Float64) =
    SM.wrench(m.model, x, env, t)
SM.wrench_caching!(m::WrappedAero, x::SM.StateSample, env::SM.EnvironmentSample,
        t::Float64, p::ODEParams, i::Int) = SM.wrench_caching!(m.model, x, env, t, p, i)

struct ExplicitAtmosphereEffector <: SM.AbstractForceTorqueModel end
SM.environment_requirements(::ExplicitAtmosphereEffector) = SM.EffectorEnvironmentRequirements(atmosphere=true)
SM.solver_partition(::ExplicitAtmosphereEffector) = :explicit

mutable struct CountingAtmosphere <: SM.AbstractDensityModel
    queries::Int
end
EM.density_model_history_dependent(::CountingAtmosphere) = true
function EM.getDensity(model::CountingAtmosphere, h::Float64, lat::Float64,
        lon::Float64, t::Float64, wind::Bool, p)
    model.queries += 1
    return 1.0e-11, 800.0, wind ? SVector(0.1 * model.queries, 0.0, 0.0) : ZERO3
end

mutable struct CountingControl <: SM.AbstractControlEffectorModel
    reads::Int
    updates::Int
end
function CH.calcControlForceTorque(model::CountingControl, u, p, i::Int64, t::Float64)
    model.reads += 1
    return ZERO3, ZERO3
end
CH.calcControlEffect!(model::CountingControl, u, p, t::Float64, i::Int64) =
    (model.updates += 1; nothing)

Base.@noinline function configuration(; model=ExponentialAtmosphereModel(1e-11, 550e3, 50e3;
        temperature_k=800.0), orientation=false, per_link=false, controls=(),
        mode=:tsit5, wrapped=false)
    planet = Earth()
    root = Link(root=true, m=500.0, ref_area=12.0)
    links = per_link ? [root, Link(root=false, m=1.0, ref_area=1.0,
        r=MVector{3, Float64}(0.0, 0.0, 2.0))] : [root]
    ic = InitialCondition(ra=planet.Rp_e + 550e3, rp=planet.Rp_e + 550e3,
        i=53.0, ω=0.0, Ω=10.0, ν=0.0,
        q=SVector{4, Float64}(0.0, sind(20.0), 0.0, cosd(20.0)))
    craft = SpacecraftModel(Joint[], links, root, true, sum(b.m for b in links),
        0.0, root.inertia, 0, 0, ic, 1)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=0.2, orientation_sim=orientation, data_rate=0.1),
        environment_model=EnvironmentModel(planet=planet, EI=300.0,
            density_model=model, thermal_model=MaxwellianHeat(
                thermal_accomodation_factor=1.0, planet=planet),
            wind=model isa CountingAtmosphere, topography=false,
            ephemerides_model=SimpleEphemeridesModel()),
        dynamics_model=DynamicsModel([craft], (InverseSquaredGravityModel(),
            (wrapped ? WrappedAero(AerodynamicCoefficientConstant()) :
                per_link ? AerodynamicCoefficientfM(per_link_atmosphere=true) : AerodynamicCoefficientConstant()))),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=controls,
            control_rates=fill(0.1, length(controls))),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9,
            abstol_orbit=1e-9, dt_max_orbit=0.1),
        solver_config=SolverConfig(solver_mode=mode, parallel=false))
end

Base.@noinline function parameters(args)
    u = SE.build_initial_conditions(args)
    p = ODEParams(n_sats=1, args=args)
    SE._initialize_runtime_env_config!(p)
    SE._initialize_in_atmosphere_flags!(p, u)
    SE._initialize_save_cache_buffers!(p)
    SE._initialize_heat_rate_buffers!(p)
    p.shared_buffers.et_start[] = SM.ephemerides_time_seconds(
        args.initial_time, args.environment_model.ephemerides_model)
    return u, p
end

Base.@noinline function arm!(u, p, t=0.0)
    cb = CB.get_initial_force_output_callback(p.args.dynamics_model.dynamic_effectors)
    cb.initialize(cb, u, t, (u=u, p=p, t=t))
    return nothing
end

# These tests inspect values and query counts, never timing. Keep the large
# production RHS out of the test caller's optimization unit; otherwise Julia
# can inline several copies into one testset and spend minutes in LLVM.
Base.@noinline evaluate_rhs!(du, u, p, t) =
    Base.invokelatest(SE.spacecraft_dynamics!, du, u, p, t)

# Independent constant-coefficient root-only expectation: CD=2.2 at max-drag
# incidence, F=-rho*CD*A*|v_relative|*v_relative/2. No RHS or aero wrench call.
Base.@noinline function expected_drag(u, p, t=0.0)
    env = p.args.environment_model
    frame = SE.sample_planet_frame(u.sc[1], p, 1, t)
    model = env.density_model
    rho = model.ρ_ref * exp(-(frame.alt_m - model.h_ref) / model.H)
    area = p.args.dynamics_model.spacecraft[1].root.ref_area
    return frame.l_pi' * (-0.5 * rho * 2.2 * area * norm(frame.vel_pp) * frame.vel_pp)
end

Base.@noinline function test_initial_publication()
    args = configuration()
    u, p = Base.invokelatest(parameters, args)
    stale = SVector{3, Float64}(41.0, -7.0, 12.0)
    fill!(p.save_cache.drag_cache, stale)
    # A calibration-like RHS before callback initialization publishes nothing.
    evaluate_rhs!(zero(u), u, p, 0.0)
    @test !p.save_cache.initial_force_output_pending[]
    arm!(u, p)
    view = (p=p,)
    drag = CB._save_drag(1, u, 0.0, view)
    lift = CB._save_lift(1, u, 0.0, view)
    cross_force = CB._save_cross(1, u, 0.0, view)
    @test all(isnan, only(drag))
    # Neither a different time nor a finite-difference Jacobian state may fill it.
    evaluate_rhs!(zero(u), u, p, 0.01)
    @test all(isnan, only(drag))
    perturbed = copy(u)
    perturbed.sc[1].pos[1] += 1.0
    evaluate_rhs!(zero(u), perturbed, p, 0.0)
    @test all(isnan, only(drag))
    evaluate_rhs!(zero(u), u, p, 0.0)
    @test only(drag) ≈ expected_drag(u, p) rtol=5e-14
    @test norm(only(drag)) > 0.0
    @test only(lift) == ZERO3 && only(cross_force) == ZERO3
    @test !p.save_cache.initial_force_output_pending[]
    @test isempty(p.save_cache.initial_force_output_destinations)
    @test p.save_cache.initial_force_output_state[] === nothing
    saved_drag = copy(drag)
    fill!(p.save_cache.drag_cache, stale)
    @test drag == saved_drag # no alias to later cache values
    # Reinitialization at a nonzero checkpoint time does not rewrite old arrays.
    arm!(u, p, 3.0)
    next_drag = CB._save_drag(1, u, 3.0, view)
    @test all(isnan, only(next_drag))
    evaluate_rhs!(zero(u), u, p, 3.0)
    @test only(next_drag) ≈ expected_drag(u, p, 3.0) rtol=5e-14
    @test drag == saved_drag
    # Abort before a matching RHS leaves explicitly unavailable output,
    # independent of calibration contents, rather than reporting a false zero.
    arm!(u, p, 4.0)
    aborted = CB._save_drag(1, u, 4.0, view)
    @test all(isnan, only(aborted))
end

Base.@noinline function test_atmospheric_wrapper()
    args = configuration(wrapped=true)
    u, p = Base.invokelatest(parameters, args)
    arm!(u, p)
    drag = CB._save_drag(1, u, 0.0, (p=p,))
    evaluate_rhs!(zero(u), u, p, 0.0)
    @test only(drag) ≈ expected_drag(u, p) rtol=5e-14
end

Base.@noinline function test_specialized_driver_scope()
    effectors = configuration().dynamics_model.dynamic_effectors
    marker = CB.get_initial_force_output_callback(effectors)
    for mode in (:multirate, :gravity_backbone_split)
        @test CB.get_initial_force_output_callback(effectors; solver_mode=mode) === nothing
        # Explicit SolverConfig must win over a conflicting process default.
        withenv("SPACEAGORA_SOLVER_MODE" => "tsit5") do
            args = configuration(mode=mode)
            callbacks = CB.get_callbacks(1, effectors, args)
            @test !any(cb -> typeof(cb.initialize) === typeof(marker.initialize),
                callbacks.discrete_callbacks)
        end
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split") do
        args = configuration(mode=:tsit5)
        callbacks = CB.get_callbacks(1, effectors, args)
        @test any(cb -> typeof(cb.initialize) === typeof(marker.initialize),
            callbacks.discrete_callbacks)
    end
    mixed = (effectors..., ExplicitAtmosphereEffector())
    @test CB.get_initial_force_output_callback(mixed; solver_mode=:split_imex) === nothing
    @test CB.get_initial_force_output_callback(mixed; solver_mode=:tsit5) !== nothing
end

Base.@noinline function test_stateful_queries()
    for freeze in ("off", "auto"), per_link in (false, true)
        withenv("SPACEAGORA_DENSITY_FREEZE_PER_STEP" => freeze,
                "SPACEAGORA_RHS_EXECUTION_MODE" => "satellite",
                "SPACEAGORA_INNER_THREAD_BUDGET" => "1") do
            reference_model, observed_model = CountingAtmosphere(0), CountingAtmosphere(0)
            reference_control, observed_control = CountingControl(0, 0), CountingControl(0, 0)
            reference_args = configuration(model=reference_model, orientation=true,
                per_link=per_link, controls=(reference_control,))
            observed_args = configuration(model=observed_model, orientation=true,
                per_link=per_link, controls=(observed_control,))
            ur, pr = Base.invokelatest(parameters, reference_args)
            u, p = Base.invokelatest(parameters, observed_args)
            # Model the already-run physical density initializer without a new
            # density query when output is armed. Freeze mode uses this sample.
            for params in (pr, p)
                CB._write_density_buffers!(params, 1, 1e-11, 800.0, ZERO3, 0.0)
            end
            arm!(u, p)
            drag = CB._save_drag(1, u, 0.0, (p=p,))
            lift = CB._save_lift(1, u, 0.0, (p=p,))
            cross_force = CB._save_cross(1, u, 0.0, (p=p,))
            @test observed_model.queries == observed_control.reads == observed_control.updates == 0
            dur, du = zero(ur), zero(u)
            evaluate_rhs!(dur, ur, pr, 0.0)
            evaluate_rhs!(du, u, p, 0.0)
            @test du == dur
            @test drag == pr.save_cache.drag_cache
            @test lift == pr.save_cache.lift_cache
            @test cross_force == pr.save_cache.cross_cache
            @test observed_model.queries == reference_model.queries
            @test observed_control.reads == reference_control.reads
            @test observed_control.updates == reference_control.updates == 0
            if freeze == "auto"
                @test observed_model.queries == 0
            else
                @test observed_model.queries > 0
            end
        end
    end
end

Base.@noinline function test_solver_paths()
    withenv("SPACEAGORA_RHS_CALIBRATE" => "off",
            "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0",
            "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0") do
        for mode in (:tsit5, :rodas5p, :split_imex), recorder_first in (false, true)
            args = configuration(mode=mode)
            u, p = Base.invokelatest(parameters, args)
            expected = expected_drag(u, p)
            rec = TrajectoryRecorder(args)
            saved = SavedValues(Float64, SM.SaveData)
            standard = CB.get_data_saving_callback(1, args, default_save_fields(args), saved)
            recorder = get_trajectory_recorder_callback(rec)
            custom = TrajectoryRecorder(args; save_fields=[SaveField(:drag,
                (u, t, integrator) -> [SVector(9.0, 8.0, 7.0)]; per_satellite=true)])
            callbacks = recorder_first ? (recorder, standard) : (standard, recorder)
            result = Base.invokelatest(run_simulation, args; isolate_state=false,
                return_solver_metadata=true,
                extra_callbacks=(callbacks..., get_trajectory_recorder_callback(custom)))
            @test result.retcode == "Success"
            @test first(saved.t) == first(trajectory_times(rec)) == 0.0
            @test saved.saveval[1][:drag][1] ≈ expected rtol=5e-14
            @test trajectory_field(rec, :drag)[:, 1, 1] ≈ expected rtol=5e-14
            @test saved.saveval[1][:drag][1] == SVector{3, Float64}(trajectory_field(rec, :drag)[:, 1, 1])
            @test all(iszero, trajectory_field(rec, :lift)[:, :, 1])
            @test all(iszero, trajectory_field(rec, :cross)[:, :, 1])
            @test trajectory_field(custom, :drag)[1] == [SVector(9.0, 8.0, 7.0)]
        end
    end
end

println("initial_force_output: initial force publication uses actual initial RHS only")
flush(stdout)
@testset "initial force publication uses actual initial RHS only" begin
    Base.invokelatest(test_initial_publication)
end
flush(stdout)

println("initial_force_output: atmospheric wrappers also publish initial force output")
flush(stdout)
@testset "atmospheric wrappers also publish initial force output" begin
    Base.invokelatest(test_atmospheric_wrapper)
end
flush(stdout)

println("initial_force_output: specialized drivers preserve their existing output path")
flush(stdout)
@testset "specialized drivers preserve their existing output path" begin
    Base.invokelatest(test_specialized_driver_scope)
end
flush(stdout)

println("initial_force_output: initial force publication adds no stateful queries")
flush(stdout)
@testset "initial force publication adds no stateful queries" begin
    Base.invokelatest(test_stateful_queries)
end
flush(stdout)

println("initial_force_output: both initial savers and first-order solver paths")
flush(stdout)
@testset "both initial savers and first-order solver paths" begin
    Base.invokelatest(test_solver_paths)
end
flush(stdout)
end
