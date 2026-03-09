using Test
using StaticArrays

const _TARGET_CALLBACKS = SimulationModel.SimulationCallbacks
const _TARGET_ENV = SimulationModel.EnvironmentModels

struct CoverageBatchDensityModel <: SimulationModel.AbstractDensityModel
end

function SimulationModel.EnvironmentModels.getDensity(
    ::CoverageBatchDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return h * 1e-9, 180.0 + lat + lon + 0.0 * el_time, SVector{3, Float64}(lat, lon, 0.0)
end

struct CoverageForceEffector <: SimulationModel.AbstractForceTorqueModel
    force::SVector{3, Float64}
    torque::SVector{3, Float64}
end

function SimulationModel.calcForceTorque(
    model::CoverageForceEffector,
    x,
    p,
    i::Int64
)
    return model.force, model.torque
end

struct CoverageNoGramBase
end

struct CoverageGramOnlyBase
end

struct CoverageGramThrowsBase
end

struct CoverageGramNoTrajectoryBase
end

struct CoverageGramTrajectoryBase
end

struct CoverageBridgeArgs
    cnf
    solution
end

struct CoverageIndexArgs
    values::Dict{Symbol, Any}
end

Base.propertynames(::CoverageNoGramBase, private::Bool=false) = ()

Base.propertynames(::CoverageGramOnlyBase, private::Bool=false) = (:gram,)
Base.getproperty(::CoverageGramOnlyBase, name::Symbol) = name === :gram ? (; id=1) : throw(ErrorException("unsupported property"))

Base.propertynames(::CoverageGramThrowsBase, private::Bool=false) = (:gram, :gram_atmosphere)
function Base.getproperty(::CoverageGramThrowsBase, name::Symbol)
    if name === :gram
        throw(ErrorException("coverage gram throw"))
    elseif name === :gram_atmosphere
        return :atm
    end
    throw(ErrorException("unsupported property"))
end

Base.propertynames(::CoverageGramNoTrajectoryBase, private::Bool=false) = (:gram, :gram_atmosphere)
function Base.getproperty(::CoverageGramNoTrajectoryBase, name::Symbol)
    if name === :gram
        return (;)
    elseif name === :gram_atmosphere
        return :atm
    end
    throw(ErrorException("unsupported property"))
end

function _coverage_generate_trajectory(
    gram_atmosphere;
    initial_height,
    initial_latitude,
    initial_longitude,
    initial_elapsed_time,
    delta_height,
    delta_latitude,
    delta_longitude,
    delta_elapsed_time,
    n_points,
    update_initial_perturbations=true
)
    step_count = max(2, n_points)
    return [(
        position=(
            height=initial_height + (k - 1) * delta_height,
            latitude=initial_latitude + (k - 1) * delta_latitude,
            longitude=initial_longitude + (k - 1) * delta_longitude,
            elapsedTime=initial_elapsed_time + (k - 1) * delta_elapsed_time
        ),
        dynamics=(density=1.23e-6, temperature=190.0),
        winds=(perturbedEWWind=0.0, perturbedNSWind=0.0, perturbedVerticalWind=0.0)
    ) for k in 1:step_count]
end

Base.propertynames(::CoverageGramTrajectoryBase, private::Bool=false) = (:planet_name, :gram, :gram_atmosphere)
function Base.getproperty(::CoverageGramTrajectoryBase, name::Symbol)
    if name === :planet_name
        return "earth"
    elseif name === :gram
        return (; generate_trajectory=_coverage_generate_trajectory)
    elseif name === :gram_atmosphere
        return :atm
    end
    throw(ErrorException("unsupported property"))
end

function SimulationModel.EnvironmentModels._gram_point_density(
    ::CoverageGramTrajectoryBase,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return 1.23e-6, 190.0, SVector{3, Float64}(0.0, 0.0, 0.0)
end

Base.getindex(args::CoverageIndexArgs, name::Symbol) = args.values[name]

@testset "Coverage Targeted >=90 Probes" begin
    @testset "from_env override branch" begin
        prev_overrides = SimulationEngine._engine_active_overrides_ref[]
        try
            SimulationEngine._engine_active_overrides_ref[] = Dict{String, String}("SPACEAGORA_UNIT_TEST_KEY" => "override_value")
            @test SimulationEngine._engine_env_get("SPACEAGORA_UNIT_TEST_KEY", "fallback") == "override_value"
            @test SimulationEngine._engine_env_get("SPACEAGORA_MISSING_KEY", "fallback") == "fallback"
        finally
            SimulationEngine._engine_active_overrides_ref[] = prev_overrides
        end
    end

    @testset "from_env parser and restore branches" begin
        @test SimulationEngine._env_bool(true) == "1"
        @test SimulationEngine._env_bool(false) == "0"
        @test SimulationEngine._parse_bool(nothing, true) === true
        @test SimulationEngine._parse_bool(" yes ", false) === true
        @test SimulationEngine._parse_bool(" OFF ", true) === false
        @test SimulationEngine._parse_bool("maybe", false) === false

        env_invalid = Dict{String, String}(
            "SPACEAGORA_PARALLEL_PROFILE" => "R7",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "yes",
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "off",
            "SPACEAGORA_EFFECTOR_PARALLEL" => "on",
            "SPACEAGORA_RHS_BATCH_PARALLEL" => "off",
            "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "serial",
            "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "parallel",
            "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => "auto",
            "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
            "SPACEAGORA_SOLVER_MAXITERS" => "not_an_int",
            "SPACEAGORA_SPLIT_IMEX_SOLVER" => "kencarp58",
            "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "not_a_float",
            "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "not_an_int",
            "SPACEAGORA_MULTIRATE_SLOW_DT_S" => "   ",
            "SPACEAGORA_MULTIRATE_SLOW_SOLVER" => "vern9",
            "SPACEAGORA_MULTIRATE_FAST_SOLVER" => "rodas5p",
            "SPACEAGORA_WARN_NORMALIZE" => "maybe",
            "SPACEAGORA_ALLOW_TYPED_NORMALIZE" => "1",
            "SPACEAGORA_GRAM_PER_SAT_INSTANCES" => "0",
            "SPACEAGORA_SRP_EPHEMERIS_CACHE" => "false",
            "SPACEAGORA_NBODY_EPHEMERIS_CACHE" => "true",
            "SPACEAGORA_PLANET_FRAME_CACHE" => "no",
            "SPACEAGORA_SPICE_RHS_MEMO" => "yes",
            "SPACEAGORA_SAVE_BUNDLE" => "off",
            "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "on"
        )
        cfg_invalid = SimulationEngine.simulation_engine_config_from_env(env_invalid)
        @test cfg_invalid.parallel.profile == "R7"
        @test cfg_invalid.parallel.outer_parallel_active === true
        @test cfg_invalid.parallel.parallel_policy_adaptive === false
        @test cfg_invalid.parallel.effector_parallel_mode == "on"
        @test cfg_invalid.solver.mode == "gravity_backbone_split"
        @test cfg_invalid.solver.maxiters === nothing
        @test cfg_invalid.solver.split_imex_solver == "kencarp58"
        @test cfg_invalid.solver.gravity_backbone_dt_s === nothing
        @test cfg_invalid.solver.multirate_fast_substeps == 8
        @test cfg_invalid.solver.multirate_slow_dt_s === nothing
        @test cfg_invalid.solver.multirate_slow_solver == "vern9"
        @test cfg_invalid.solver.multirate_fast_solver == "rodas5p"
        @test cfg_invalid.runtime_policy.warn_normalize === true
        @test cfg_invalid.runtime_policy.allow_typed_normalize === true
        @test cfg_invalid.runtime_policy.gram_per_sat_instances === false
        @test cfg_invalid.runtime_policy.srp_ephemeris_cache === false
        @test cfg_invalid.runtime_policy.nbody_ephemeris_cache === true
        @test cfg_invalid.runtime_policy.planet_frame_cache === false
        @test cfg_invalid.runtime_policy.spice_rhs_memo === true
        @test cfg_invalid.artifacts.save_bundle === false
        @test cfg_invalid.artifacts.warn_deprecated_config === true

        env_valid = Dict{String, String}(
            "SPACEAGORA_SOLVER_MAXITERS" => "321",
            "SPACEAGORA_MULTIRATE_SLOW_DT_S" => "17.5",
            "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "4.25",
            "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "12"
        )
        cfg_valid = SimulationEngine.simulation_engine_config_from_env(env_valid)
        @test cfg_valid.solver.maxiters == 321
        @test cfg_valid.solver.multirate_slow_dt_s == 17.5
        @test cfg_valid.solver.gravity_backbone_dt_s == 4.25
        @test cfg_valid.solver.multirate_fast_substeps == 12

        withenv("SPACEAGORA_UNIT_TEST_ENV_KEY" => "present") do
            prev_overrides = SimulationEngine._engine_active_overrides_ref[]
            try
                SimulationEngine._engine_active_overrides_ref[] = nothing
                @test SimulationEngine._engine_env_get("SPACEAGORA_UNIT_TEST_ENV_KEY", "fallback") == "present"
                @test SimulationEngine._engine_env_haskey("SPACEAGORA_UNIT_TEST_ENV_KEY") === true
                @test SimulationEngine._engine_env_haskey("SPACEAGORA_UNIT_TEST_MISSING_KEY") === false
            finally
                SimulationEngine._engine_active_overrides_ref[] = prev_overrides
            end
        end

        config = SimulationEngine.SimulationEngineConfig(
            parallel=SimulationEngine.ParallelConfig(
                profile="R9",
                outer_parallel_active=true,
                parallel_policy_adaptive=true,
                effector_parallel_mode="on",
                rhs_batch_parallel_mode="off",
                density_callback_parallel_mode="serial",
                control_callback_parallel_mode="parallel",
                thermal_callback_parallel_mode="manual"
            ),
            solver=SimulationEngine.SolverConfig(
                mode="gravity_backbone_split",
                maxiters=123,
                split_imex_solver="kencarp47",
                gravity_backbone_dt_s=6.5,
                multirate_fast_substeps=11,
                multirate_slow_dt_s=14.0,
                multirate_slow_solver="vern9",
                multirate_fast_solver="rodas5p"
            ),
            runtime_policy=SimulationEngine.RuntimePolicyConfig(
                warn_normalize=false,
                allow_typed_normalize=true,
                gram_per_sat_instances=true,
                srp_ephemeris_cache=false,
                nbody_ephemeris_cache=false,
                planet_frame_cache=false,
                spice_rhs_memo=false
            ),
            artifacts=SimulationEngine.ArtifactConfig(
                save_bundle=false,
                warn_deprecated_config=false
            ),
            env_overrides=Dict(
                "SPACEAGORA_SOLVER_MODE" => "from_env_override",
                "SPACEAGORA_CUSTOM_TEST_KEY" => "custom"
            )
        )
        overrides = SimulationEngine._engine_env_overrides(config)
        @test overrides["SPACEAGORA_PARALLEL_PROFILE"] == "R9"
        @test overrides["SPACEAGORA_OUTER_PARALLEL_ACTIVE"] == "1"
        @test overrides["SPACEAGORA_PARALLEL_POLICY_ADAPTIVE"] == "1"
        @test overrides["SPACEAGORA_WARN_NORMALIZE"] == "0"
        @test overrides["SPACEAGORA_ALLOW_TYPED_NORMALIZE"] == "1"
        @test overrides["SPACEAGORA_GRAVITY_BACKBONE_DT_S"] == "6.5"
        @test overrides["SPACEAGORA_MULTIRATE_SLOW_DT_S"] == "14.0"
        @test overrides["SPACEAGORA_SOLVER_MAXITERS"] == "123"
        @test overrides["SPACEAGORA_SOLVER_MODE"] == "from_env_override"
        @test overrides["SPACEAGORA_CUSTOM_TEST_KEY"] == "custom"

        withenv(
            "SPACEAGORA_SOLVER_MODE" => "previous_mode",
            "SPACEAGORA_CUSTOM_TEST_KEY" => nothing
        ) do
            prev_config = SimulationEngine._engine_active_config_ref[]
            prev_overrides = SimulationEngine._engine_active_overrides_ref[]
            result = SimulationEngine._with_engine_env_overrides(() -> begin
                @test SimulationEngine._engine_active_config_ref[] === config
                @test SimulationEngine._engine_active_overrides_ref[] isa Dict{String, String}
                @test ENV["SPACEAGORA_SOLVER_MODE"] == "from_env_override"
                @test ENV["SPACEAGORA_CUSTOM_TEST_KEY"] == "custom"
                @test SimulationEngine._engine_env_get("SPACEAGORA_SOLVER_MODE") == "from_env_override"
                return :ok
            end, config)
            @test result === :ok
            @test ENV["SPACEAGORA_SOLVER_MODE"] == "previous_mode"
            @test !haskey(ENV, "SPACEAGORA_CUSTOM_TEST_KEY")
            @test SimulationEngine._engine_active_config_ref[] === prev_config
            @test SimulationEngine._engine_active_overrides_ref[] === prev_overrides
        end

        withenv("SPACEAGORA_SAVE_BUNDLE" => "1") do
            prev_config = SimulationEngine._engine_active_config_ref[]
            prev_overrides = SimulationEngine._engine_active_overrides_ref[]
            @test_throws ErrorException SimulationEngine._with_engine_env_overrides(config) do
                @test ENV["SPACEAGORA_SAVE_BUNDLE"] == "0"
                error("from_env coverage throw")
            end
            @test ENV["SPACEAGORA_SAVE_BUNDLE"] == "1"
            @test SimulationEngine._engine_active_config_ref[] === prev_config
            @test SimulationEngine._engine_active_overrides_ref[] === prev_overrides
        end
    end

    @testset "thermal contact branch" begin
        model = MaxwellianHeat(
            thermal_accomodation_factor=1.0,
            planet=EARTH,
            thermal_contact=true
        )
        heat = SimulationModel.getHeatRate(model, 2.0, 250.0, 1e-5, 7_500.0, 0.2)
        @test isfinite(heat)
        @test heat > 0.0
    end

    @testset "density model batch branches" begin
        args = build_config(
            spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=120.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),)
        )
        p = ODEParams{1}(args=args)

        hs = [120_000.0, 150_000.0]
        lats = [0.1, 0.2]
        lons = [0.3, 0.4]
        rhos = zeros(2)
        Ts = zeros(2)
        winds = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:2]

        poly_empty = _TARGET_ENV.PolynomialFitAtmosphereModel(Float64[])
        _TARGET_ENV.getDensityBatch!(rhos, Ts, winds, poly_empty, hs, lats, lons, 0.0, true, p)
        @test all(rhos .== 1.0)
        @test all(Ts .== args.environment_model.planet.T_ref)

        generic_model = CoverageBatchDensityModel()
        _TARGET_ENV.getDensityBatch!(rhos, Ts, winds, generic_model, hs, lats, lons, 0.0, true, p)
        @test all(rhos .> 0.0)
        @test all(isfinite.(Ts))
        @test all(winds .== [SVector{3, Float64}(lats[i], lons[i], 0.0) for i in eachindex(hs)])
    end

    @testset "dynamics rhs targeted branches" begin
        f = MVector{3, Float64}(0.0, 0.0, 0.0)
        τ = MVector{3, Float64}(0.0, 0.0, 0.0)
        effectors = (
            CoverageForceEffector(SVector{3, Float64}(1.0, 2.0, 3.0), SVector{3, Float64}(0.1, 0.2, 0.3)),
            CoverageForceEffector(SVector{3, Float64}(2.0, 3.0, 4.0), SVector{3, Float64}(0.4, 0.5, 0.6)),
        )
        decision = (use_threads=true, allotment=2, policy_applied=false, mode=:on)
        SimulationEngine._accumulate_dynamic_effectors!(f, τ, ComponentVector((x=0.0,)), nothing, 1, effectors, decision)
        @test f == MVector{3, Float64}(3.0, 5.0, 7.0)
        @test isapprox(norm(τ - MVector{3, Float64}(0.5, 0.7, 0.9)), 0.0; atol=1e-12, rtol=0.0)

        du_view = ComponentVector((q = zeros(4), ω = zeros(3)))
        sc_view = ComponentVector((q = [1.0, 0.0, 0.0, 0.0], ω = [0.01, -0.02, 0.03]))
        inertia = Matrix{Float64}(I, 3, 3)
        SimulationEngine._assign_orientation_rhs!(
            du_view,
            sc_view,
            inertia,
            SVector{3, Float64}(0.0, 0.0, 0.0);
            propagate_quaternion=false,
            include_gyroscopic=false
        )
        @test all(du_view.q .== 0.0)
        @test all(isfinite.(du_view.ω))
    end

    @testset "density callback helper branches" begin
        withenv("SPACEAGORA_DENSITY_BATCH_PARALLEL" => "off") do
            @test _TARGET_CALLBACKS._density_batch_enabled(8) == false
        end
        withenv("SPACEAGORA_DENSITY_BATCH_PARALLEL" => "on") do
            @test _TARGET_CALLBACKS._density_batch_enabled(1) == true
        end

        surrogate_no_gram = _TARGET_ENV.GRAMAtmosphereModelSurrogate(CoverageNoGramBase(), "unused", nothing)
        surrogate_gram_only = _TARGET_ENV.GRAMAtmosphereModelSurrogate(CoverageGramOnlyBase(), "unused", nothing)
        surrogate_gram_throw = _TARGET_ENV.GRAMAtmosphereModelSurrogate(CoverageGramThrowsBase(), "unused", nothing)
        surrogate_no_traj = _TARGET_ENV.GRAMAtmosphereModelSurrogate(CoverageGramNoTrajectoryBase(), "unused", nothing)
        surrogate_traj = _TARGET_ENV.GRAMAtmosphereModelSurrogate(CoverageGramTrajectoryBase(), "unused", nothing)
        @test _TARGET_CALLBACKS._gram_track_trajectory_supported(surrogate_no_gram) == false
        @test _TARGET_CALLBACKS._gram_track_trajectory_supported(surrogate_gram_only) == false
        @test _TARGET_CALLBACKS._gram_track_trajectory_supported(surrogate_gram_throw) == false
        @test _TARGET_CALLBACKS._gram_track_trajectory_supported(surrogate_no_traj) == false
        @test _TARGET_CALLBACKS._gram_track_trajectory_supported(surrogate_traj) == true

        args_batch = build_config_multi(
            spacecraft=[
                make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=170.0),
                make_spacecraft(ra_alt_m=520e3, rp_alt_m=520e3, ν_deg=165.0),
            ],
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=60.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),)
        )
        p_batch = ODEParams{2}(args=args_batch)
        empty!(p_batch.shared_buffers.density_models)
        push!(p_batch.shared_buffers.density_models, surrogate_no_traj)
        @test _TARGET_CALLBACKS._density_batch_model_for_callback(p_batch, 2) === nothing
        empty!(p_batch.shared_buffers.density_models)
        push!(p_batch.shared_buffers.density_models, surrogate_no_traj)
        push!(p_batch.shared_buffers.density_models, surrogate_traj)
        @test _TARGET_CALLBACKS._density_batch_model_for_callback(p_batch, 2) === nothing
        empty!(p_batch.shared_buffers.density_models)
        common_model = surrogate_traj
        push!(p_batch.shared_buffers.density_models, common_model)
        push!(p_batch.shared_buffers.density_models, common_model)
        @test _TARGET_CALLBACKS._density_batch_model_for_callback(p_batch, 2) === common_model

        withenv("SPACEAGORA_ENTRY_TARGET_COUNT" => "bad") do
            @test_throws ArgumentError _TARGET_CALLBACKS._entry_target_count()
        end

        args_entry = build_config(
            spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
            density_model=_TARGET_ENV.ExponentialAtmosphereModel(1.225, 0.0, 8_500.0),
            orientation_sim=false,
            mission_time=120.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),)
        )
        withenv("SPACEAGORA_ENTRY_TARGET_COUNT" => "1") do
            cbs = _TARGET_CALLBACKS.get_callbacks(1, args_entry.dynamics_model.dynamic_effectors, args_entry)
            @test cbs isa CallbackSet
        end

        gram_model = _TARGET_ENV.GRAMAtmosphereModel(planet_name="earth")
        args_gram = build_config(
            spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
            density_model=gram_model,
            orientation_sim=false,
            mission_time=120.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=false
        )
        p_gram = ODEParams{1}(args=args_gram)
        withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "on") do
            empty!(p_gram.shared_buffers.density_models)
            @test _TARGET_CALLBACKS._gram_isolated_pool_batch_model_for_callback(p_gram, 1) === gram_model
        end

        models_0, locks_0 = _TARGET_CALLBACKS._ensure_gram_isolated_pool!(p_gram, gram_model, 0)
        @test models_0 === p_gram.shared_buffers.gram_isolated_pool_models
        @test locks_0 === p_gram.shared_buffers.gram_isolated_pool_locks

        empty!(p_gram.shared_buffers.gram_isolated_pool_models)
        empty!(p_gram.shared_buffers.gram_isolated_pool_locks)
        push!(p_gram.shared_buffers.gram_isolated_pool_models, gram_model)
        push!(p_gram.shared_buffers.gram_isolated_pool_locks, ReentrantLock())
        push!(p_gram.shared_buffers.gram_isolated_pool_locks, ReentrantLock())
        _TARGET_CALLBACKS._ensure_gram_isolated_pool!(p_gram, gram_model, 1)
        @test length(p_gram.shared_buffers.gram_isolated_pool_models) == 1
        @test p_gram.shared_buffers.gram_isolated_pool_models[1] isa _TARGET_ENV.GRAMAtmosphereModel

        ρ_poly, T_poly, wind_poly = _TARGET_CALLBACKS._gram_isolated_pool_density_state(
            gram_model,
            500_000.0,
            0.0,
            0.0,
            0.0,
            true,
            p_gram,
            ReentrantLock()
        )
        @test isfinite(ρ_poly)
        @test isfinite(T_poly)
        @test wind_poly isa SVector{3, Float64}

        ρ_vac, T_vac, wind_vac = _TARGET_CALLBACKS._gram_isolated_pool_density_state(
            gram_model,
            2_100_000.0,
            0.0,
            0.0,
            0.0,
            true,
            p_gram,
            ReentrantLock()
        )
        @test ρ_vac == 0.0
        @test T_vac == p_gram.args.environment_model.planet.T_ref
        @test wind_vac == SVector{3, Float64}(0.0, 0.0, 0.0)

        args_cache = build_config(
            spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
            density_model=surrogate_traj,
            orientation_sim=false,
            mission_time=120.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),)
        )
        p_cache = ODEParams{1}(args=args_cache)
        u_cache = build_initial_conditions(args_cache)
        pos = SVector{3, Float64}(u_cache.sc[1].pos)
        vel = SVector{3, Float64}(u_cache.sc[1].vel)
        rp, _ = r_intor_p!(pos, vel, p_cache.args.environment_model.planet)
        alt, lat, lon = rtolatlong(rp, p_cache.args.environment_model.planet)
        cache = _TARGET_CALLBACKS.GramTrackCache()
        cache.valid = true
        cache.t0 = 0.0
        cache.t1 = 10.0
        cache.index_hint = 1
        cache.times = [0.0, 10.0]
        cache.alts = [alt, alt]
        cache.lats = [lat, lat]
        cache.lons = [lon, lon]
        cache.rhos = [1.23e-6, 1.23e-6]
        cache.Ts = [190.0, 190.0]
        cache.winds = [SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)]
        p_cache.shared_buffers.gram_density_cache[1] = cache
        integrator = (p=p_cache, u=u_cache, t=0.0, sol=(prob=(tspan=(0.0, 120.0),),))
        withenv(
            "SPACEAGORA_GRAM_PROFILE" => "1",
            "SPACEAGORA_GRAM_TRACK_CACHE" => "on",
            "SPACEAGORA_DENSITY_BATCH_PARALLEL" => "off",
            "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off"
        ) do
            density_cb = _TARGET_CALLBACKS.get_density_callback(1, args_cache.dynamics_model.dynamic_effectors, args_cache)
            density_cb.affect!(integrator)
        end
        @test isapprox(p_cache.shared_buffers.densities[1], 1.23e-6; atol=1e-12, rtol=0.0)
        @test isapprox(p_cache.shared_buffers.temperatures[1], 190.0; atol=1e-12, rtol=0.0)
    end
end

@testset "GNC Bridge Helper Probes" begin
    gh = SimulationModel.GuidanceHooks
    ch = SimulationModel.ControlHooks

    cnf = (time_switch_1 = 1.0, time_switch_2 = 2.0)
    solution = (orientation = (; time = Float64[]),)
    args = (; cnf = cnf, solution = solution, EI = 120.0, dry_mass = 1.0)
    object_args = CoverageBridgeArgs(cnf, solution)

    @test gh._bridge_get_cnf(args) === cnf
    @test gh._bridge_get_solution(args; cnf=cnf) === solution
    @test ch._bridge_get_cnf(args) === cnf
    @test ch._bridge_get_solution(args; cnf=cnf) === solution
    @test gh._bridge_get_cnf(nothing; cnf=cnf) === cnf
    @test gh._bridge_get_solution(nothing; solution=solution) === solution
    @test gh._bridge_get_cnf(object_args) === cnf
    @test gh._bridge_get_solution(object_args) === solution
    @test gh._bridge_get_solution((;); cnf=(solution=solution,)) === solution
    runtime_context = gh._make_aerobraking_runtime_context(
        mission=:mission,
        index_phase_aerobraking=2,
        ip=:ip,
        aerobraking_phase=3,
        t_prev=4.0,
        date_initial=:date,
        time_0=5.0,
        args=args,
        initial_state=:state,
        gram_atmosphere=:atm,
        gram=:gram,
        cnf=cnf,
        solution=solution,
    )
    @test runtime_context.cnf === cnf
    @test gh._with_control_gain(runtime_context, 1.5).control_gain == 1.5
    @test gh._with_time_switch(runtime_context, 7.5).time_switch == 7.5
    @test gh.CONTROL_BRIDGE_STATE_LOCK isa ReentrantLock
    @test ch.CONTROL_BRIDGE_STATE_LOCK isa ReentrantLock

    spacecraft_args = (;
        simulation_settings=(verbose=true,),
        environment_model=(topography=true, EI=111.0),
        dynamics_model=(spacecraft=[(; dry_mass=66.0)],),
    )
    mission = (
        aerodynamics=(heat_rate_limit=55.0,),
        body=(dry_mass=77.0,),
        engines=(ϕ=0.75,),
    )
    dict_args = Dict{Symbol, Any}(
        :EI => 120.0,
        :AE => 135.0,
        :body_shape => "Capsule",
        :heat_load_sol => 2,
        :max_heat_rate => 88.0,
        :srp => true,
        :control_mode => 3,
        :struct_ctrl => true,
        :dry_mass => 44.0,
        :phi => 0.25,
        :control_in_loop => true,
        :integrator => "Custom",
        :drag_passage => true,
        :topography_model => "Spherical Harmonics",
    )
    index_args = CoverageIndexArgs(Dict{Symbol, Any}(:AE => 146.0))

    withenv("SPACEAGORA_DEBUG_LEGACY_CONTROL" => "1") do
        @test gh._bridge_verbose_enabled() === true
    end
    @test gh._bridge_verbose_enabled(spacecraft_args) === true
    @test gh._bridge_verbose_enabled((; verbose=true)) === true
    @test gh._bridge_required_field(dict_args, :EI) == 120.0
    @test gh._bridge_optional_field(dict_args, :missing, 9.0) == 9.0
    @test gh._bridge_optional_field(index_args, :AE, 0.0) == 146.0
    @test gh._bridge_aerobraking_topography_enabled(spacecraft_args) === true
    @test gh._bridge_aerobraking_topography_enabled(dict_args) === true
    @test gh._bridge_aerobraking_entry_interface_m(spacecraft_args) == 111_000.0
    @test gh._bridge_aerobraking_entry_interface_m(dict_args) == 120_000.0
    @test gh._bridge_aerobraking_exit_interface_m((; EI=120.0)) == 120_000.0
    @test gh._bridge_aerobraking_exit_interface_m(dict_args) == 135_000.0
    @test gh._bridge_aerobraking_body_shape(dict_args) == "Capsule"
    @test gh._bridge_aerobraking_heat_load_solution(dict_args) == 2
    @test gh._bridge_aerobraking_max_heat_rate((;), mission) == 55.0
    @test gh._bridge_aerobraking_max_heat_rate(dict_args, mission) == 88.0
    @test gh._bridge_aerobraking_srp_enabled(dict_args) === true
    @test gh._bridge_aerobraking_control_mode(dict_args) == 3
    @test gh._bridge_aerobraking_struct_control_enabled(dict_args) === true
    @test gh._bridge_aerobraking_dry_mass(spacecraft_args, mission) == 66.0
    @test gh._bridge_aerobraking_dry_mass((;), mission) == 77.0
    @test gh._bridge_aerobraking_dry_mass(dict_args, mission) == 77.0
    @test gh._bridge_aerobraking_dry_mass(dict_args, nothing) == 44.0
    @test gh._bridge_aerobraking_thrust_phi((;), mission) == 0.75
    @test gh._bridge_aerobraking_thrust_phi(dict_args, mission) == 0.25
    @test gh._bridge_aerobraking_control_in_loop(dict_args) === true
    @test gh._bridge_aerobraking_integrator_name(dict_args) == "Custom"
    @test gh._bridge_aerobraking_drag_passage(dict_args) === true

    settings = gh._make_aerobraking_runtime_settings(dict_args, mission)
    @test settings.topography_enabled === true
    @test settings.entry_interface_m == 120_000.0
    @test settings.exit_interface_m == 135_000.0
    @test settings.dry_mass == 77.0
    @test settings.thrust_phi == 0.25

    runtime_context_mission = gh._make_aerobraking_runtime_context(
        mission=mission,
        index_phase_aerobraking=1,
        ip=:ip,
        aerobraking_phase=2,
        t_prev=3.0,
        date_initial=:date,
        time_0=4.0,
        args=dict_args,
        initial_state=:state,
        gram_atmosphere=:atm,
        gram=:gram,
    )
    @test runtime_context_mission.settings.max_heat_rate == 88.0
    @test runtime_context_mission.settings.control_mode == 3

    @test_throws ArgumentError gh._bridge_get_cnf((;))
    @test_throws ArgumentError gh._bridge_get_solution((;); cnf=nothing)
    @test_throws ArgumentError ch._bridge_get_cnf((;))
    @test_throws ArgumentError ch._bridge_get_solution((;); cnf=nothing)
end

@testset "Point Mass Dynamics Branch Probes" begin
    translational = SimulationModel.DynamicsTranslational

    zero_acc = translational.acceleration_from_force(SVector(1.0, -2.0, 3.0), 0.0)
    @test zero_acc == SVector{3, Float64}(0.0, 0.0, 0.0)

    inf_acc = translational.acceleration_from_force(SVector(1.0, -2.0, 3.0), Inf)
    @test inf_acc == SVector{3, Float64}(0.0, 0.0, 0.0)

    @test translational.mass_derivative(NaN) == 0.0
end

@testset "Low-Coverage Utility Branch Probes" begin
    @test SimulationModel.NoGramPresets.make_no_gram_planet(:venus) isa Venus
    @test_throws ArgumentError SimulationModel.NoGramPresets.make_no_gram_planet(:pluto)
    @test_throws ArgumentError SimulationModel.NoGramPresets.make_no_gram_density_model(Earth(), :mystery)

    @test TelemetryVerification._planet_from_name("mars") isa Mars
    @test TelemetryVerification._planet_from_name("venus") isa Venus
    @test TelemetryVerification._nbody_primary_name("mars") == "Mars"
    @test TelemetryVerification._nbody_primary_name("venus") == "Venus"
    @test TelemetryVerification._nbody_primary_name("titan") == "Titan"

    args_public = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=5.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        ephemerides_model=SimpleEphemeridesModel()
    )
    @test SimulationEngine._require_simulation_configuration(args_public) === args_public
    invoke(SimulationEngine.run_simulation, Tuple{Any}, args_public)

    cli = SpaceAGORA.SpaceAGORACLI
    @test cli._normalize_example_name("AGORA_Earth_NoGRAM") == "AGORA_Earth_NoGRAM.jl"
    temp_example, temp_io = mktemp()
    close(temp_io)
    @test cli._resolve_example_path(temp_example) == abspath(temp_example)

    script_path, script_io = mktemp()
    write(script_io, "println(\"cli_probe_ok\")\n")
    close(script_io)
    out = IOBuffer()
    err = IOBuffer()
    @test cli._run_subprocess(script_path, String[]; io=out, errio=err) == 0
    @test occursin("cli_probe_ok", String(take!(out)))

    assets_text = sprint(io -> @test SpaceAGORA.run_cli(["assets", "check"]; io=io, errio=io) == 0)
    @test occursin("no_gram_mode", assets_text)
    @test_throws ArgumentError SpaceAGORA.run_cli(["run", "--bad-flag"])
    @test_throws ArgumentError SpaceAGORA.run_cli(["telemetry", "--bogus"])
    @test_throws ArgumentError SpaceAGORA.run_cli(["benchmark", "bogus-mode"])
    @test_throws ArgumentError SpaceAGORA.run_cli(["assets", "bogus"])
    @test_throws ArgumentError SpaceAGORA.run_cli(["bogus"])
end

@testset "Precompile Workload Probe" begin
    precompile_probe = Module(:PrecompileCoverageProbe)
    Core.eval(precompile_probe, quote
        const SimulationModel = Main.SimulationModel
        const TelemetryVerification = Main.TelemetryVerification
        const parse_parallel_profile = Main.SpaceAGORA.parse_parallel_profile
        const simulation_engine_config_from_env = Main.SimulationEngine.simulation_engine_config_from_env
        const run_simulation = Main.SimulationEngine.run_simulation

        macro setup_workload(ex)
            return esc(ex)
        end

        macro compile_workload(ex)
            return esc(ex)
        end

        Base.include(@__MODULE__, joinpath(Main.REPO_ROOT, "src", "precompile_workload.jl"))
    end)

    @test isdefined(precompile_probe, :_run_spaceagora_precompile_workload)
    @test isdefined(precompile_probe, :_spaceagora_precompile_args)
end

@testset "Effector Sampling Helper Branch Probes" begin
    pos_vec, vel_vec = SimulationEngine._extract_sample_pos_vel([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    @test pos_vec == SVector{3, Float64}(1.0, 2.0, 3.0)
    @test vel_vec == SVector{3, Float64}(4.0, 5.0, 6.0)
    @test SimulationEngine._extract_sample_mass_kg((mass_kg=5.0,)) == 5.0
    @test SimulationEngine._extract_sample_mass_kg([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]) == 7.0

    args_ephem = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        ephemerides_model=SpiceEphemeridesModel()
    )
    p_ephem = ODEParams{1}(args=args_ephem)
    u_ephem = build_initial_conditions(args_ephem)
    sc_ephem = u_ephem.sc[1]
    sample_ephem = build_state_sample(sc_ephem, args_ephem.dynamics_model.spacecraft[1], false)

    p_ephem.shared_buffers.densities[1] = 9.0
    p_ephem.shared_buffers.temperatures[1] = 250.0
    p_ephem.shared_buffers.winds[1] = SVector{3, Float64}(1.0, 2.0, 3.0)
    buffered = SimulationEngine.sample_buffered_atmosphere(sample_ephem, p_ephem, 1, 0.0)
    @test buffered.rho_kg_m3 == 9.0
    @test buffered.temperature_k == 250.0
    @test buffered.wind_pp == SVector{3, Float64}(1.0, 2.0, 3.0)

    p_ephem.shared_buffers.srp_sun_ephemeris_cache[] = nothing
    solar = SimulationEngine.sample_solar_ephemeris(sample_ephem, p_ephem, 1, 0.0)
    @test all(isfinite, solar.sun_pos_ii)

    p_ephem.shared_buffers.nbody_ephemeris_cache[] = nothing
    nbody_model = NBodyGravityModel(body_names=("Moon",), primary_body_name="Earth", planet=EARTH)
    third_body = SimulationEngine.sample_third_body_ephemerides(nbody_model, sample_ephem, p_ephem, 1, 0.0)
    @test third_body.names == ("Moon",)
    @test length(third_body.positions_ii) == 1
    @test all(isfinite, third_body.positions_ii[1])
end

@testset "Density Runtime Helper Branch Probes" begin
    args_density = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=520e3, rp_alt_m=520e3, ν_deg=165.0),
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        ephemerides_model=SimpleEphemeridesModel(),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_density = ODEParams{2}(args=args_density)
    u_density = build_initial_conditions(args_density)
    cache_cfg = _TARGET_CALLBACKS._gram_track_cache_config()

    p_density.shared_buffers.solve_segment_end_time[] = 25.0
    @test _TARGET_CALLBACKS._density_segment_end_t(p_density, 5.0, cache_cfg) == 25.0
    p_density.shared_buffers.solve_segment_end_time[] = 0.0
    @test _TARGET_CALLBACKS._density_segment_end_t(p_density, 5.0, cache_cfg) == 60.0
    p_density_inf = (
        shared_buffers=(solve_segment_end_time=Ref(0.0),),
        args=(mission_configuration=(mission_time=Inf,),),
    )
    @test _TARGET_CALLBACKS._density_segment_end_t(p_density_inf, 5.0, cache_cfg) == 5.0 + cache_cfg.orbit_horizon_s

    p_density.shared_buffers.densities[1] = 1.5
    p_density.shared_buffers.temperatures[1] = 222.0
    p_density.shared_buffers.winds[1] = SVector{3, Float64}(4.0, 5.0, 6.0)
    buffered_state = _TARGET_CALLBACKS._buffered_stage_environment_state(u_density.sc[1], p_density, 1, 0.0)
    @test buffered_state.rho == 1.5
    @test buffered_state.T == 222.0
    @test buffered_state.wind == SVector{3, Float64}(4.0, 5.0, 6.0)
    @test isfinite(buffered_state.alt)

    _TARGET_CALLBACKS._write_density_buffers!(p_density, 1, 2.5, 333.0, SVector{3, Float64}(7.0, 8.0, 9.0))
    @test p_density.shared_buffers.densities[1] == 2.5
    @test p_density.shared_buffers.temperatures[1] == 333.0
    @test p_density.shared_buffers.winds[1] == SVector{3, Float64}(7.0, 8.0, 9.0)

    surrogate_density = _TARGET_ENV.GRAMAtmosphereModel(planet_name="earth")
    empty!(p_density.shared_buffers.density_models)
    push!(p_density.shared_buffers.density_models, surrogate_density)
    push!(p_density.shared_buffers.density_models, surrogate_density)
    integrator_density = (p=p_density, u=u_density, t=0.0, sol=(prob=(tspan=(0.0, 60.0),),))
    withenv(
        "SPACEAGORA_DENSITY_BATCH_PARALLEL" => "off",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off"
    ) do
        density_cb = _TARGET_CALLBACKS.get_density_callback(2, args_density.dynamics_model.dynamic_effectors, args_density)
        density_cb.affect!(integrator_density)
    end
    @test all(p_density.shared_buffers.densities .>= 0.0)
    @test all(isfinite.(p_density.shared_buffers.temperatures))

    withenv(
        "SPACEAGORA_DENSITY_BATCH_PARALLEL" => "on",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_GRAM_ISOLATED_POOL" => "off"
    ) do
        density_cb_batch_models = _TARGET_CALLBACKS.get_density_callback(2, args_density.dynamics_model.dynamic_effectors, args_density)
        density_cb_batch_models.affect!(integrator_density)
    end
    @test all(p_density.shared_buffers.densities .>= 0.0)
    @test all(isfinite.(p_density.shared_buffers.temperatures))

    args_density_batch = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=480e3, rp_alt_m=480e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=510e3, rp_alt_m=510e3, ν_deg=165.0),
        ],
        density_model=_TARGET_ENV.GRAMAtmosphereModel(planet_name="earth"),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        ephemerides_model=SimpleEphemeridesModel(),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_density_batch = ODEParams{2}(args=args_density_batch)
    u_density_batch = build_initial_conditions(args_density_batch)
    empty!(p_density_batch.shared_buffers.density_models)
    integrator_density_batch = (p=p_density_batch, u=u_density_batch, t=0.0, sol=(prob=(tspan=(0.0, 60.0),),))
    withenv(
        "SPACEAGORA_DENSITY_BATCH_PARALLEL" => "on",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_GRAM_ISOLATED_POOL" => "on"
    ) do
        density_cb_batch = _TARGET_CALLBACKS.get_density_callback(2, args_density_batch.dynamics_model.dynamic_effectors, args_density_batch)
        density_cb_batch.affect!(integrator_density_batch)
    end
    @test all(p_density_batch.shared_buffers.densities .>= 0.0)
    @test all(isfinite.(p_density_batch.shared_buffers.temperatures))
    @test all(all(isfinite, wind) for wind in p_density_batch.shared_buffers.winds)

    withenv(
        "SPACEAGORA_DENSITY_BATCH_PARALLEL" => "off",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_GRAM_ISOLATED_POOL" => "on"
    ) do
        density_cb_pool = _TARGET_CALLBACKS.get_density_callback(2, args_density_batch.dynamics_model.dynamic_effectors, args_density_batch)
        density_cb_pool.affect!(integrator_density_batch)
    end
    @test all(p_density_batch.shared_buffers.densities .>= 0.0)
    @test all(isfinite.(p_density_batch.shared_buffers.temperatures))
end

@testset "Aerodynamic Helper Branch Probes" begin
    dyn = SimulationModel.DynamicEffectors
    workspace = dyn._make_aero_scratch_workspace(1)
    dyn._ensure_aero_workspace_capacity!(workspace, 3)
    @test length(workspace.thread_force) == 3
    @test length(workspace.thread_cl) == 3
    @test dyn.AerodynamicEffectors._constant_drag_coefficient(0.1) > 0.8

    args_aero = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        ephemerides_model=SimpleEphemeridesModel(),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_aero = ODEParams{1}(args=args_aero)
    fresh_workspace = dyn._aero_workspace_for_sat!(p_aero, 5, 2)
    @test length(fresh_workspace.thread_force) == 2

    state = SimulationModel.StateSample(
        SVector{3, Float64}(EARTH.Rp_e + 400e3, 0.0, 0.0),
        SVector{3, Float64}(0.0, 7_500.0, 0.0),
        500.0;
        spacecraft=args_aero.dynamics_model.spacecraft[1]
    )
    planet_frame = SimulationModel.PlanetFrameSample(
        SMatrix{3, 3, Float64, 9}(I),
        SVector{3, Float64}(EARTH.Rp_e + 400e3, 0.0, 0.0),
        SVector{3, Float64}(0.0, 7_500.0, 0.0),
        400e3,
        0.0,
        0.0
    )

    env_bad_density = SimulationModel.EnvironmentSample(
        EARTH;
        planet_frame=planet_frame,
        atmosphere=SimulationModel.AtmosphereSample(0.0, 200.0, SVector{3, Float64}(0.0, 0.0, 0.0))
    )
    force_bad_density, torque_bad_density = SimulationModel.wrench(AerodynamicCoefficientfM(), state, env_bad_density, 0.0)
    @test force_bad_density == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test torque_bad_density == SVector{3, Float64}(0.0, 0.0, 0.0)

    env_zero_h = SimulationModel.EnvironmentSample(
        EARTH;
        planet_frame=SimulationModel.PlanetFrameSample(
            SMatrix{3, 3, Float64, 9}(I),
            SVector{3, Float64}(EARTH.Rp_e + 400e3, 0.0, 0.0),
            SVector{3, Float64}(1.0, 0.0, 0.0),
            400e3,
            0.0,
            0.0
        ),
        atmosphere=SimulationModel.AtmosphereSample(1e-6, 220.0, SVector{3, Float64}(0.0, 0.0, 0.0))
    )
    force_zero_h, torque_zero_h = SimulationModel.wrench(AerodynamicCoefficientConstant(), state, env_zero_h, 0.0)
    @test force_zero_h == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test torque_zero_h == SVector{3, Float64}(0.0, 0.0, 0.0)

    uD, uN, uE = latlongtoNED((planet_frame.alt_m, planet_frame.lat_rad, planet_frame.lon_rad))
    desired_wind_pp = -planet_frame.vel_pp
    wind_components = SVector{3, Float64}(
        dot(desired_wind_pp, uE),
        dot(desired_wind_pp, uN),
        -dot(desired_wind_pp, uD),
    )
    env_zero_rel = SimulationModel.EnvironmentSample(
        EARTH;
        planet_frame=planet_frame,
        atmosphere=SimulationModel.AtmosphereSample(1e-6, 220.0, wind_components)
    )
    force_zero_rel, torque_zero_rel = SimulationModel.wrench(AerodynamicCoefficientNoBallisticFlight(), state, env_zero_rel, 0.0)
    @test force_zero_rel == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test torque_zero_rel == SVector{3, Float64}(0.0, 0.0, 0.0)

    args_legacy_aero = build_config(
        spacecraft=make_spacecraft(ra_alt_m=450e3, rp_alt_m=420e3, ν_deg=160.0),
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        ephemerides_model=SimpleEphemeridesModel(),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_legacy_aero = ODEParams{1}(args=args_legacy_aero)
    u_legacy_aero = build_initial_conditions(args_legacy_aero)
    p_legacy_aero.shared_buffers.current_time[] = 0.0
    force_legacy, torque_legacy = SimulationModel.calcForceTorque(
        AerodynamicCoefficientfM(),
        u_legacy_aero.sc[1],
        p_legacy_aero,
        1
    )
    @test all(isfinite, force_legacy)
    @test all(isfinite, torque_legacy)
    @test norm(force_legacy) > 0.0
end

@testset "Dynamics RHS Batch and Helper Branch Probes" begin
    du_heat = zeros(2)
    SimulationEngine._assign_heat_rate_derivative!(du_heat, [1.0, 2.0, 3.0])
    @test du_heat == [1.0, 2.0]

    q0 = SVector{4, Float64}(1.0, 0.0, 0.0, 0.0)
    ω0 = SVector{3, Float64}(0.01, -0.02, 0.03)
    sc_rhs = [
        make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=170.0, orientation_state=(q0, ω0)),
        make_spacecraft(ra_alt_m=520e3, rp_alt_m=520e3, ν_deg=165.0, orientation_state=(q0, ω0)),
    ]
    args_rhs = build_config_multi(
        spacecraft=sc_rhs,
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=true,
        mission_time=5.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            AerodynamicCoefficientfM(),
            ConstantTorqueModel(SVector{3, Float64}(0.05, -0.02, 0.01)),
        ),
        control_effectors=(TimedTangentialThrusterModel(0.2, 1.0, 0.0, 5.0),),
        control_rates=[1.0],
        ephemerides_model=SimpleEphemeridesModel(),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    u_rhs = build_initial_conditions(args_rhs)
    p_rhs = ODEParams{2}(args=args_rhs)
    _initialize_heat_rate_buffers!(p_rhs)
    p_rhs.is_active[2] = false

    withenv("SPACEAGORA_RHS_BATCH_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_PARALLEL" => "off") do
        du_full = copy(u_rhs)
        du_slow = copy(u_rhs)
        du_implicit = copy(u_rhs)
        du_explicit = copy(u_rhs)
        du_fast = copy(u_rhs)
        du_full .= 0.0
        du_slow .= 0.0
        du_implicit .= 0.0
        du_explicit .= 0.0
        du_fast .= 0.0

        SimulationEngine.spacecraft_dynamics!(du_full, u_rhs, p_rhs, 0.0)
        SimulationEngine.spacecraft_dynamics_slow!(du_slow, u_rhs, p_rhs, 0.0)
        spacecraft_dynamics_implicit_atmosphere!(du_implicit, u_rhs, p_rhs, 0.0)
        spacecraft_dynamics_explicit_remainder!(du_explicit, u_rhs, p_rhs, 0.0)
        SimulationEngine.spacecraft_dynamics_fast_control!(du_fast, u_rhs, p_rhs, 0.0)

        @test norm(SVector{3, Float64}(du_full.sc[1].vel)) > 0.0
        @test norm(SVector{3, Float64}(du_slow.sc[1].vel)) > 0.0
        @test all(isfinite, du_implicit.sc[1].vel)
        @test norm(SVector{3, Float64}(du_explicit.sc[1].pos)) > 0.0
        @test norm(SVector{3, Float64}(du_fast.sc[1].vel)) > 0.0
        @test du_implicit.sc[1].heat_loads == zeros(length(du_implicit.sc[1].heat_loads))
        @test all(==(0.0), du_full.sc[2].pos)
        @test all(==(0.0), du_explicit.sc[2].pos)
        @test all(==(0.0), du_fast.sc[2].pos)
    end

    args_backbone_batch = build_config_multi(
        spacecraft=[
            make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=170.0),
            make_single_link_spacecraft(ra_alt_m=520e3, rp_alt_m=520e3, ν_deg=165.0),
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=5.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        ephemerides_model=SimpleEphemeridesModel(),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    u_backbone_batch = build_initial_conditions(args_backbone_batch)
    p_backbone_batch = ODEParams{2}(args=args_backbone_batch)
    p_backbone_batch.is_active[2] = false
    q_backbone, dq_backbone = _gravity_backbone_initial_states(u_backbone_batch, args_backbone_batch)
    ddu_backbone = copy(dq_backbone)
    ddu_backbone .= 0.0
    withenv("SPACEAGORA_RHS_BATCH_PARALLEL" => "on") do
        spacecraft_dynamics_gravity_backbone!(ddu_backbone, dq_backbone, q_backbone, p_backbone_batch, 0.0)
    end
    @test norm(SVector{3, Float64}(ddu_backbone.sc[1].vel)) > 0.0
    @test ddu_backbone.sc[2].vel == zeros(length(ddu_backbone.sc[2].vel))
end
