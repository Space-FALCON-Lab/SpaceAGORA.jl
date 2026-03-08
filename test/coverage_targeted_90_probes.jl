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
