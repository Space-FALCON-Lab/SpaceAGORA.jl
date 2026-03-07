@testset "Solver/Env Helper Parsing Coverage" begin
    withenv("SPACEAGORA_SOLVER_MODE" => nothing) do
        @test _solver_policy_mode() == :tsit5
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "auto") do
        @test _solver_policy_mode() == :auto_stiff
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "rodas") do
        @test _solver_policy_mode() == :rodas5p
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "split_imex") do
        @test _solver_policy_mode() == :split_imex
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "unsupported-mode") do
        @test_throws ArgumentError _solver_policy_mode()
    end

    withenv("SPACEAGORA_GRAM_PER_SAT_INSTANCES" => "on") do
        @test _gram_per_sat_instances_enabled() == true
    end
    withenv("SPACEAGORA_GRAM_PER_SAT_INSTANCES" => "off") do
        @test _gram_per_sat_instances_enabled() == false
    end
    withenv("SPACEAGORA_GRAM_PER_SAT_INSTANCES" => "maybe") do
        @test_throws ArgumentError _gram_per_sat_instances_enabled()
    end

    withenv("SPACEAGORA_EFFECTOR_PARALLEL" => "off") do
        @test _effector_parallel_mode() == :off
    end
    withenv("SPACEAGORA_EFFECTOR_PARALLEL" => "on") do
        @test _effector_parallel_mode() == :on
    end
    withenv("SPACEAGORA_EFFECTOR_PARALLEL" => "auto") do
        @test _effector_parallel_mode() == :auto
    end
    withenv("SPACEAGORA_EFFECTOR_PARALLEL" => "invalid") do
        @test_throws ArgumentError _effector_parallel_mode()
    end
    withenv("SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "3") do
        @test _effector_thread_threshold() == 3
    end
    withenv("SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "oops") do
        @test_throws ArgumentError _effector_thread_threshold()
    end
    withenv("SPACEAGORA_EFFECTOR_MAX_THREADS" => "3") do
        @test _effector_max_threads() == 3
    end
    withenv("SPACEAGORA_EFFECTOR_MAX_THREADS" => "oops") do
        @test_throws ArgumentError _effector_max_threads()
    end
    withenv("SPACEAGORA_EFFECTOR_LONG_MISSION_THRESHOLD_S" => "10.0") do
        @test _effector_long_mission_threshold_s() == 10.0
    end
    withenv("SPACEAGORA_EFFECTOR_LONG_MISSION_THRESHOLD_S" => "0.0") do
        @test_throws ArgumentError _effector_long_mission_threshold_s()
    end
    withenv("SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT" => "12345.0") do
        @test _effector_cost_ns_per_item_default() == 12345.0
    end
    withenv("SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT" => "0.0") do
        @test_throws ArgumentError _effector_cost_ns_per_item_default()
    end
    withenv("SPACEAGORA_EFFECTOR_COST_EMA_ALPHA" => "0.3") do
        @test _effector_cost_ema_alpha() == 0.3
    end
    withenv("SPACEAGORA_EFFECTOR_COST_EMA_ALPHA" => "1.5") do
        @test_throws ArgumentError _effector_cost_ema_alpha()
    end
    withenv("SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD" => "25000.0") do
        @test _effector_work_ns_per_worker_threshold() == 25000.0
    end
    withenv("SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD" => "0.0") do
        @test_throws ArgumentError _effector_work_ns_per_worker_threshold()
    end

    args_eff_single = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), InverseSquaredJ2GravityModel()),
        keplerian=true
    )
    p_eff_single = ODEParams{1}(args=args_eff_single)
    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
        "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
        "SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY" => "1",
        "SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT" => "120000.0",
        "SPACEAGORA_EFFECTOR_COST_MIN_SAMPLES" => "999",
        "SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD" => "1000.0",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0"
    ) do
        decision_single = _dynamic_effector_thread_decision(args_eff_single, p_eff_single, args_eff_single.dynamics_model.dynamic_effectors, 1)
        if Threads.nthreads() > 1
            @test decision_single.use_threads == true
            @test decision_single.allotment >= 2
            @test decision_single.allotment <= min(Threads.nthreads(), 4)
        else
            @test decision_single.use_threads == false
        end
    end

    args_eff_multi = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), InverseSquaredJ2GravityModel()),
        keplerian=true
    )
    p_eff_multi = ODEParams{1}(args=args_eff_multi)
    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
        "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
        "SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY" => "1",
        "SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT" => "1.0",
        "SPACEAGORA_EFFECTOR_COST_MIN_SAMPLES" => "999",
        "SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD" => "1.0e9",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0"
    ) do
        decision_multi = _dynamic_effector_thread_decision(args_eff_multi, p_eff_multi, args_eff_multi.dynamics_model.dynamic_effectors, 1)
        @test decision_multi.use_threads == false
    end

    args_eff_constellation = args_eff_single
    p_eff_constellation = ODEParams{1}(args=args_eff_constellation)
    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
        "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
        "SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY" => "1",
        "SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT" => "120000.0",
        "SPACEAGORA_EFFECTOR_COST_MIN_SAMPLES" => "999",
        "SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD" => "1000.0",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_EFFECTOR_MAX_THREADS" => string(Threads.nthreads())
    ) do
        decision_constellation = _dynamic_effector_thread_decision(
            args_eff_constellation,
            p_eff_constellation,
            args_eff_constellation.dynamics_model.dynamic_effectors,
            4
        )
        share_budget = max(1, fld(max(1, Threads.nthreads()), min(4, max(1, Threads.nthreads()))))
        inner_floor = Threads.nthreads() > 1 ? min(2, Threads.nthreads()) : 1
        expected_cap = min(Threads.nthreads(), max(share_budget, inner_floor))
        if expected_cap > 1
            @test decision_constellation.use_threads == true
            @test decision_constellation.allotment <= expected_cap
        else
            @test decision_constellation.use_threads == false
            @test decision_constellation.allotment == 1
        end
    end

    @test _retcode_is_stiff_symptom(:Unstable)
    @test _retcode_is_stiff_symptom("DtLessThanMin")
    @test _retcode_is_stiff_symptom(:InitialFailure)
    @test !_retcode_is_stiff_symptom(:Success)

    withenv("SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        @test _solver_maxiters() === nothing
    end
    withenv("SPACEAGORA_SOLVER_MAXITERS" => "2500") do
        @test _solver_maxiters() == 2500
    end
    withenv("SPACEAGORA_SOLVER_MAXITERS" => "0") do
        @test_throws ArgumentError _solver_maxiters()
    end
    withenv("SPACEAGORA_SOLVER_MAXITERS" => "not-an-int") do
        @test_throws ArgumentError _solver_maxiters()
    end

    withenv("SPACEAGORA_SPLIT_IMEX_SOLVER" => nothing) do
        split_solver = _split_imex_solver_spec()
        @test split_solver.label == "KenCarp4"
    end
    withenv("SPACEAGORA_SPLIT_IMEX_SOLVER" => "kencarp47") do
        split_solver = _split_imex_solver_spec()
        @test split_solver.label == "KenCarp47"
    end
    withenv("SPACEAGORA_SPLIT_IMEX_SOLVER" => "kencarp58") do
        split_solver = _split_imex_solver_spec()
        @test split_solver.label == "KenCarp58"
    end
    withenv("SPACEAGORA_SPLIT_IMEX_SOLVER" => "unsupported-split-solver") do
        @test_throws ArgumentError _split_imex_solver_spec()
    end

    struct DummyNoAlgChoice end
    struct DummyAlgChoice
        alg_choice::Vector{Int}
    end
    @test _auto_stiff_switched(DummyNoAlgChoice()) == false
    @test _auto_stiff_switched(DummyAlgChoice(Int[])) == false
    @test _auto_stiff_switched(DummyAlgChoice([1, 1, 1])) == false
    @test _auto_stiff_switched(DummyAlgChoice([1, 2, 1])) == true

    solver_args = args_eff_single
    prob_simple = ODEProblem(
        (du, u, p, t) -> begin
            du[1] = -u[1]
        end,
        [1.0],
        (0.0, 1.0)
    )

    withenv("SPACEAGORA_SOLVER_MODE" => "tsit5", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol, meta = _solve_with_solver_policy(prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "Tsit5"
        @test meta.fallback_used == false
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "rodas5p", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol, meta = _solve_with_solver_policy(prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "Rodas5P"
        @test meta.fallback_used == false
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "auto_stiff", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol, meta = _solve_with_solver_policy(prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "AutoTsit5(Rodas5P)"
        @test meta.initial_solver == "AutoTsit5"
    end
    split_prob_simple = SplitODEProblem(
        (du, u, p, t) -> begin
            du[1] = -u[1]
        end,
        (du, u, p, t) -> begin
            du[1] = -2u[1]
        end,
        [1.0],
        (0.0, 1.0)
    )
    withenv("SPACEAGORA_SOLVER_MODE" => "split_imex", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol, meta = _solve_with_solver_policy(split_prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "KenCarp4(IMEX)"
        @test meta.initial_solver == "KenCarp4"
    end
    withenv(
        "SPACEAGORA_SOLVER_MODE" => "split_imex",
        "SPACEAGORA_SPLIT_IMEX_SOLVER" => "kencarp47",
        "SPACEAGORA_SOLVER_MAXITERS" => nothing
    ) do
        sol, meta = _solve_with_solver_policy(split_prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "KenCarp47(IMEX)"
        @test meta.initial_solver == "KenCarp47"
    end

    withenv("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "6") do
        @test _multirate_fast_substeps() == 6
    end
    withenv("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "bad") do
        @test_throws ArgumentError _multirate_fast_substeps()
    end
    withenv("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "0") do
        @test_throws ArgumentError _multirate_fast_substeps()
    end

    withenv("SPACEAGORA_MULTIRATE_SLOW_DT_S" => nothing) do
        @test _multirate_slow_dt_s(solver_args) == min(solver_args.integration_tolerances.dt_max_orbit, 2.0)
    end
    withenv("SPACEAGORA_MULTIRATE_SLOW_DT_S" => "0.4") do
        @test _multirate_slow_dt_s(solver_args) == 0.4
    end
    withenv("SPACEAGORA_MULTIRATE_SLOW_DT_S" => "bad") do
        @test_throws ArgumentError _multirate_slow_dt_s(solver_args)
    end
    withenv("SPACEAGORA_MULTIRATE_SLOW_DT_S" => "0.0") do
        @test_throws ArgumentError _multirate_slow_dt_s(solver_args)
    end

    withenv("SPACEAGORA_MULTIRATE_SLOW_SOLVER" => "auto") do
        spec = _multirate_slow_solver_spec()
        @test spec.label == "AutoTsit5(Rodas5P)"
        @test spec.auto_switch_capable == true
    end
    withenv("SPACEAGORA_MULTIRATE_SLOW_SOLVER" => "rodas5p") do
        spec = _multirate_slow_solver_spec()
        @test spec.label == "Rodas5P"
        @test spec.auto_switch_capable == false
    end
    withenv("SPACEAGORA_MULTIRATE_FAST_SOLVER" => "kencarp4") do
        spec = _multirate_fast_solver_spec()
        @test spec.label == "KenCarp4"
        @test spec.auto_switch_capable == false
    end
    withenv("SPACEAGORA_MULTIRATE_FAST_SOLVER" => "unsupported") do
        @test_throws ArgumentError _multirate_fast_solver_spec()
    end

    multirate_subprob = _split_subproblem(split_prob_simple, split_prob_simple.f.f1, [1.0], (0.0, 0.5))
    @test multirate_subprob.tspan == (0.0, 0.5)

    @test_throws ArgumentError _solve_with_multirate_solver(prob_simple, solver_args, 1e-8, 1e-8)

    split_prob_zero = SplitODEProblem(
        (du, u, p, t) -> begin
            du[1] = -u[1]
        end,
        (du, u, p, t) -> begin
            du[1] = -2u[1]
        end,
        [1.0],
        (1.0, 1.0)
    )
    withenv("SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol_zero, meta_zero = _solve_with_multirate_solver(split_prob_zero, solver_args, 1e-8, 1e-8)
        @test string(sol_zero.retcode) == "Success"
        @test meta_zero.macro_steps == 0
        @test meta_zero.fast_substeps == 0
        @test meta_zero.slow_solver == "Tsit5"
        @test meta_zero.fast_solver == "Tsit5"
    end

    withenv(
        "SPACEAGORA_SOLVER_MAXITERS" => nothing,
        "SPACEAGORA_MULTIRATE_SLOW_SOLVER" => "tsit5",
        "SPACEAGORA_MULTIRATE_FAST_SOLVER" => "tsit5",
        "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "4",
        "SPACEAGORA_MULTIRATE_SLOW_DT_S" => "0.2"
    ) do
        sol_mr, meta_mr = _solve_with_multirate_solver(split_prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol_mr.retcode) == "Success"
        @test meta_mr.macro_steps >= 1
        @test meta_mr.fast_substeps == 4
        @test meta_mr.slow_solver == "Tsit5"
        @test meta_mr.fast_solver == "Tsit5"
        @test isapprox(meta_mr.slow_dt_s, 0.2; atol=0.0, rtol=0.0)
        @test isapprox(meta_mr.fast_dt_s, 0.05; atol=0.0, rtol=0.0)
        @test meta_mr.auto_switch_events == 0
    end

    withenv(
        "SPACEAGORA_SOLVER_MODE" => "multirate",
        "SPACEAGORA_SOLVER_MAXITERS" => nothing,
        "SPACEAGORA_MULTIRATE_SLOW_SOLVER" => "tsit5",
        "SPACEAGORA_MULTIRATE_FAST_SOLVER" => "tsit5",
        "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "4",
        "SPACEAGORA_MULTIRATE_SLOW_DT_S" => "0.2"
    ) do
        sol, meta = _solve_with_solver_policy(split_prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test occursin("Multirate(Strang;", meta.solver)
        @test meta.initial_solver == "Tsit5"
    end

    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES" => "oops") do
        @test_throws ArgumentError _ephemeris_reuse_max_entries()
    end

    reuse_cache = Dict{Any, SRPSunEphemerisCache}()
    reuse_value_a = SRPSunEphemerisCache([0.0], SVector{3, Float64}[SVector{3, Float64}(1.0, 0.0, 0.0)])
    reuse_value_b = SRPSunEphemerisCache([0.0], SVector{3, Float64}[SVector{3, Float64}(2.0, 0.0, 0.0)])
    @test _ephemeris_reuse_store!(reuse_cache, :k1, reuse_value_a, 0) === reuse_value_a
    @test !haskey(reuse_cache, :k1)
    _ephemeris_reuse_store!(reuse_cache, :k1, reuse_value_a, 2)
    @test _ephemeris_reuse_store!(reuse_cache, :k1, reuse_value_b, 2) === reuse_value_a
    _ephemeris_reuse_store!(reuse_cache, :k2, reuse_value_b, 1)
    @test !haskey(reuse_cache, :k1)
    @test haskey(reuse_cache, :k2)

    withenv("SPACEAGORA_SOLVER_MAXITERS" => "1000") do
        sol = _solve_with_explicit_solver(prob_simple, solver_args, Tsit5(), 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test !isempty(sol.t)
    end

    sc_nan_inertia = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    sc_nan_inertia.inertia_tensor = SMatrix{3, 3, Float64}(1.0, 0.0, 0.0, 0.0, NaN, 0.0, 0.0, 0.0, 1.0)
    args_nan_inertia = build_config(
        spacecraft=sc_nan_inertia,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    @test_throws ArgumentError _validate_orientation_inertia!(args_nan_inertia)

    sc_nonsym_inertia = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    sc_nonsym_inertia.inertia_tensor = SMatrix{3, 3, Float64}(2.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0)
    args_nonsym_inertia = build_config(
        spacecraft=sc_nonsym_inertia,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    @test_throws ArgumentError _validate_orientation_inertia!(args_nonsym_inertia)

    struct ThermalModelNoHeatRate <: SimulationModel.AbstractThermalModel end
    @test_throws ArgumentError _validate_thermal_model_support!((environment_model=(thermal_model=ThermalModelNoHeatRate(),),))

    withenv("SPACEAGORA_SRP_EPHEMERIS_CACHE_DT_S" => "12.5") do
        @test _srp_ephemeris_cache_dt_s() == 12.5
    end
    withenv("SPACEAGORA_SRP_EPHEMERIS_CACHE_DT_S" => "bad") do
        @test_throws ArgumentError _srp_ephemeris_cache_dt_s()
    end
    withenv("SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES" => "1234") do
        @test _srp_ephemeris_cache_max_samples() == 1234
    end
    withenv("SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES" => "bad") do
        @test_throws ArgumentError _srp_ephemeris_cache_max_samples()
    end
    withenv("SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S" => "15.0") do
        @test _nbody_ephemeris_cache_dt_s() == 15.0
    end
    withenv("SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S" => "bad") do
        @test_throws ArgumentError _nbody_ephemeris_cache_dt_s()
    end
    withenv("SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES" => "1234") do
        @test _nbody_ephemeris_cache_max_samples() == 1234
    end
    withenv("SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES" => "bad") do
        @test_throws ArgumentError _nbody_ephemeris_cache_max_samples()
    end
    withenv("SPACEAGORA_PLANET_FRAME_CACHE_DT_S" => "2.5") do
        @test _planet_frame_cache_dt_s() == 2.5
    end
    withenv("SPACEAGORA_PLANET_FRAME_CACHE_DT_S" => "bad") do
        @test_throws ArgumentError _planet_frame_cache_dt_s()
    end
    withenv("SPACEAGORA_PLANET_FRAME_CACHE_MAX_SAMPLES" => "4321") do
        @test _planet_frame_cache_max_samples() == 4321
    end
    withenv("SPACEAGORA_PLANET_FRAME_CACHE_MAX_SAMPLES" => "bad") do
        @test_throws ArgumentError _planet_frame_cache_max_samples()
    end
    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "1") do
        @test _ephemeris_reuse_enabled() == true
    end
    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "0") do
        @test _ephemeris_reuse_enabled() == false
    end
    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "bad") do
        @test_throws ArgumentError _ephemeris_reuse_enabled()
    end
    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES" => "7") do
        @test _ephemeris_reuse_max_entries() == 7
    end
    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES" => "-1") do
        @test_throws ArgumentError _ephemeris_reuse_max_entries()
    end
    withenv("SPACEAGORA_EFFECTOR_LONG_ORBIT_THRESHOLD" => "9") do
        @test _effector_long_orbit_threshold() == 9
        args_orbit_mission = (
            mission_configuration=(
                mission_type=SimulationModel.MissionOrbits,
                number_of_orbits=9,
                mission_time=1.0
            ),
        )
        @test _mission_is_long_for_effector_threads(args_orbit_mission)
    end
    @test _has_active_srp_effector((SolarRadiationPressureModel(1.2, 10.0),))

    args_eff_unsupported = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), ConstantForceModel(SVector{3, Float64}(1.0, 0.0, 0.0))),
        keplerian=true
    )
    decision_unsupported = _dynamic_effector_thread_decision(args_eff_unsupported, args_eff_unsupported.dynamics_model.dynamic_effectors, 1)
    @test decision_unsupported.use_threads == false
    @test decision_unsupported.policy_applied == false

    p_workspace_resize = ODEParams{1}(args=args_eff_single)
    resize!(p_workspace_resize.shared_buffers.harmonics_workspaces, 0)
    resize!(p_workspace_resize.shared_buffers.nbody_workspaces, 0)
    resize!(p_workspace_resize.shared_buffers.aero_workspaces, 0)
    _initialize_harmonics_workspace_buffers!(p_workspace_resize)
    _initialize_nbody_workspace_buffers!(p_workspace_resize)
    _initialize_aero_workspace_buffers!(p_workspace_resize)
    @test length(p_workspace_resize.shared_buffers.harmonics_workspaces) == 1
    @test length(p_workspace_resize.shared_buffers.nbody_workspaces) == 1
    @test length(p_workspace_resize.shared_buffers.aero_workspaces) == 1
    @test all(isnothing, p_workspace_resize.shared_buffers.harmonics_workspaces)
    @test all(isnothing, p_workspace_resize.shared_buffers.nbody_workspaces)
    @test all(isnothing, p_workspace_resize.shared_buffers.aero_workspaces)

    gram_model_instances = SimulationModel.GRAMAtmosphereModel(planet_name="earth")
    args_density_instances = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=gram_model_instances,
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    p_density_instances = ODEParams{1}(args=args_density_instances)
    withenv("SPACEAGORA_GRAM_PER_SAT_INSTANCES" => "on") do
        SimulationEngine._initialize_density_model_instances!(p_density_instances)
    end
    @test length(p_density_instances.shared_buffers.density_models) == 1
    @test p_density_instances.shared_buffers.density_models[1] isa GRAMAtmosphereModel
    @test p_density_instances.shared_buffers.density_models[1] !== gram_model_instances
    SimulationEngine._initialize_density_cache_buffers!(p_density_instances)
    @test length(p_density_instances.shared_buffers.gram_density_cache) == 1
    @test p_density_instances.shared_buffers.gram_density_cache[1] === nothing
    cache_typed = @inferred SimulationModel.SimulationCallbacks._gram_density_cache_for_sat!(
        p_density_instances.shared_buffers.gram_density_cache,
        1
    )
    @test cache_typed isa GramTrackCache
    typed_density_models = GRAMAtmosphereModel[gram_model_instances]
    selected_density_model = @inferred SimulationModel.SimulationCallbacks._density_model_for_sat(
        typed_density_models,
        gram_model_instances,
        1
    )
    @test selected_density_model === gram_model_instances
    pool_models_typed, pool_locks_typed = @inferred SimulationModel.SimulationCallbacks._ensure_gram_isolated_pool!(
        p_density_instances,
        gram_model_instances,
        1
    )
    @test pool_models_typed === p_density_instances.shared_buffers.gram_isolated_pool_models
    @test pool_locks_typed === p_density_instances.shared_buffers.gram_isolated_pool_locks
    @test length(pool_models_typed) == 1
    @test pool_models_typed[1] isa GRAMAtmosphereModel

    gram_surrogate_instances = SimulationModel.GRAMAtmosphereModelSurrogate(gram_model_instances, "unused", nothing)
    args_density_surrogate = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=gram_surrogate_instances,
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    p_density_surrogate = ODEParams{1}(args=args_density_surrogate)
    withenv("SPACEAGORA_GRAM_PER_SAT_INSTANCES" => "on") do
        SimulationEngine._initialize_density_model_instances!(p_density_surrogate)
    end
    @test length(p_density_surrogate.shared_buffers.density_models) == 1
    @test p_density_surrogate.shared_buffers.density_models[1] isa GRAMAtmosphereModelSurrogate
    @test p_density_surrogate.shared_buffers.density_models[1] !== gram_surrogate_instances

    args_srp = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), SolarRadiationPressureModel(1.2, 12.0)),
        keplerian=true
    )
    p_srp = ODEParams{1}(args=args_srp)
    _initialize_srp_sun_cache_buffer!(p_srp)
    withenv("SPACEAGORA_SRP_EPHEMERIS_CACHE" => "1") do
        _initialize_srp_sun_ephemeris_cache!(p_srp, 0.0, 0.0)
    end
    @test p_srp.shared_buffers.srp_sun_ephemeris_cache[] === nothing
    withenv(
        "SPACEAGORA_SRP_EPHEMERIS_CACHE" => "1",
        "SPACEAGORA_SRP_EPHEMERIS_CACHE_DT_S" => "1.0",
        "SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES" => "2"
    ) do
        @test_logs (:warn, r"SRP ephemeris cache disabled") _initialize_srp_sun_ephemeris_cache!(p_srp, 0.0, 10.0)
    end
    @test p_srp.shared_buffers.srp_sun_ephemeris_cache[] === nothing

    _clear_ephemeris_reuse_cache!()
    p_srp_reuse_a = ODEParams{1}(args=args_srp)
    p_srp_reuse_b = ODEParams{1}(args=args_srp)
    _initialize_srp_sun_cache_buffer!(p_srp_reuse_a)
    _initialize_srp_sun_cache_buffer!(p_srp_reuse_b)
    _reset_spice_runtime_counters!(p_srp_reuse_a)
    _reset_spice_runtime_counters!(p_srp_reuse_b)
    withenv(
        "SPACEAGORA_SRP_EPHEMERIS_CACHE" => "1",
        "SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "1",
        "SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES" => "16",
        "SPACEAGORA_SRP_EPHEMERIS_CACHE_DT_S" => "10.0",
        "SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES" => "100"
    ) do
        _initialize_srp_sun_ephemeris_cache!(p_srp_reuse_a, 0.0, 10.0)
        cache_a = p_srp_reuse_a.shared_buffers.srp_sun_ephemeris_cache[]
        @test cache_a isa SRPSunEphemerisCache
        @test p_srp_reuse_a.shared_buffers.spice_runtime_counters.srp_spkpos_cache_build_calls[] == 2
        _initialize_srp_sun_ephemeris_cache!(p_srp_reuse_b, 0.0, 10.0)
        cache_b = p_srp_reuse_b.shared_buffers.srp_sun_ephemeris_cache[]
        @test cache_b === cache_a
        @test p_srp_reuse_b.shared_buffers.spice_runtime_counters.srp_spkpos_cache_build_calls[] == 0
    end
    _clear_ephemeris_reuse_cache!()

    args_nbody = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), NBodyGravityModel(body_names=("moon",), primary_body_name="Earth")),
        keplerian=true
    )
    p_nbody = ODEParams{1}(args=args_nbody)
    _initialize_nbody_ephemeris_cache_buffer!(p_nbody)
    withenv("SPACEAGORA_NBODY_EPHEMERIS_CACHE" => "1") do
        _initialize_nbody_ephemeris_cache!(p_nbody, 0.0, 0.0)
    end
    @test p_nbody.shared_buffers.nbody_ephemeris_cache[] === nothing
    withenv(
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE" => "1",
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S" => "1.0",
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES" => "2"
    ) do
        @test_logs (:warn, r"N-body ephemeris cache disabled") _initialize_nbody_ephemeris_cache!(p_nbody, 0.0, 10.0)
    end
    @test p_nbody.shared_buffers.nbody_ephemeris_cache[] === nothing

    p_planet_frame = ODEParams{1}(args=args_srp)
    _initialize_planet_frame_cache_buffer!(p_planet_frame)
    withenv("SPACEAGORA_PLANET_FRAME_CACHE" => "1") do
        _initialize_planet_frame_ephemeris_cache!(p_planet_frame, 0.0, 0.0)
    end
    @test p_planet_frame.shared_buffers.planet_frame_ephemeris_cache[] === nothing
    withenv(
        "SPACEAGORA_PLANET_FRAME_CACHE" => "1",
        "SPACEAGORA_PLANET_FRAME_CACHE_DT_S" => "1.0",
        "SPACEAGORA_PLANET_FRAME_CACHE_MAX_SAMPLES" => "2"
    ) do
        @test_logs (:warn, r"Planet frame cache disabled") _initialize_planet_frame_ephemeris_cache!(p_planet_frame, 0.0, 10.0)
    end
    @test p_planet_frame.shared_buffers.planet_frame_ephemeris_cache[] === nothing

    args_nbody_srp = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            NBodyGravityModel(body_names=("moon", "sun"), primary_body_name="Earth"),
            SolarRadiationPressureModel(1.2, 12.0)
        ),
        keplerian=true
    )
    p_nbody_srp = ODEParams{1}(args=args_nbody_srp)
    _initialize_nbody_ephemeris_cache_buffer!(p_nbody_srp)
    _initialize_srp_sun_cache_buffer!(p_nbody_srp)
    _reset_spice_runtime_counters!(p_nbody_srp)
    _reset_spice_rhs_memo!(p_nbody_srp)
    p_nbody_srp.shared_buffers.et_start[] = 0.0
    p_nbody_srp.shared_buffers.current_time[] = 123.0
    nbody_model = args_nbody_srp.dynamics_model.dynamic_effectors[2]
    srp_model = args_nbody_srp.dynamics_model.dynamic_effectors[3]
    x_a = Float64[7000e3, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0]
    x_b = Float64[6999e3, 10.0, 0.0, 0.0, 0.0, 0.0, 200.0]
    SimulationModel.calcForceTorque(nbody_model, x_a, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls[] == 2
    SimulationModel.calcForceTorque(nbody_model, x_b, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls[] == 2
    SimulationModel.calcForceTorque(srp_model, x_a, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls[] == 0
    p_nbody_srp.shared_buffers.current_time[] = 124.0
    SimulationModel.calcForceTorque(srp_model, x_a, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls[] == 1
    p_nbody_srp.shared_buffers.spice_rhs_memo_enabled[] = false
    _reset_spice_runtime_counters!(p_nbody_srp)
    _reset_spice_rhs_memo!(p_nbody_srp)
    p_nbody_srp.shared_buffers.current_time[] = 223.0
    SimulationModel.calcForceTorque(nbody_model, x_a, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls[] == 2
    SimulationModel.calcForceTorque(nbody_model, x_b, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls[] == 4
    SimulationModel.calcForceTorque(srp_model, x_a, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls[] == 1

    dyn = SimulationModel.DynamicEffectors
    harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
    @test_throws ArgumentError GravitationalHarmonicsModel(-1, 0, harmonics_file, EARTH)
    @test_throws ArgumentError GravitationalHarmonicsModel(10_000, 0, harmonics_file, EARTH)

    nbody_ws = dyn._make_nbody_scratch_workspace(1)
    dyn._ensure_nbody_workspace_capacity!(nbody_ws, 3, 4)
    @test length(nbody_ws.pos_primary_k_all) == 3
    @test length(nbody_ws.thread_force) == 4
    nbody_ws_typed = @inferred dyn._nbody_workspace_for_sat!(p_nbody_srp, 1, 2, 2)
    @test nbody_ws_typed isa NBodyScratchWorkspace
    @test p_nbody_srp.shared_buffers.nbody_workspaces[1] === nbody_ws_typed
    nbody_ws_oob = @inferred dyn._nbody_workspace_for_sat!(p_nbody_srp, 5, 2, 2)
    @test length(nbody_ws_oob.pos_primary_k_all) == 2
    @test length(nbody_ws_oob.thread_force) == 2

    aero_ws_typed = @inferred dyn._aero_workspace_for_sat!(p_workspace_resize, 1, 2)
    @test aero_ws_typed isa AeroScratchWorkspace
    @test p_workspace_resize.shared_buffers.aero_workspaces[1] === aero_ws_typed

    harmonics_model = GravitationalHarmonicsModel(8, 8, harmonics_file, EARTH)
    harmonics_ws_typed = @inferred dyn._harmonics_workspace_for_sat!(harmonics_model, p_workspace_resize, 1)
    @test harmonics_ws_typed isa HarmonicsScratchWorkspace
    @test p_workspace_resize.shared_buffers.harmonics_workspaces[1] isa Dict{UInt, HarmonicsScratchWorkspace}

    nbody_positions = reshape(
        SVector{3, Float64}[
            SVector{3, Float64}(1.0, 0.0, 0.0),
            SVector{3, Float64}(2.0, 0.0, 0.0),
            SVector{3, Float64}(3.0, 0.0, 0.0),
            SVector{3, Float64}(4.0, 0.0, 0.0)
        ],
        4,
        1
    )
    nbody_cache = NBodyEphemerisCache(
        "earth_barycenter",
        ["moon_barycenter"],
        Dict("moon_barycenter" => 1),
        [0.0, 5.0, 10.0, 15.0],
        nbody_positions
    )
    @test dyn._nbody_body_position_from_cache_j2000(nbody_cache, -1.0, "moon_barycenter", "earth_barycenter") === nothing
    @test dyn._nbody_body_position_from_cache_j2000(nbody_cache, NaN, "moon_barycenter", "earth_barycenter") isa SVector{3, Float64}
    @test dyn._nbody_body_position_from_cache_j2000(nbody_cache, 15.0, "moon_barycenter", "earth_barycenter") == nbody_positions[4, 1]
    @test dyn._nbody_body_position_from_cache_j2000(nbody_cache, 2.5, "moon_barycenter", "earth_barycenter") == SVector{3, Float64}(1.5, 0.0, 0.0)
    nbody_interp_catmull = dyn._nbody_body_position_from_cache_j2000(nbody_cache, 7.5, "moon_barycenter", "earth_barycenter")
    @test nbody_interp_catmull isa SVector{3, Float64}

    srp_cache = SRPSunEphemerisCache(
        [0.0, 5.0, 10.0, 15.0],
        SVector{3, Float64}[
            SVector{3, Float64}(1.0, 0.0, 0.0),
            SVector{3, Float64}(2.0, 0.0, 0.0),
            SVector{3, Float64}(3.0, 0.0, 0.0),
            SVector{3, Float64}(4.0, 0.0, 0.0)
        ]
    )
    @test dyn._srp_sun_position_from_cache_j2000(srp_cache, -1.0) === nothing
    @test dyn._srp_sun_position_from_cache_j2000(srp_cache, NaN) isa SVector{3, Float64}
    @test dyn._srp_sun_position_from_cache_j2000(srp_cache, 15.0) == srp_cache.positions_j2000[end]
    @test dyn._srp_sun_position_from_cache_j2000(srp_cache, 2.5) == SVector{3, Float64}(1.5, 0.0, 0.0)
    srp_interp_catmull = dyn._srp_sun_position_from_cache_j2000(srp_cache, 7.5)
    @test srp_interp_catmull isa SVector{3, Float64}

    force_zero_srp, torque_zero_srp = calcForceTorque(SolarRadiationPressureModel(1.2, 0.0), x_a, p_nbody_srp, 1)
    @test force_zero_srp == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test torque_zero_srp == SVector{3, Float64}(0.0, 0.0, 0.0)

    p_srp.shared_buffers.srp_sun_ephemeris_cache[] = SRPSunEphemerisCache(
        [0.0, 10.0],
        SVector{3, Float64}[SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)]
    )
    p_srp.shared_buffers.et_start[] = 0.0
    p_srp.shared_buffers.current_time[] = 5.0
    x_zero_dist = Float64[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0]
    force_zero_dist, torque_zero_dist = calcForceTorque(SolarRadiationPressureModel(1.2, 12.0), x_zero_dist, p_srp, 1)
    @test force_zero_dist == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test torque_zero_dist == SVector{3, Float64}(0.0, 0.0, 0.0)

    @test dyn.eclipse_area_calc(SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(1.0, 0.0, 0.0), EARTH.Rp_e) == 1.0
    @test dyn.eclipse_area_calc(
        SVector{3, Float64}(-5.271937279754128e6, -4.185218527153555e6, -1.0434271238606143e6),
        SVector{3, Float64}(-1.1107937751409042e10, 1.6143690900498734e11, 8.75872644289004e10),
        EARTH.Rp_e
    ) == 0.0
    annular_ratio = dyn.eclipse_area_calc(
        SVector{3, Float64}(0.0, 0.0, 1.0e10),
        SVector{3, Float64}(0.0, 0.0, -1.5e11),
        EARTH.Rp_e
    )
    @test 0.0 < annular_ratio < 1.0
    partial_ratio = dyn.eclipse_area_calc(
        SVector{3, Float64}(6.930027129188876e6, -2.6352977471555886e6, -3.363004422597388e6),
        SVector{3, Float64}(-9.438128429326639e10, -6.696979722657439e10, 1.7822072441008075e11),
        EARTH.Rp_e
    )
    @test 0.0 < partial_ratio < 1.0
    none_ratio = dyn.eclipse_area_calc(
        SVector{3, Float64}(1.2249535847697716e7, -5.782145543435082e6, -7.299266925237677e6),
        SVector{3, Float64}(-4.937007687846062e10, -7.233778731136734e10, 1.7733640551002377e11),
        EARTH.Rp_e
    )
    @test isapprox(none_ratio, 1.0; atol=1e-12, rtol=0.0)

    @testset "SRP Regression Contracts" begin
        # On/off delta: SRP-enabled model must produce non-zero force while area=0 is zero force.
        p_srp.shared_buffers.srp_sun_ephemeris_cache[] = SRPSunEphemerisCache(
            [0.0, 10.0],
            SVector{3, Float64}[
                SVector{3, Float64}(149_597_870.7, 0.0, 0.0),
                SVector{3, Float64}(149_597_870.7, 0.0, 0.0)
            ]
        )
        p_srp.shared_buffers.et_start[] = 0.0
        p_srp.shared_buffers.current_time[] = 5.0
        x_srp_regression = Float64[7_000_000.0, 100.0, 50.0, 0.0, 0.0, 0.0, 200.0]

        force_on, _ = calcForceTorque(SolarRadiationPressureModel(1.2, 12.0), x_srp_regression, p_srp, 1)
        force_off, _ = calcForceTorque(SolarRadiationPressureModel(1.2, 0.0), x_srp_regression, p_srp, 1)
        @test norm(force_on) > 0.0
        @test force_off == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test norm(force_on - force_off) > 0.0

        # Eclipse attenuation: partial eclipse force norm must scale by eclipse ratio.
        pos_partial = SVector{3, Float64}(6.930027129188876e6, -2.6352977471555886e6, -3.363004422597388e6)
        sun_partial = SVector{3, Float64}(-9.438128429326639e10, -6.696979722657439e10, 1.7822072441008075e11)
        partial_ratio_regression = dyn.eclipse_area_calc(pos_partial, sun_partial, EARTH.Rp_e)
        force_partial = dyn.srp_cannonball_accel(pos_partial, sun_partial, EARTH.Rp_e, 4.56e-6, 1.2, 12.0, 200.0)
        force_no_eclipse = dyn.srp_cannonball_accel(pos_partial, sun_partial, 0.0, 4.56e-6, 1.2, 12.0, 200.0)
        @test 0.0 < partial_ratio_regression < 1.0
        @test norm(force_no_eclipse) > 0.0
        @test isapprox(norm(force_partial), partial_ratio_regression * norm(force_no_eclipse); atol=0.0, rtol=1e-12)

        # Solver parity: legacy solver helper (`srp`) must match canonical SRP kernel for same geometry.
        planet = args_srp.environment_model.planet
        pos_solver = SVector{3, Float64}(x_a[1:3])
        et_solver = 123.0
        primary_body_name = dyn._spice_query_name(planet.name)
        sun_j2000 = lock(SimulationModel.SPICE_LOCK) do
            SVector{3, Float64}(spkpos("sun", et_solver, "J2000", "none", primary_body_name)[1])
        end
        sun_pci = SVector{3, Float64}(planet.J2000_to_pci * sun_j2000 * 1e3)
        force_kernel = dyn.srp_cannonball_accel(pos_solver, sun_pci, planet.Rp_e, 4.56e-6, 1.2, 12.0, 200.0)
        force_solver = dyn.srp(planet, 4.56e-6, 1.2, 12.0, 200.0, pos_solver, et_solver)
        @test isapprox(force_solver, force_kernel; atol=1e-12, rtol=1e-10)
    end

    default_ckpt_settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    args_default_ckpt = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=default_ckpt_settings
    )
    @test _checkpoint_directory(args_default_ckpt) == joinpath("output", "checkpoints")

    mktempdir() do tmp
        bad_ckpt_settings = SimulationSettings(
            results=true,
            verbose=false,
            generate_plots=false,
            results_directory="output",
            save_csv=true,
            normalize=false,
            checkpoint_directory=tmp
        )
        args_bad_ckpt = build_config(
            spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=10.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true,
            simulation_settings=bad_ckpt_settings
        )
        paths_bad_ckpt = _checkpoint_paths(args_bad_ckpt)
        mkpath(dirname(paths_bad_ckpt.data))
        open(paths_bad_ckpt.data, "w") do io
            serialize(io, Dict(:invalid => true))
        end
        @test_throws ArgumentError _load_checkpoint(args_bad_ckpt)
    end

    @test _find_sample_value([nothing, nothing]) === nothing
    results_df = DataFrame()
    _append_series_columns!(results_df, "meta", [(a=1, b=2), (a=3, b=4)])
    _append_series_columns!(results_df, "dictmeta", [Dict(:z => 1, :a => 2), Dict(:z => 3, :a => 4)])
    struct BoxedValue
        value::Int
    end
    _append_series_columns!(results_df, "boxed", [BoxedValue(1), BoxedValue(2)])
    @test results_df.meta_a == [1, 3]
    @test results_df.meta_b == [2, 4]
    @test results_df.dictmeta_a == [2, 4]
    @test results_df.dictmeta_z == [1, 3]
    @test length(results_df.boxed) == 2

    mktempdir() do tmp
        out_path = joinpath(tmp, "artifact.txt")
        @test_throws ErrorException _atomic_write_file(out_path, tmp_path -> begin
            write(tmp_path, "tmp-data")
            throw(ErrorException("forced writer failure"))
        end)
        @test !isfile(out_path)
        @test isempty(readdir(tmp))
    end

    args_heat_copy = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    u_heat_copy = build_initial_conditions(args_heat_copy)
    du_heat_copy = copy(u_heat_copy)
    du_heat_copy .= 0.0
    p_heat_copy = ODEParams{1}(args=args_heat_copy)
    _initialize_heat_rate_buffers!(p_heat_copy)
    p_heat_copy.shared_buffers.heat_rates[1] = Float64[1.25, 9.0]
    spacecraft_dynamics!(du_heat_copy, u_heat_copy, p_heat_copy, 0.0)
    @test du_heat_copy.sc[1].heat_loads[1] == 1.25

    @test_throws ArgumentError _resolve_component_tolerance(-1.0, 1.0, "unit_test_tol")
end

@testset "Rigid-Body Angular Dynamics Uses Inertia Tensor" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.02, -0.03, 0.015)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    inertia_tensor = SMatrix{3, 3, Float64}(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0)
    sc.inertia_tensor = inertia_tensor
    applied_torque = SVector{3, Float64}(0.12, -0.08, 0.05)

    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(ConstantTorqueModel(applied_torque),),
        keplerian=true
    )

    u0 = build_initial_conditions(args)
    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams{1}(args=args)
    spacecraft_dynamics!(du0, u0, p, 0.0)

    ω = SVector{3, Float64}(u0.sc[1].ω)
    qdot_expected = SimulationModel.DynamicsRotational.quaternion_derivative(ω, SVector{4, Float64}(u0.sc[1].q))
    ωdot_expected = SimulationModel.DynamicsRotational.angular_acceleration(
        ω,
        inertia_tensor,
        applied_torque;
        include_gyroscopic=true
    )
    ωdot_legacy = inertia_tensor \ (applied_torque - cross(ω, inertia_tensor * ω))
    @test isapprox(SVector{4, Float64}(du0.sc[1].q), qdot_expected; atol=1e-12, rtol=1e-10)
    @test isapprox(ωdot_expected, ωdot_legacy; atol=1e-12, rtol=1e-10)
    @test isapprox(SVector{3, Float64}(du0.sc[1].ω), ωdot_expected; atol=1e-12, rtol=1e-10)
end

@testset "Orientation Simulation Rejects Invalid Inertia Tensor" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    sc.inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))

    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )

    @test_throws ArgumentError run_simulation(args)
end

@testset "Heat Loads Are Not Coupled To Force Magnitude" begin
    sc = make_single_link_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(ConstantForceModel(SVector{3, Float64}(1.0e5, -2.0e5, 3.0e5)),),
        keplerian=true
    )
    args.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))

    u0 = build_initial_conditions(args)
    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams{1}(args=args)
    spacecraft_dynamics!(du0, u0, p, 0.0)

    @test norm(SVector{3, Float64}(du0.sc[1].vel)) > 0.0
    @test all(==(0.0), du0.sc[1].heat_loads)
end

@testset "Heat Load Derivative Uses Thermal Model" begin
    sc = make_single_link_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    args.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))

    u0 = build_initial_conditions(args)
    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams{1}(args=args)
    _initialize_heat_rate_buffers!(p)
    p.shared_buffers.densities[1] = 1.0e-6
    p.shared_buffers.temperatures[1] = 250.0
    p.shared_buffers.winds[1] = SVector{3, Float64}(0.0, 0.0, 0.0)
    thermal_cb = SimulationModel.SimulationCallbacks.get_thermal_callback(1, args)
    thermal_cb.affect!((p=p, u=u0, t=0.0))
    spacecraft_dynamics!(du0, u0, p, 0.0)

    @test all(isfinite, p.shared_buffers.heat_rates[1])
    @test any(>(0.0), p.shared_buffers.heat_rates[1])
    @test all(isfinite, du0.sc[1].heat_loads)
    @test any(>(0.0), du0.sc[1].heat_loads)
end

@testset "Drag Dissipates Specific Orbital Energy" begin
    sc = make_spacecraft(
        ra_alt_m=220e3,
        rp_alt_m=100e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=180.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )

    df = run_case(args)
    eps = specific_energy(df, EARTH.μ)
    @test last(eps) < first(eps) - 1e5
end

@testset "AGORA Earth Regression (Golden)" begin
    golden_path = joinpath(REPO_ROOT, "test", "golden", "agora_earth_regression.csv")
    @test isfile(golden_path)
    golden = CSV.read(golden_path, DataFrame)

    sc = make_agora_earth_spacecraft()
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=6.0 * 3600.0,
        EI_km=300.0,
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        keplerian=true,
        initial_time=SimulationModel.InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=15.0
        )
    )

    df = run_case(args)
    @test nrow(df) > 1000
    times = Vector{Float64}(df.time)

    for row in eachrow(golden)
        t = Float64(row.time)
        pos_atol = Float64(row.pos_atol_m)
        vel_atol = Float64(row.vel_atol_mps)

        pos1 = interp_linear(times, df.sc1_pos_1, t)
        pos2 = interp_linear(times, df.sc1_pos_2, t)
        pos3 = interp_linear(times, df.sc1_pos_3, t)
        vel1 = interp_linear(times, df.sc1_vel_1, t)
        vel2 = interp_linear(times, df.sc1_vel_2, t)
        vel3 = interp_linear(times, df.sc1_vel_3, t)

        @test isapprox(pos1, Float64(row.pos1); atol=pos_atol, rtol=0.0)
        @test isapprox(pos2, Float64(row.pos2); atol=pos_atol, rtol=0.0)
        @test isapprox(pos3, Float64(row.pos3); atol=pos_atol, rtol=0.0)
        @test isapprox(vel1, Float64(row.vel1); atol=vel_atol, rtol=0.0)
        @test isapprox(vel2, Float64(row.vel2); atol=vel_atol, rtol=0.0)
        @test isapprox(vel3, Float64(row.vel3); atol=vel_atol, rtol=0.0)
    end
end

@testset "N-Body Gravity Typed API Smoke" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), NBodyGravityModel(["Sun"], "Earth", SPICE_PATH)),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df = run_case(args)
    @test nrow(df) > 10
    @test all(isfinite, df.sc1_pos_1)
    @test all(isfinite, df.sc1_vel_1)
end

@testset "Two-Spacecraft Isolation vs Single-Craft Baselines" begin
    sc_a = make_spacecraft(ra_alt_m=520e3, rp_alt_m=430e3, ν_deg=150.0)
    sc_b = make_spacecraft(ra_alt_m=700e3, rp_alt_m=650e3, ν_deg=110.0)

    args_multi = build_config_multi(
        spacecraft=[sc_a, sc_b],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    args_a = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=430e3, ν_deg=150.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    args_b = build_config(
        spacecraft=make_spacecraft(ra_alt_m=700e3, rp_alt_m=650e3, ν_deg=110.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df_multi = run_case_silent(args_multi)
    df_a = run_case_silent(args_a)
    df_b = run_case_silent(args_b)

    @test nrow(df_multi) > 10
    sample_idxs = round.(Int, range(1, nrow(df_multi), length=8))
    for idx in sample_idxs
        t = Float64(df_multi.time[idx])

        pa_m = SVector{3, Float64}(Float64(df_multi.sc1_pos_1[idx]), Float64(df_multi.sc1_pos_2[idx]), Float64(df_multi.sc1_pos_3[idx]))
        va_m = SVector{3, Float64}(Float64(df_multi.sc1_vel_1[idx]), Float64(df_multi.sc1_vel_2[idx]), Float64(df_multi.sc1_vel_3[idx]))
        pa_s = SVector{3, Float64}(
            interp_linear(df_a.time, df_a.sc1_pos_1, t),
            interp_linear(df_a.time, df_a.sc1_pos_2, t),
            interp_linear(df_a.time, df_a.sc1_pos_3, t)
        )
        va_s = SVector{3, Float64}(
            interp_linear(df_a.time, df_a.sc1_vel_1, t),
            interp_linear(df_a.time, df_a.sc1_vel_2, t),
            interp_linear(df_a.time, df_a.sc1_vel_3, t)
        )
        # Multi-satellite adaptive stepping can introduce modest trajectory differences vs single-body runs.
        @test norm(pa_m - pa_s) < 200.0
        @test norm(va_m - va_s) < 0.2

        pb_m = SVector{3, Float64}(Float64(df_multi.sc2_pos_1[idx]), Float64(df_multi.sc2_pos_2[idx]), Float64(df_multi.sc2_pos_3[idx]))
        vb_m = SVector{3, Float64}(Float64(df_multi.sc2_vel_1[idx]), Float64(df_multi.sc2_vel_2[idx]), Float64(df_multi.sc2_vel_3[idx]))
        pb_s = SVector{3, Float64}(
            interp_linear(df_b.time, df_b.sc1_pos_1, t),
            interp_linear(df_b.time, df_b.sc1_pos_2, t),
            interp_linear(df_b.time, df_b.sc1_pos_3, t)
        )
        vb_s = SVector{3, Float64}(
            interp_linear(df_b.time, df_b.sc1_vel_1, t),
            interp_linear(df_b.time, df_b.sc1_vel_2, t),
            interp_linear(df_b.time, df_b.sc1_vel_3, t)
        )
        @test norm(pb_m - pb_s) < 200.0
        @test norm(vb_m - vb_s) < 0.2
    end
end

@testset "Single-Link Drag Dissipates Specific Orbital Energy" begin
    sc = make_single_link_spacecraft(
        ra_alt_m=220e3,
        rp_alt_m=100e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=180.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )

    df = run_case(args)
    eps = specific_energy(df, EARTH.μ)
    @test last(eps) < first(eps) - 1e5
end
