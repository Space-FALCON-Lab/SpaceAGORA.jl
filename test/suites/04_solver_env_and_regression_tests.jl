@testset "Solver/Env Helper Parsing Coverage" begin
    withenv("SPACEAGORA_SOLVER_MODE" => nothing) do
        @test _solver_policy_mode() == :tsit5
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "auto") do
        @test _solver_policy_mode() == :auto_stiff
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "symplectic") do
        @test _solver_policy_mode() == :symplectic
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split") do
        @test _solver_policy_mode() == :gravity_backbone_split
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
    withenv("SPACEAGORA_RHS_EXECUTION_MODE" => "flat") do
        @test _rhs_execution_mode_env() == :flat_constellation_effector_queue
    end
    withenv("SPACEAGORA_RHS_EXECUTION_MODE" => "satellite") do
        @test _rhs_execution_mode_env() == :satellite_batch
    end
    withenv("SPACEAGORA_RHS_EXECUTION_MODE" => "per_satellite") do
        @test _rhs_execution_mode_env() == :per_satellite_effector_reduce
    end
    withenv("SPACEAGORA_RHS_EXECUTION_MODE" => "bad") do
        @test_throws ArgumentError _rhs_execution_mode_env()
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
    p_eff_single = ODEParams(n_sats=1, args=args_eff_single)
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
    p_eff_multi = ODEParams(n_sats=1, args=args_eff_multi)
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
    p_eff_constellation = ODEParams(n_sats=1, args=args_eff_constellation)
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

    args_flat_rhs = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3 + 10e3 * i, rp_alt_m=480e3 + 5e3 * i, ν_deg=150.0 + i)
            for i in 1:4
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            InverseSquaredJ2GravityModel(),
            InverseSquaredGravityModel(),
            InverseSquaredJ2GravityModel(),
        ),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    u_flat_rhs = build_initial_conditions(args_flat_rhs)
    du_serial_rhs = copy(u_flat_rhs)
    du_flat_rhs = copy(u_flat_rhs)
    du_serial_rhs .= 0.0
    du_flat_rhs .= 0.0
    p_serial_rhs = ODEParams(n_sats=4, args=args_flat_rhs)
    p_flat_rhs = ODEParams(n_sats=4, args=args_flat_rhs)
    _initialize_heat_rate_buffers!(p_serial_rhs)
    _initialize_heat_rate_buffers!(p_flat_rhs)
    withenv(
        "SPACEAGORA_RHS_EXECUTION_MODE" => "serial",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off"
    ) do
        spacecraft_dynamics!(du_serial_rhs, u_flat_rhs, p_serial_rhs, 0.0)
    end
    withenv(
        "SPACEAGORA_RHS_EXECUTION_MODE" => "flat",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(max(1, Threads.nthreads()))
    ) do
        plan_flat_rhs = _rhs_execution_plan(
            args_flat_rhs,
            p_flat_rhs,
            args_flat_rhs.dynamics_model.dynamic_effectors,
            4
        )
        if Threads.nthreads() > 1
            @test plan_flat_rhs.mode == :flat_constellation_effector_queue
        end
        spacecraft_dynamics!(du_flat_rhs, u_flat_rhs, p_flat_rhs, 0.0)
    end
    for i in 1:4
        @test isapprox(du_flat_rhs.sc[i].pos, du_serial_rhs.sc[i].pos; atol=1e-12, rtol=1e-12)
        @test isapprox(du_flat_rhs.sc[i].vel, du_serial_rhs.sc[i].vel; atol=1e-10, rtol=1e-10)
        @test isapprox(du_flat_rhs.sc[i].mass, du_serial_rhs.sc[i].mass; atol=1e-12, rtol=1e-12)
        @test isapprox(du_flat_rhs.sc[i].heat_loads, du_serial_rhs.sc[i].heat_loads; atol=1e-12, rtol=1e-12)
    end

    harmonics_file_flat = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
    harmonics_flat_model = GravitationalHarmonicsModel(8, 8, harmonics_file_flat, EARTH)
    args_harmonics_flat = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=520e3 + 8e3 * i, rp_alt_m=500e3 + 4e3 * i, ν_deg=120.0 + 2.0 * i)
            for i in 1:4
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(harmonics_flat_model,),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    u_harmonics_flat = build_initial_conditions(args_harmonics_flat)
    du_harmonics_serial = copy(u_harmonics_flat)
    du_harmonics_flat = copy(u_harmonics_flat)
    du_harmonics_serial .= 0.0
    du_harmonics_flat .= 0.0
    p_harmonics_serial = ODEParams(n_sats=4, args=args_harmonics_flat)
    p_harmonics_flat = ODEParams(n_sats=4, args=args_harmonics_flat)
    _initialize_heat_rate_buffers!(p_harmonics_serial)
    _initialize_heat_rate_buffers!(p_harmonics_flat)
    _initialize_harmonics_workspace_buffers!(p_harmonics_serial)
    _initialize_harmonics_workspace_buffers!(p_harmonics_flat)
    withenv(
        "SPACEAGORA_RHS_EXECUTION_MODE" => "serial",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off"
    ) do
        spacecraft_dynamics!(du_harmonics_serial, u_harmonics_flat, p_harmonics_serial, 0.0)
    end
    withenv(
        "SPACEAGORA_RHS_EXECUTION_MODE" => "flat",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off",
        "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "1",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(max(1, Threads.nthreads()))
    ) do
        plan_harmonics_flat = _rhs_execution_plan(
            args_harmonics_flat,
            p_harmonics_flat,
            args_harmonics_flat.dynamics_model.dynamic_effectors,
            4
        )
        if Threads.nthreads() > 1
            @test plan_harmonics_flat.mode == :flat_constellation_effector_queue
        end
        spacecraft_dynamics!(du_harmonics_flat, u_harmonics_flat, p_harmonics_flat, 0.0)
    end
    for i in 1:4
        @test isapprox(du_harmonics_flat.sc[i].pos, du_harmonics_serial.sc[i].pos; atol=1e-12, rtol=1e-12)
        @test isapprox(du_harmonics_flat.sc[i].vel, du_harmonics_serial.sc[i].vel; atol=1e-8, rtol=1e-9)
        @test isapprox(du_harmonics_flat.sc[i].mass, du_harmonics_serial.sc[i].mass; atol=1e-12, rtol=1e-12)
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

    withenv("SPACEAGORA_SYMPLECTIC_DT_S" => nothing) do
        @test _symplectic_fixed_dt_s(args_eff_single) == args_eff_single.integration_tolerances.dt_max_orbit
    end
    withenv("SPACEAGORA_SYMPLECTIC_DT_S" => "4.0") do
        @test _symplectic_fixed_dt_s(args_eff_single) == 4.0
    end
    withenv("SPACEAGORA_SYMPLECTIC_DT_S" => "bad") do
        @test_throws ArgumentError _symplectic_fixed_dt_s(args_eff_single)
    end
    withenv("SPACEAGORA_SYMPLECTIC_DT_S" => "0.0") do
        @test_throws ArgumentError _symplectic_fixed_dt_s(args_eff_single)
    end

    withenv("SPACEAGORA_GRAVITY_BACKBONE_DT_S" => nothing) do
        @test _gravity_backbone_fixed_dt_s(args_eff_single) == args_eff_single.integration_tolerances.dt_max_orbit
    end
    withenv("SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "5.0") do
        @test _gravity_backbone_fixed_dt_s(args_eff_single) == 5.0
    end
    withenv("SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "bad") do
        @test_throws ArgumentError _gravity_backbone_fixed_dt_s(args_eff_single)
    end
    withenv("SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "0.0") do
        @test_throws ArgumentError _gravity_backbone_fixed_dt_s(args_eff_single)
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
    solver_args_atmo = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=true
    )

    prob_simple = ODEProblem(
        (du, u, p, t) -> begin
            du[1] = -u[1]
        end,
        [1.0],
        (0.0, 1.0)
    )

    withenv("SPACEAGORA_SOLVER_MODE" => "tsit5", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol, meta = _solve_with_solver_policy(prob_simple, _active_solver_config(), solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "Tsit5"
        @test meta.fallback_used == false
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "rodas5p", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol, meta = _solve_with_solver_policy(prob_simple, _active_solver_config(), solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "Rodas5P"
        @test meta.fallback_used == false
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "auto_stiff", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol, meta = _solve_with_solver_policy(prob_simple, _active_solver_config(), solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "Tsit5"
        @test meta.initial_solver == "Tsit5"
        @test meta.fallback_used == false
    end
    withenv(
        "SPACEAGORA_SOLVER_MODE" => "auto_stiff",
        "SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5" => "0",
        "SPACEAGORA_SOLVER_MAXITERS" => nothing
    ) do
        sol, meta = _solve_with_solver_policy(prob_simple, _active_solver_config(), solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "AutoTsit5(Rodas5P)"
        @test meta.initial_solver == "AutoTsit5"
    end
    withenv(
        "SPACEAGORA_SOLVER_MODE" => "auto_stiff",
        "SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5" => "1",
        "SPACEAGORA_SOLVER_MAXITERS" => nothing
    ) do
        sol, meta = _solve_with_solver_policy(prob_simple, _active_solver_config(), solver_args_atmo, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "AutoTsit5(Rodas5P)"
        @test meta.initial_solver == "AutoTsit5"
    end

    u_alloc_probe = build_initial_conditions(solver_args)
    p_alloc_probe = ODEParams(n_sats=1, args=solver_args)
    sc_alloc_probe = u_alloc_probe.sc[1]
    gravity_alloc_model = InverseSquaredGravityModel()
    j2_alloc_model = InverseSquaredJ2GravityModel()
    for _ in 1:1000
        calcForceTorque(gravity_alloc_model, sc_alloc_probe, p_alloc_probe, 1)
        calcForceTorque(j2_alloc_model, sc_alloc_probe, p_alloc_probe, 1)
        SimulationEngine._extract_sample_pos_vel(sc_alloc_probe)
    end
    @test (@allocated calcForceTorque(gravity_alloc_model, sc_alloc_probe, p_alloc_probe, 1)) == 0
    @test (@allocated calcForceTorque(j2_alloc_model, sc_alloc_probe, p_alloc_probe, 1)) == 0
    @test (@allocated SimulationEngine._extract_sample_pos_vel(sc_alloc_probe)) == 0

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
        sol, meta = _solve_with_solver_policy(split_prob_simple, _active_solver_config(), solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "KenCarp4(IMEX)"
        @test meta.initial_solver == "KenCarp4"
    end
    withenv(
        "SPACEAGORA_SOLVER_MODE" => "split_imex",
        "SPACEAGORA_SPLIT_IMEX_SOLVER" => "kencarp47",
        "SPACEAGORA_SOLVER_MAXITERS" => nothing
    ) do
        sol, meta = _solve_with_solver_policy(split_prob_simple, _active_solver_config(), solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "KenCarp47(IMEX)"
        @test meta.initial_solver == "KenCarp47"
    end

    args_split_gravity = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    u0_split_gravity = build_initial_conditions(args_split_gravity)
    p_split_gravity = ODEParams(n_sats=1, args=args_split_gravity)
    _initialize_heat_rate_buffers!(p_split_gravity)

    withenv("SPACEAGORA_SOLVER_MODE" => "split_imex") do
        prob_split = _build_typed_solver_problem(u0_split_gravity, (0.0, 1.0), p_split_gravity, CallbackSet(), _solver_policy_mode())
        @test hasproperty(prob_split.f, :f1)
        @test hasproperty(prob_split.f, :f2)
        du_implicit = copy(u0_split_gravity)
        du_explicit = copy(u0_split_gravity)
        du_implicit .= 0.0
        du_explicit .= 0.0
        prob_split.f.f1(du_implicit, u0_split_gravity, p_split_gravity, 0.0)
        prob_split.f.f2(du_explicit, u0_split_gravity, p_split_gravity, 0.0)
        @test all(==(0.0), du_implicit.sc[1].pos)
        @test all(==(0.0), du_implicit.sc[1].vel)
        @test du_implicit.sc[1].mass == 0.0
        @test isapprox(
            SVector{3, Float64}(du_explicit.sc[1].pos),
            SVector{3, Float64}(u0_split_gravity.sc[1].vel);
            atol=1e-12,
            rtol=1e-10
        )
        @test norm(SVector{3, Float64}(du_explicit.sc[1].vel)) > 0.0
    end

    withenv("SPACEAGORA_SOLVER_MODE" => "multirate") do
        prob_multirate = _build_typed_solver_problem(u0_split_gravity, (0.0, 1.0), p_split_gravity, CallbackSet(), _solver_policy_mode())
        @test hasproperty(prob_multirate.f, :f1)
        @test hasproperty(prob_multirate.f, :f2)
        du_slow = copy(u0_split_gravity)
        du_fast = copy(u0_split_gravity)
        du_slow .= 0.0
        du_fast .= 0.0
        prob_multirate.f.f1(du_slow, u0_split_gravity, p_split_gravity, 0.0)
        prob_multirate.f.f2(du_fast, u0_split_gravity, p_split_gravity, 0.0)
        @test isapprox(
            SVector{3, Float64}(du_slow.sc[1].pos),
            SVector{3, Float64}(u0_split_gravity.sc[1].vel);
            atol=1e-12,
            rtol=1e-10
        )
        @test norm(SVector{3, Float64}(du_slow.sc[1].vel)) > 0.0
        @test all(==(0.0), du_fast.sc[1].pos)
        @test all(==(0.0), du_fast.sc[1].vel)
        @test du_fast.sc[1].mass == 0.0
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

    args_symplectic = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel()
    )
    @test _symplectic_conservative_eligible(args_symplectic) == true
    @test _symplectic_conservative_eligible(args_eff_single) == false
    @test SimulationModel.gravity_backbone_structure(ConstantGravityModel()) == :position_only_static_gravity
    @test SimulationModel.gravity_backbone_structure(InverseSquaredGravityModel()) == :position_only_static_gravity
    @test SimulationModel.gravity_backbone_structure(InverseSquaredJ2GravityModel()) == :position_only_static_gravity
    @test _gravity_backbone_structure_validated(InverseSquaredGravityModel()) == :position_only_static_gravity
    @test_throws ArgumentError _gravity_backbone_structure_validated(InvalidBackboneStructureModel())
    @test SimulationModel.gravity_backbone_kick_structure(SolarRadiationPressureModel(1.2, 12.0)) == :velocity_kick_explicit
    @test SimulationModel.gravity_backbone_kick_structure(NBodyGravityModel(body_names=("moon",), primary_body_name="Earth")) == :velocity_kick_explicit
    @test _gravity_backbone_kick_structure_validated(SolarRadiationPressureModel(1.2, 12.0)) == :velocity_kick_explicit
    @test_throws ArgumentError _gravity_backbone_kick_structure_validated(InvalidBackboneKickStructureModel())

    kepler_second_order! = function (ddu, du, u, p, t)
        r = norm(u)
        @. ddu = -p.μ * u / r^3
    end
    r0_sym = [EARTH.Rp_e + 500e3, 0.0, 0.0]
    v0_sym = [0.0, sqrt(EARTH.μ / r0_sym[1]), 0.0]
    prob_symplectic = SecondOrderODEProblem(kepler_second_order!, v0_sym, r0_sym, (0.0, 120.0), (μ=EARTH.μ,))

    withenv(
        "SPACEAGORA_SOLVER_MODE" => "symplectic",
        "SPACEAGORA_SYMPLECTIC_DT_S" => "2.0",
        "SPACEAGORA_SOLVER_MAXITERS" => nothing
    ) do
        sol, meta = _solve_with_solver_policy(prob_symplectic, _active_solver_config(), args_symplectic, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "KahanLi8(Symplectic)"
        @test meta.initial_solver == "KahanLi8"
        @test meta.fallback_used == false
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "symplectic", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        @test_throws ArgumentError _solve_with_solver_policy(prob_simple, _active_solver_config(), args_symplectic, 1e-8, 1e-8)
    end

    args_backbone = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel()
    )
    @test _gravity_backbone_eligible(args_backbone) == true
    @test isnothing(_gravity_backbone_reject_reason(args_backbone))

    args_backbone_orient = build_config(
        spacecraft=make_spacecraft(
            ra_alt_m=500e3,
            rp_alt_m=500e3,
            orientation_state=(normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9)), SVector{3, Float64}(0.0, 0.0, 0.0))
        ),
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel()
    )
    @test occursin("orientation_sim=false", _gravity_backbone_reject_reason(args_backbone_orient))

    args_backbone_control = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(TimedTangentialThrusterModel(1.0, 1.0, 0.0, 10.0),),
        control_rates=[1.0],
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel()
    )
    @test occursin("no control effectors", _gravity_backbone_reject_reason(args_backbone_control))

    args_backbone_guidance = SimulationConfiguration(
        simulation_settings=args_backbone.simulation_settings,
        mission_configuration=args_backbone.mission_configuration,
        environment_model=args_backbone.environment_model,
        dynamics_model=args_backbone.dynamics_model,
        guidance_model=GuidanceModel(guidance_effectors=(CountingGuidanceModel([0]),), guidance_rates=[1.0]),
        navigation_model=args_backbone.navigation_model,
        control_model=args_backbone.control_model,
        initial_time=args_backbone.initial_time,
        integration_tolerances=args_backbone.integration_tolerances
    )
    @test occursin("no guidance effectors", _gravity_backbone_reject_reason(args_backbone_guidance))

    args_backbone_navigation = SimulationConfiguration(
        simulation_settings=args_backbone.simulation_settings,
        mission_configuration=args_backbone.mission_configuration,
        environment_model=args_backbone.environment_model,
        dynamics_model=args_backbone.dynamics_model,
        guidance_model=args_backbone.guidance_model,
        navigation_model=NavigationModel(navigation_effectors=(CountingNavigationModel([0]),), navigation_rates=[1.0]),
        control_model=args_backbone.control_model,
        initial_time=args_backbone.initial_time,
        integration_tolerances=args_backbone.integration_tolerances
    )
    @test occursin("no navigation effectors", _gravity_backbone_reject_reason(args_backbone_navigation))

    args_backbone_bad_core_solar = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(SolarBackboneModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel()
    )
    @test occursin("solar/SRP-dependent gravity core", _gravity_backbone_reject_reason(args_backbone_bad_core_solar))

    args_backbone_srp = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), SolarRadiationPressureModel(1.2, 12.0; direct=true, albedo=true, ir=true)),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel()
    )
    @test _gravity_backbone_eligible(args_backbone_srp) == true

    args_backbone_nbody = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), NBodyGravityModel(body_names=("moon",), primary_body_name="Earth")),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel()
    )
    @test _gravity_backbone_eligible(args_backbone_nbody) == true

    args_backbone_srp_nbody = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            NBodyGravityModel(body_names=("moon",), primary_body_name="Earth"),
            SolarRadiationPressureModel(1.2, 12.0; direct=true, albedo=true, ir=true),
        ),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel()
    )
    @test _gravity_backbone_eligible(args_backbone_srp_nbody) == true

    args_backbone_bad_kick_capability = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), PlanetFrameKickModel()),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel()
    )
    @test occursin("planet-frame-dependent perturbation kicks", _gravity_backbone_reject_reason(args_backbone_bad_kick_capability))

    args_backbone_custom = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=2.0,
        EI_km=120.0,
        dynamic_effectors=(BackboneCustomGravityModel(SVector{3, Float64}(1.0e-6, 0.0, 0.0)),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel()
    )
    @test _gravity_backbone_eligible(args_backbone_custom) == true

    u0_backbone = build_initial_conditions(args_backbone)
    p_backbone = ODEParams(n_sats=1, args=args_backbone)
    _initialize_heat_rate_buffers!(p_backbone)
    callbacks_backbone = CallbackSet()
    withenv("SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split") do
        prob_backbone = _build_typed_solver_problem(u0_backbone, (0.0, 30.0), p_backbone, callbacks_backbone, _solver_policy_mode())
        @test getproperty(prob_backbone, :problem_type) isa SecondOrderODEProblem
        q0_backbone, dq0_backbone = _gravity_backbone_initial_states(u0_backbone, args_backbone)
        ddu_backbone = copy(dq0_backbone)
        ddu_backbone .= 0.0
        spacecraft_dynamics_gravity_backbone!(ddu_backbone, dq0_backbone, q0_backbone, p_backbone, 0.0)
        @test norm(SVector{3, Float64}(ddu_backbone.sc[1].vel)) > 0.0
    end

    withenv(
        "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
        "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "2.0",
        "SPACEAGORA_SOLVER_MAXITERS" => nothing
    ) do
        prob_backbone = _build_typed_solver_problem(u0_backbone, (0.0, 30.0), p_backbone, callbacks_backbone, _solver_policy_mode())
        sol_backbone, meta_backbone = _solve_with_solver_policy(prob_backbone, _active_solver_config(), args_backbone, 1e-8, 1e-8)
        @test string(sol_backbone.retcode) == "Success"
        @test meta_backbone.solver == "KahanLi8(GravityBackbone)"
        @test _is_gravity_backbone_state(sol_backbone.u[end])
    end

    args_backbone_kicks = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=8.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            NBodyGravityModel(body_names=("moon",), primary_body_name="Earth"),
            SolarRadiationPressureModel(1.2, 12.0; direct=true, albedo=true, ir=true),
        ),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel()
    )
    u0_backbone_kicks = build_initial_conditions(args_backbone_kicks)
    p_backbone_kicks = ODEParams(n_sats=1, args=args_backbone_kicks)
    _initialize_heat_rate_buffers!(p_backbone_kicks)
    _initialize_nbody_ephemeris_cache_buffer!(p_backbone_kicks)
    _initialize_srp_sun_cache_buffer!(p_backbone_kicks)
    _initialize_planet_frame_cache_buffer!(p_backbone_kicks)
    _initialize_spice_rhs_memo_mode!(p_backbone_kicks)
    _reset_spice_rhs_memo!(p_backbone_kicks)
    _reset_spice_runtime_counters!(p_backbone_kicks)
    p_backbone_kicks.shared_buffers.et_start[] = SimulationModel.ephemerides_time_seconds(
        args_backbone_kicks.initial_time,
        args_backbone_kicks.environment_model.ephemerides_model,
    )
    _initialize_nbody_ephemeris_cache!(p_backbone_kicks, p_backbone_kicks.shared_buffers.et_start[], args_backbone_kicks.mission_configuration.mission_time)
    _initialize_srp_sun_ephemeris_cache!(p_backbone_kicks, p_backbone_kicks.shared_buffers.et_start[], args_backbone_kicks.mission_configuration.mission_time)

    withenv(
        "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
        "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "2.0",
        "SPACEAGORA_SOLVER_MAXITERS" => nothing
    ) do
        callbacks_backbone_kicks = CallbackSet()
        prob_backbone_kicks = _build_typed_solver_problem(u0_backbone_kicks, (0.0, 2.0), p_backbone_kicks, callbacks_backbone_kicks, _solver_policy_mode())
        sol_backbone_kicks, meta_backbone_kicks = _solve_with_solver_policy(prob_backbone_kicks, _active_solver_config(), args_backbone_kicks, 1e-8, 1e-8)
        @test string(sol_backbone_kicks.retcode) == "Success"
        @test meta_backbone_kicks.solver == "KahanLi8(GravityBackbone+Kicks)"
        @test length(sol_backbone_kicks.t) == 2
        @test _is_gravity_backbone_state(sol_backbone_kicks.u[end])

        manual_state = deepcopy(prob_backbone_kicks.u0)
        _gravity_backbone_half_kick!(manual_state, p_backbone_kicks, 0.0, 1.0)
        core_prob = remake(prob_backbone_kicks; u0=manual_state, tspan=(0.0, 2.0))
        core_sol = solve(core_prob, KahanLi8(); dt=2.0)
        @test string(core_sol.retcode) == "Success"
        manual_state = deepcopy(core_sol.u[end])
        _gravity_backbone_half_kick!(manual_state, p_backbone_kicks, 2.0, 1.0)

        final_state = sol_backbone_kicks.u[end]
        @test _state_position_ii(final_state, 1) ≈ _state_position_ii(manual_state, 1) atol=1e-8 rtol=1e-8
        @test _state_velocity_ii(final_state, 1) ≈ _state_velocity_ii(manual_state, 1) atol=1e-10 rtol=1e-8
    end

    args_backbone_kicks_run = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=6.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            NBodyGravityModel(["moon"], "Earth", SPICE_PATH),
            SolarRadiationPressureModel(1.2, 12.0),
        ),
        keplerian=true,
        ephemerides_model=SpiceEphemeridesModel(),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, save_csv=false, normalize=false)
    )

    withenv(
        "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
        "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "2.0"
    ) do
        backbone_result = run_simulation(args_backbone_kicks_run; return_solution=true, return_solver_metadata=true)
        @test backbone_result.solution !== nothing
        @test backbone_result.solver_mode == "gravity_backbone_split"
        @test length(backbone_result.solution.t) >= 2
        @test _is_gravity_backbone_state(backbone_result.solution.u[end])
        @test backbone_result.solver_trace[end].solver == "KahanLi8(GravityBackbone+Kicks)"
        @test backbone_result.solution.t[end] ≈ 6.0 atol=1e-8 rtol=0.0
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
        sol, meta = _solve_with_solver_policy(split_prob_simple, _active_solver_config(), solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test occursin("Multirate(Strang;", meta.solver)
        @test meta.initial_solver == "Tsit5"
    end

    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES" => "oops") do
        @test_throws ArgumentError _ephemeris_reuse_max_entries()
    end

    reuse_cache = Dict{Symbol, SRPSunEphemerisCache}()
    reuse_value_a = SRPSunEphemerisCache([0.0], SVector{3, Float64}[SVector{3, Float64}(1.0, 0.0, 0.0)])
    reuse_value_b = SRPSunEphemerisCache([0.0], SVector{3, Float64}[SVector{3, Float64}(2.0, 0.0, 0.0)])
    @test _ephemeris_reuse_store!(reuse_cache, :k1, reuse_value_a, 0) === reuse_value_a
    @test !haskey(reuse_cache, :k1)
    _ephemeris_reuse_store!(reuse_cache, :k1, reuse_value_a, 2)
    @test _ephemeris_reuse_store!(reuse_cache, :k1, reuse_value_b, 2) === reuse_value_a
    _ephemeris_reuse_store!(reuse_cache, :k2, reuse_value_b, 1)
    @test !haskey(reuse_cache, :k1)
    @test haskey(reuse_cache, :k2)

    @test SimulationEngine._srp_ephemeris_reuse_key("earth", 0.0, 10.0, 1.0) isa SimulationEngine.SRPEphemerisReuseKey
    @test SimulationEngine._nbody_ephemeris_reuse_key("earth", ["moon", "sun"], 0.0, 10.0, 1.0) isa SimulationEngine.NBodyEphemerisReuseKey
    @test SimulationEngine._planet_frame_ephemeris_reuse_key(EARTH, SpiceEphemeridesModel(), 0.0, 10.0, 1.0) isa SimulationEngine.PlanetFrameEphemerisReuseKey
    @test _wrench_method_available(WrenchOnlyForceModel(SVector{3, Float64}(1.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)))
    @test !_wrench_method_available(ConstantForceModel(SVector{3, Float64}(1.0, 0.0, 0.0)))
    @test SimulationModel.solver_partition(AerodynamicCoefficientConstant()) == :implicit
    @test SimulationModel.solver_partition(AerodynamicCoefficientfM()) == :implicit
    @test SimulationModel.solver_partition(AerodynamicCoefficientNoBallisticFlight()) == :implicit
    @test SimulationModel.solver_partition(InverseSquaredGravityModel()) == :explicit
    @test SimulationModel.solver_partition(SolarRadiationPressureModel(1.2, 12.0, 500.0)) == :explicit
    @test SimulationModel.solver_partition(NBodyGravityModel(body_names=("moon",), primary_body_name="Earth")) == :explicit
    @test _solver_partition_validated(AerodynamicCoefficientfM()) == :implicit
    @test _solver_partition_validated(ConstantForceModel(SVector{3, Float64}(1.0, 0.0, 0.0))) == :explicit
    @test_throws ArgumentError _solver_partition_validated(InvalidPartitionForceModel())

    q0_wrench = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    ω0_wrench = SVector{3, Float64}(0.0, 0.0, 0.0)
    sc_wrench = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0,
        orientation_state=(q0_wrench, ω0_wrench)
    )
    sc_wrench.inertia_tensor = SMatrix{3, 3, Float64}(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0)
    legacy_force = SVector{3, Float64}(1.0, 2.0, 3.0)
    typed_force = SVector{3, Float64}(4.0, 5.0, 6.0)
    typed_torque = SVector{3, Float64}(0.2, -0.1, 0.3)
    args_wrench_mix = build_config(
        spacecraft=sc_wrench,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(
            ConstantForceModel(legacy_force),
            WrenchOnlyForceModel(typed_force, typed_torque),
        ),
        keplerian=true
    )
    u0_wrench_mix = build_initial_conditions(args_wrench_mix)
    du0_wrench_mix = copy(u0_wrench_mix)
    du0_wrench_mix .= 0.0
    p_wrench_mix = ODEParams(n_sats=1, args=args_wrench_mix)
    spacecraft_dynamics!(du0_wrench_mix, u0_wrench_mix, p_wrench_mix, 0.0)
    expected_force_mix = legacy_force + typed_force
    @test isapprox(SVector{3, Float64}(du0_wrench_mix.sc[1].vel), expected_force_mix / u0_wrench_mix.sc[1].mass; atol=1e-12, rtol=1e-10)
    @test isapprox(
        SVector{3, Float64}(du0_wrench_mix.sc[1].ω),
        args_wrench_mix.dynamics_model.spacecraft[1].inertia_tensor \ typed_torque;
        atol=1e-12,
        rtol=1e-10
    )

    env_empty = sample_environment(
        EffectorEnvironmentRequirements(),
        WrenchOnlyForceModel(typed_force, typed_torque),
        u0_wrench_mix.sc[1],
        p_wrench_mix,
        1,
        0.0;
        write_buffers=false,
    )
    @test env_empty.planet === EARTH
    @test env_empty.planet_frame === nothing
    @test env_empty.atmosphere === nothing
    @test env_empty.solar === nothing
    @test env_empty.third_bodies === nothing

    args_probe = build_config(
        spacecraft=make_spacecraft(ra_alt_m=300e3, rp_alt_m=300e3, ν_deg=180.0),
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(AtmosphereProbeWrenchModel(),),
        keplerian=true
    )
    u0_probe = build_initial_conditions(args_probe)
    p_probe = ODEParams(n_sats=1, args=args_probe)
    state_probe = build_state_sample(u0_probe.sc[1], args_probe.dynamics_model.spacecraft[1], false)
    req_probe = environment_requirements(AtmosphereProbeWrenchModel())
    env_probe = sample_environment(req_probe, AtmosphereProbeWrenchModel(), u0_probe.sc[1], p_probe, 1, 0.0; write_buffers=false)
    force_probe, torque_probe = _evaluate_dynamic_effector(AtmosphereProbeWrenchModel(), u0_probe.sc[1], state_probe, p_probe, 1, 0.0)
    @test force_probe == SVector{3, Float64}(env_probe.atmosphere.rho_kg_m3, env_probe.planet_frame.alt_m, 0.0)
    @test torque_probe == SVector{3, Float64}(0.0, 0.0, 0.0)

    SimulationEngine._ensure_rhs_flat_effector_scratch!(p_probe.shared_buffers, 1, 1)
    SimulationEngine._prefill_environment_samples!(p_probe, 0.0, u0_probe.sc; atmosphere=true)
    p_probe.shared_buffers.rhs_planet_frame_prefilled[] = true
    p_probe.shared_buffers.rhs_atmosphere_prefilled[] = true
    env_probe_reused = SimulationEngine.sample_environment_with_reusable_buffers(
        req_probe,
        AtmosphereProbeWrenchModel(),
        u0_probe.sc[1],
        p_probe,
        1,
        0.0,
    )
    @test env_probe_reused.planet_frame == SimulationEngine.sample_buffered_planet_frame(p_probe, 1)
    @test env_probe_reused.atmosphere == SimulationEngine.sample_buffered_atmosphere(u0_probe.sc[1], p_probe, 1, 0.0)
    p_probe.shared_buffers.rhs_planet_frame_prefilled[] = false
    p_probe.shared_buffers.rhs_atmosphere_prefilled[] = false

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
    @test _has_active_srp_effector((SolarRadiationPressureModel(1.2, 10.0; direct=false, albedo=true),))
    @test !_has_active_srp_effector((SolarRadiationPressureModel(1.2, 10.0; direct=false, albedo=false, ir=true),))

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

    p_workspace_resize = ODEParams(n_sats=1, args=args_eff_single)
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

    gram_model_instances = SimulationModel.GRAMAtmosphereModel(Ref{Any}(nothing))
    args_density_instances = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=gram_model_instances,
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    p_density_instances = ODEParams(n_sats=1, args=args_density_instances)
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
    p_density_surrogate = ODEParams(n_sats=1, args=args_density_surrogate)
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
    p_srp = ODEParams(n_sats=1, args=args_srp)
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
    p_srp_reuse_a = ODEParams(n_sats=1, args=args_srp)
    p_srp_reuse_b = ODEParams(n_sats=1, args=args_srp)
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
        primary_body_name = SimulationModel.DynamicEffectors._spice_query_name(args_srp.environment_model.planet.name)
        raw_sun_j2000_km = lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
            SVector{3, Float64}(spkpos("sun", cache_a.ets[1], "J2000", "none", primary_body_name)[1])
        end
        @test cache_a.positions_j2000_m[1] == raw_sun_j2000_km * 1.0e3
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
    p_nbody = ODEParams(n_sats=1, args=args_nbody)
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
    _reset_spice_runtime_counters!(p_nbody)
    withenv(
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE" => "1",
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S" => "10.0",
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES" => "100"
    ) do
        _initialize_nbody_ephemeris_cache!(p_nbody, 0.0, 10.0)
        cache = p_nbody.shared_buffers.nbody_ephemeris_cache[]
        @test cache isa NBodyEphemerisCache
        @test p_nbody.shared_buffers.spice_runtime_counters.nbody_spkpos_cache_build_calls[] == 2
        raw_moon_j2000_km = lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
            SVector{3, Float64}(spkpos(cache.body_query_names[1], cache.ets[1], "J2000", "none", cache.primary_body_name)[1])
        end
        @test cache.positions_j2000_m[1, 1] == raw_moon_j2000_km * 1.0e3
    end

    et_start_nbody = SimulationModel.ephemerides_time_seconds(args_nbody.initial_time, args_nbody.environment_model.ephemerides_model)
    _clear_ephemeris_reuse_cache!()
    withenv(
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S" => "10.0",
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES" => "100",
        "SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "0"
    ) do
        cache_prewarmed = SpaceAGORA.prewarm_nbody_ephemeris_cache(args_nbody; dt_s=10.0, mission_end_s=10.0)
        @test cache_prewarmed.primary_body_name == "earth"
        @test cache_prewarmed.body_query_names == ["moon"]
        @test length(cache_prewarmed.ets) == 2
        @test size(cache_prewarmed.positions_j2000_m) == (2, 1)
        @test cache_prewarmed.ets[1] == et_start_nbody
    end
    _clear_ephemeris_reuse_cache!()

    mktempdir() do tmp
        cache_path = joinpath(tmp, "nbody_ephemeris_cache.bin")
        withenv(
            "SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S" => "10.0",
            "SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES" => "100",
            "SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "0"
        ) do
            cache_saved = SpaceAGORA.prewarm_nbody_ephemeris_cache(args_nbody; dt_s=10.0, mission_end_s=10.0, save_path=cache_path)
            @test isfile(cache_path)
            _clear_ephemeris_reuse_cache!()

            cache_loaded = SpaceAGORA.load_nbody_ephemeris_cache!(cache_path)
            @test cache_loaded.primary_body_name == cache_saved.primary_body_name
            @test cache_loaded.body_query_names == cache_saved.body_query_names
            @test cache_loaded.ets == cache_saved.ets
            @test cache_loaded.positions_j2000_m == cache_saved.positions_j2000_m
        end

        bad_cache_path = joinpath(tmp, "bad_nbody_ephemeris_cache.bin")
        open(bad_cache_path, "w") do io
            serialize(io, Dict(:invalid => true))
        end
        @test_throws ArgumentError SpaceAGORA.load_nbody_ephemeris_cache!(bad_cache_path)
    end
    _clear_ephemeris_reuse_cache!()

    p_planet_frame = ODEParams(n_sats=1, args=args_srp)
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
    p_nbody_srp = ODEParams(n_sats=1, args=args_nbody_srp)
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
    @test harmonics_model.coefficient_normalization == :full
    @test SimulationModel.solver_partition(harmonics_model) == :explicit

    @testset "Harmonics Normalization Conversion And Reference Parity" begin
        mktempdir() do tmp
            max_degree = 4
            max_order = 4
            base_coeffs = Dict(
                (2, 0) => (-4.8416945732e-4, 0.0),
                (2, 1) => (-3.1034310672e-10, 1.4107575094e-9),
                (2, 2) => (2.4393734159e-6, -1.4002940118e-6),
                (3, 0) => (9.5716475834e-7, 0.0),
                (3, 1) => (2.0304620109e-6, 2.4821193046e-7),
                (4, 0) => (5.3996586664e-7, 0.0),
                (4, 2) => (3.5098919076e-7, -1.8722151127e-7)
            )
            coeff_at(l, m) = get(base_coeffs, (l, m), (0.0, 0.0))

            coeff_full_path = joinpath(tmp, "harmonics_full.csv")
            coeff_schmidt_path = joinpath(tmp, "harmonics_schmidt.csv")
            coeff_unnorm_path = joinpath(tmp, "harmonics_unnorm.csv")

            write_gravity_harmonics_csv(coeff_full_path, max_degree, max_order) do l, m
                coeff_at(l, m)
            end
            write_gravity_harmonics_csv(coeff_schmidt_path, max_degree, max_order) do l, m
                c_full, s_full = coeff_at(l, m)
                scale = sqrt(2.0 * l + 1.0)
                return c_full * scale, s_full * scale
            end
            write_gravity_harmonics_csv(coeff_unnorm_path, max_degree, max_order) do l, m
                c_full, s_full = coeff_at(l, m)
                scale = harmonics_full_normalization_factor(l, m)
                return c_full * scale, s_full * scale
            end

            model_full = GravitationalHarmonicsModel(max_degree, max_order, coeff_full_path, EARTH; coefficient_normalization=:full)
            model_full_alias = GravitationalHarmonicsModel(max_degree, max_order, coeff_full_path, EARTH; coefficient_normalization=:fully_normalized)
            model_schmidt = GravitationalHarmonicsModel(max_degree, max_order, coeff_schmidt_path, EARTH; coefficient_normalization=:schmidt)
            model_unnorm = GravitationalHarmonicsModel(max_degree, max_order, coeff_unnorm_path, EARTH; coefficient_normalization=:unnormalized)

            @test model_full.coefficient_normalization == :full
            @test model_full_alias.coefficient_normalization == :full
            @test model_schmidt.coefficient_normalization == :schmidt
            @test model_unnorm.coefficient_normalization == :unnormalized
            @test_throws ArgumentError GravitationalHarmonicsModel(max_degree, max_order, coeff_full_path, EARTH; coefficient_normalization=:bogus)
            @test isapprox(model_full_alias.C, model_full.C; atol=0.0, rtol=0.0)
            @test isapprox(model_full_alias.S, model_full.S; atol=0.0, rtol=0.0)
            @test isapprox(model_schmidt.C, model_full.C; atol=1e-15, rtol=0.0)
            @test isapprox(model_schmidt.S, model_full.S; atol=1e-15, rtol=0.0)
            @test isapprox(model_unnorm.C, model_full.C; atol=1e-15, rtol=0.0)
            @test isapprox(model_unnorm.S, model_full.S; atol=1e-15, rtol=0.0)

            coeff_gfc_path = joinpath(tmp, "harmonics_reference.gfc")
            write_icgem_from_harmonics_model(coeff_gfc_path, model_full; model_name="SpaceAGORA_Test_Harmonics")
            ref_model = SatelliteToolboxGravityModels.GravityModels.load(SatelliteToolboxGravityModels.IcgemFile, coeff_gfc_path)
            @test SatelliteToolboxGravityModels.GravityModels.coefficient_norm(ref_model) == :full

            args_harmonics_compare = build_config(
                spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
                density_model=NoAtmosphereModel(),
                orientation_sim=false,
                mission_time=60.0,
                EI_km=120.0,
                dynamic_effectors=(InverseSquaredGravityModel(),),
                keplerian=true,
                ephemerides_model=SimpleEphemeridesModel()
            )
            p_harmonics_compare = ODEParams(n_sats=1, args=args_harmonics_compare)
            _initialize_harmonics_workspace_buffers!(p_harmonics_compare)

            sample_positions_pp = (
                SVector{3, Float64}(EARTH.Rp_e + 450e3, 0.0, 0.0),
                SVector{3, Float64}(4.8e6, 2.1e6, 1.7e6),
                SVector{3, Float64}(-3.9e6, 5.2e6, 2.4e6)
            )
            for pos_pp in sample_positions_pp
                state = ComponentVector(pos=collect(pos_pp), vel=[0.0, 0.0, 0.0], mass=200.0)
                force_pp, torque_pp = calcForceTorque(model_full, state, p_harmonics_compare, 1)
                acc_pp = force_pp / state.mass
                ref_total_pp = SVector{3, Float64}(
                    SatelliteToolboxGravityModels.GravityModels.gravitational_acceleration(
                        ref_model,
                        collect(pos_pp);
                        max_degree=model_full.L,
                        max_order=model_full.M
                    )
                )

                @test torque_pp == SVector{3, Float64}(0.0, 0.0, 0.0)
                @test isapprox(acc_pp, ref_total_pp; atol=1e-12, rtol=1e-9)
            end
        end
    end

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
    @test dyn._nbody_body_position_from_cache_j2000_m(nbody_cache, -1.0, "moon_barycenter", "earth_barycenter") === nothing
    @test dyn._nbody_body_position_from_cache_j2000_m(nbody_cache, NaN, "moon_barycenter", "earth_barycenter") isa SVector{3, Float64}
    @test dyn._nbody_body_position_from_cache_j2000_m(nbody_cache, 15.0, "moon_barycenter", "earth_barycenter") == nbody_positions[4, 1]
    @test dyn._nbody_body_position_from_cache_j2000_m(nbody_cache, 2.5, "moon_barycenter", "earth_barycenter") == SVector{3, Float64}(1.5, 0.0, 0.0)
    nbody_interp_catmull = dyn._nbody_body_position_from_cache_j2000_m(nbody_cache, 7.5, "moon_barycenter", "earth_barycenter")
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
    @test dyn._srp_sun_position_from_cache_j2000_m(srp_cache, -1.0) === nothing
    @test dyn._srp_sun_position_from_cache_j2000_m(srp_cache, NaN) isa SVector{3, Float64}
    @test dyn._srp_sun_position_from_cache_j2000_m(srp_cache, 15.0) == srp_cache.positions_j2000_m[end]
    @test dyn._srp_sun_position_from_cache_j2000_m(srp_cache, 2.5) == SVector{3, Float64}(1.5, 0.0, 0.0)
    srp_interp_catmull = dyn._srp_sun_position_from_cache_j2000_m(srp_cache, 7.5)
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
    annular_a = asin(695000e3 / 1.5e11)
    annular_b = asin(EARTH.Rp_e / 1.0e10)
    @test isapprox(annular_ratio, 1.0 - annular_b^2 / annular_a^2; atol=0.0, rtol=1e-12)
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
        function _shadow_geometry_for_apparent_separation(c_rad::Float64; r_sat_m::Float64=10.0 * EARTH.Rp_e, r_sun_m::Float64=149_597_870_700.0)
            return (
                SVector{3, Float64}(-r_sat_m * cos(c_rad), r_sat_m * sin(c_rad), 0.0),
                SVector{3, Float64}(r_sun_m, 0.0, 0.0),
            )
        end

        # Direct 1 AU cannonball case: force magnitude should match P0 * Cr * A when the solver
        # evaluates the SRP effector with a live spacecraft mass.
        p_srp.shared_buffers.srp_sun_ephemeris_cache[] = SRPSunEphemerisCache(
            [0.0, 10.0],
            SVector{3, Float64}[
                SVector{3, Float64}(149_604_870_700.0, 100.0, 50.0),
                SVector{3, Float64}(149_604_870_700.0, 100.0, 50.0)
            ]
        )
        p_srp.shared_buffers.et_start[] = 0.0
        p_srp.shared_buffers.current_time[] = 5.0
        x_srp_regression = Float64[7_000_000.0, 100.0, 50.0, 0.0, 0.0, 0.0, 200.0]

        force_on, _ = calcForceTorque(SolarRadiationPressureModel(1.2, 12.0), x_srp_regression, p_srp, 1)
        force_off, _ = calcForceTorque(SolarRadiationPressureModel(1.2, 0.0), x_srp_regression, p_srp, 1)
        expected_force_mag = 4.56e-6 * 1.2 * 12.0
        @test norm(force_on) > 0.0
        @test force_off == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test norm(force_on - force_off) > 0.0
        @test isapprox(norm(force_on), expected_force_mag; atol=0.0, rtol=1e-12)

        # Eclipse attenuation: partial eclipse force norm must scale by eclipse ratio.
        pos_partial = SVector{3, Float64}(6.930027129188876e6, -2.6352977471555886e6, -3.363004422597388e6)
        sun_partial = SVector{3, Float64}(-9.438128429326639e10, -6.696979722657439e10, 1.7822072441008075e11)
        partial_ratio_regression = dyn.eclipse_area_calc(pos_partial, sun_partial, EARTH.Rp_e)
        accel_partial = dyn.srp_cannonball_accel(pos_partial, sun_partial, EARTH.Rp_e, 4.56e-6, 1.2, 12.0, 200.0)
        accel_no_eclipse = dyn.srp_cannonball_accel(pos_partial, sun_partial, 0.0, 4.56e-6, 1.2, 12.0, 200.0)
        @test 0.0 < partial_ratio_regression < 1.0
        @test norm(accel_no_eclipse) > 0.0
        @test isapprox(norm(accel_partial), partial_ratio_regression * norm(accel_no_eclipse); atol=0.0, rtol=1e-12)

        # Eclipse transition geometry should move monotonically from total -> partial -> none.
        transition_r_sat = 10.0 * EARTH.Rp_e
        transition_r_sun = 149_597_870_700.0
        apparent_sun = asin(695000e3 / transition_r_sun)
        apparent_planet = asin(EARTH.Rp_e / transition_r_sat)
        transition_eps = 1.0e-4
        pos_total, sun_total = _shadow_geometry_for_apparent_separation(apparent_planet - apparent_sun - transition_eps; r_sat_m=transition_r_sat, r_sun_m=transition_r_sun)
        pos_just_partial, sun_just_partial = _shadow_geometry_for_apparent_separation(apparent_planet - apparent_sun + transition_eps; r_sat_m=transition_r_sat, r_sun_m=transition_r_sun)
        pos_just_shadowed, sun_just_shadowed = _shadow_geometry_for_apparent_separation(apparent_planet + apparent_sun - transition_eps; r_sat_m=transition_r_sat, r_sun_m=transition_r_sun)
        pos_none, sun_none = _shadow_geometry_for_apparent_separation(apparent_planet + apparent_sun + transition_eps; r_sat_m=transition_r_sat, r_sun_m=transition_r_sun)
        ratio_total = dyn.eclipse_area_calc(pos_total, sun_total, EARTH.Rp_e)
        ratio_just_partial = dyn.eclipse_area_calc(pos_just_partial, sun_just_partial, EARTH.Rp_e)
        ratio_just_shadowed = dyn.eclipse_area_calc(pos_just_shadowed, sun_just_shadowed, EARTH.Rp_e)
        ratio_none_transition = dyn.eclipse_area_calc(pos_none, sun_none, EARTH.Rp_e)
        @test ratio_total == 0.0
        @test 0.0 < ratio_just_partial < 1.0
        @test 0.0 < ratio_just_shadowed < 1.0
        @test isapprox(ratio_none_transition, 1.0; atol=1e-12, rtol=0.0)
        @test ratio_just_partial < ratio_just_shadowed

        # Planetary albedo and IR remain separate additive terms. Albedo should vanish on the
        # nightside while IR remains available there.
        pos_dayside = SVector{3, Float64}(8.0e6, 0.0, 0.0)
        pos_nightside = SVector{3, Float64}(-8.0e6, 0.0, 0.0)
        sun_dayside = SVector{3, Float64}(149_597_870_700.0, 0.0, 0.0)
        accel_direct_dayside = dyn.srp_cannonball_accel(pos_dayside, sun_dayside, EARTH.Rp_e, 4.56e-6, 1.2, 12.0, 200.0)
        accel_albedo_dayside = dyn.planetary_albedo_accel(pos_dayside, sun_dayside, EARTH.Rp_e, 4.56e-6, 1.2, 12.0, 200.0, 0.3)
        accel_albedo_nightside = dyn.planetary_albedo_accel(pos_nightside, sun_dayside, EARTH.Rp_e, 4.56e-6, 1.2, 12.0, 200.0, 0.3)
        accel_ir_nightside = dyn.planetary_ir_accel(pos_nightside, EARTH.Rp_e, 1.2, 12.0, 200.0, 237.0)
        @test norm(accel_direct_dayside) > 0.0
        @test norm(accel_albedo_dayside) > 0.0
        @test isapprox(norm(accel_albedo_nightside), 0.0; atol=1e-20, rtol=0.0)
        @test norm(accel_ir_nightside) > 0.0

        p_srp.shared_buffers.srp_sun_ephemeris_cache[] = SRPSunEphemerisCache(
            [0.0, 10.0],
            SVector{3, Float64}[sun_dayside, sun_dayside]
        )
        x_dayside = Float64[8.0e6, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0]
        force_total, _ = calcForceTorque(
            SolarRadiationPressureModel(1.2, 12.0; direct=true, albedo=true, ir=true, planet_albedo=0.3, planet_ir_flux_w_m2=237.0),
            x_dayside,
            p_srp,
            1
        )
        expected_total_force = 200.0 * (accel_direct_dayside + accel_albedo_dayside + dyn.planetary_ir_accel(pos_dayside, EARTH.Rp_e, 1.2, 12.0, 200.0, 237.0))
        @test isapprox(force_total, expected_total_force; atol=1e-12, rtol=1e-12)

        x_nightside = Float64[-8.0e6, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0]
        force_ir_only, _ = calcForceTorque(
            SolarRadiationPressureModel(1.2, 12.0; direct=false, albedo=false, ir=true, planet_ir_flux_w_m2=237.0),
            x_nightside,
            p_srp,
            1
        )
        @test isapprox(force_ir_only, 200.0 * accel_ir_nightside; atol=1e-12, rtol=1e-12)

        # Solver parity: legacy solver helper (`srp`) must match canonical SRP kernel for same geometry.
        planet = args_srp.environment_model.planet
        pos_solver = SVector{3, Float64}(x_a[1:3])
        et_solver = 123.0
        primary_body_name = dyn._spice_query_name(planet.name)
        sun_j2000_m = SimulationModel.EphemeridesModels.spice_position_j2000_m("sun", et_solver, primary_body_name)
        accel_kernel = dyn.srp_cannonball_accel(pos_solver, sun_j2000_m, planet.Rp_e, 4.56e-6, 1.2, 12.0, 200.0)
        accel_solver = dyn.srp(planet, 4.56e-6, 1.2, 12.0, 200.0, pos_solver, et_solver)
        @test isapprox(accel_solver, accel_kernel; atol=1e-12, rtol=1e-10)
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
        density_model=ConstantDensityModel(1.0e-9, 220.0),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    u_heat_copy = build_initial_conditions(args_heat_copy)
    du_heat_copy = copy(u_heat_copy)
    du_heat_copy .= 0.0
    p_heat_copy = ODEParams(n_sats=1, args=args_heat_copy)
    _initialize_heat_rate_buffers!(p_heat_copy)
    expected_heat_rates = copy(
        SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
            p_heat_copy,
            u_heat_copy.sc[1],
            1,
            0.0;
            use_buffered_density=false,
        )
    )
    spacecraft_dynamics!(du_heat_copy, u_heat_copy, p_heat_copy, 0.0)
    @test all(isfinite, expected_heat_rates)
    @test any(>(0.0), expected_heat_rates)
    @test du_heat_copy.sc[1].heat_loads == expected_heat_rates

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
    p = ODEParams(n_sats=1, args=args)
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

@testset "Gravity-Gradient Torque Couples Into Rotational RHS" begin
    θ = deg2rad(45.0)
    q0 = normalize(SVector{4, Float64}(0.0, 0.0, sin(θ / 2), cos(θ / 2)))
    ω0 = SVector{3, Float64}(0.0, 0.0, 0.0)
    pos_ii = SVector{3, Float64}(EARTH.Rp_e + 500e3, 0.0, 0.0)
    inertia_tensor = SMatrix{3, 3, Float64}(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0)

    r_hat_body = normalize(SimulationModel.rot(q0) * pos_ii)
    τ_expected = 3.0 * EARTH.μ / norm(pos_ii)^3 * cross(r_hat_body, inertia_tensor * r_hat_body)
    ωdot_expected = inertia_tensor \ τ_expected
    @test norm(τ_expected) > 0.0

    function gravity_gradient_ωdot(effector)
        sc = make_spacecraft(
            ra_alt_m=500e3,
            rp_alt_m=500e3,
            i_deg=35.0,
            ω_deg=40.0,
            Ω_deg=10.0,
            ν_deg=175.0,
            orientation_state=(q0, ω0)
        )
        sc.inertia_tensor = inertia_tensor

        args = build_config(
            spacecraft=sc,
            density_model=NoAtmosphereModel(),
            orientation_sim=true,
            mission_time=10.0,
            EI_km=120.0,
            dynamic_effectors=(effector,),
            keplerian=true
        )

        u0 = build_initial_conditions(args)
        u0.sc[1].pos .= pos_ii
        u0.sc[1].vel .= 0.0
        u0.sc[1].q .= q0
        u0.sc[1].ω .= ω0

        du0 = copy(u0)
        du0 .= 0.0
        p = ODEParams(n_sats=1, args=args)
        spacecraft_dynamics!(du0, u0, p, 0.0)
        return SVector{3, Float64}(du0.sc[1].ω)
    end

    @test isapprox(
        gravity_gradient_ωdot(InverseSquaredGravityModel(gravity_gradient=false)),
        SVector{3, Float64}(0.0, 0.0, 0.0);
        atol=1e-12,
        rtol=1e-10
    )

    for effector in (
        ConstantGravityModel(gravity_gradient=true),
        InverseSquaredGravityModel(gravity_gradient=true),
        InverseSquaredJ2GravityModel(gravity_gradient=true),
    )
        @test isapprox(gravity_gradient_ωdot(effector), ωdot_expected; atol=1e-12, rtol=1e-10)
    end
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
    p = ODEParams(n_sats=1, args=args)
    spacecraft_dynamics!(du0, u0, p, 0.0)

    @test norm(SVector{3, Float64}(du0.sc[1].vel)) > 0.0
    @test all(==(0.0), du0.sc[1].heat_loads)
end

@testset "Atmosphere-Implicit Split RHS Recombines To Full RHS" begin
    q0_split = normalize(SVector{4, Float64}(0.15, -0.1, 0.25, 0.95))
    ω0_split = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc_split = make_spacecraft(
        ra_alt_m=220e3,
        rp_alt_m=100e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=0.0,
        orientation_state=(q0_split, ω0_split)
    )
    explicit_torque = SVector{3, Float64}(0.12, -0.08, 0.05)
    args_split = build_config(
        spacecraft=sc_split,
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=true,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            AerodynamicCoefficientfM(),
            ConstantTorqueModel(explicit_torque),
        ),
        keplerian=true
    )
    args_split.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))

    u0_split = build_initial_conditions(args_split)
    du_full = copy(u0_split)
    du_implicit = copy(u0_split)
    du_explicit = copy(u0_split)
    du_sum = copy(u0_split)
    du_full .= 0.0
    du_implicit .= 0.0
    du_explicit .= 0.0
    du_sum .= 0.0
    p_split = ODEParams(n_sats=1, args=args_split)
    _initialize_heat_rate_buffers!(p_split)
    p_split.shared_buffers.in_atmosphere[1] = false
    p_split.shared_buffers.in_atmosphere_sample_t[1] = -1.0
    @test !SimulationEngine._drag_state_buffer_current(p_split, 1, 0.0)
    @test !SimulationEngine._all_active_spacecraft_outside_atmosphere(u0_split.sc, p_split, 0.0)

    spacecraft_dynamics!(du_full, u0_split, p_split, 0.0)
    spacecraft_dynamics_implicit_atmosphere!(du_implicit, u0_split, p_split, 0.0)
    spacecraft_dynamics_explicit_remainder!(du_explicit, u0_split, p_split, 0.0)
    du_sum .= du_implicit .+ du_explicit

    @test isapprox(du_sum.sc[1].pos, du_full.sc[1].pos; atol=1e-12, rtol=1e-10)
    @test isapprox(du_sum.sc[1].vel, du_full.sc[1].vel; atol=1e-12, rtol=1e-10)
    @test isapprox(du_sum.sc[1].mass, du_full.sc[1].mass; atol=1e-12, rtol=1e-10)
    @test isapprox(du_sum.sc[1].q, du_full.sc[1].q; atol=1e-12, rtol=1e-10)
    @test isapprox(du_sum.sc[1].ω, du_full.sc[1].ω; atol=1e-12, rtol=1e-10)
    @test isapprox(du_sum.sc[1].heat_loads, du_full.sc[1].heat_loads; atol=1e-12, rtol=1e-10)

    @test all(==(0.0), du_implicit.sc[1].pos)
    @test du_implicit.sc[1].mass == 0.0
    @test all(==(0.0), du_implicit.sc[1].q)
    @test all(==(0.0), du_implicit.sc[1].heat_loads)
    @test norm(SVector{3, Float64}(du_implicit.sc[1].vel)) > 0.0
    @test norm(SVector{3, Float64}(du_implicit.sc[1].ω)) > 0.0

    @test isapprox(
        SVector{3, Float64}(du_explicit.sc[1].pos),
        SVector{3, Float64}(u0_split.sc[1].vel);
        atol=1e-12,
        rtol=1e-10
    )
    @test du_explicit.sc[1].mass == 0.0
    @test any(>(0.0), du_explicit.sc[1].heat_loads)

    args_split_no_aero = build_config(
        spacecraft=sc_split,
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=true,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            ConstantTorqueModel(explicit_torque),
        ),
        keplerian=true
    )
    args_split_no_aero.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))
    u0_split_no_aero = build_initial_conditions(args_split_no_aero)
    du_explicit_no_aero = copy(u0_split_no_aero)
    du_explicit_no_aero .= 0.0
    p_split_no_aero = ODEParams(n_sats=1, args=args_split_no_aero)
    _initialize_heat_rate_buffers!(p_split_no_aero)
    spacecraft_dynamics_explicit_remainder!(du_explicit_no_aero, u0_split_no_aero, p_split_no_aero, 0.0)
    @test isapprox(du_explicit.sc[1].ω, du_explicit_no_aero.sc[1].ω; atol=1e-12, rtol=1e-10)
end

@testset "Atmosphere-Implicit Fast Path Requires Current Outside Decision" begin
    q0_outside = normalize(SVector{4, Float64}(0.1, 0.2, -0.05, 0.97))
    ω0_outside = SVector{3, Float64}(0.01, 0.0, -0.015)
    sc_outside = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        orientation_state=(q0_outside, ω0_outside)
    )
    args_outside = build_config(
        spacecraft=sc_outside,
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=true,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            AerodynamicCoefficientfM(),
        ),
        keplerian=true
    )
    args_outside.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))

    u0_outside = build_initial_conditions(args_outside)
    du_outside = copy(u0_outside)
    du_outside .= 0.0
    p_outside = ODEParams(n_sats=1, args=args_outside)
    _initialize_heat_rate_buffers!(p_outside)

    p_outside.shared_buffers.in_atmosphere[1] = true
    p_outside.shared_buffers.in_atmosphere_sample_t[1] = -1.0
    @test !SimulationEngine._drag_state_buffer_current(p_outside, 1, 0.0)
    @test SimulationEngine._all_active_spacecraft_outside_atmosphere(u0_outside.sc, p_outside, 0.0)

    spacecraft_dynamics_implicit_atmosphere!(du_outside, u0_outside, p_outside, 0.0)
    @test all(==(0.0), du_outside.sc[1].pos)
    @test all(==(0.0), du_outside.sc[1].vel)
    @test du_outside.sc[1].mass == 0.0
    @test all(==(0.0), du_outside.sc[1].q)
    @test all(==(0.0), du_outside.sc[1].ω)
    @test all(==(0.0), du_outside.sc[1].heat_loads)

    p_outside.shared_buffers.in_atmosphere[1] = false
    p_outside.shared_buffers.in_atmosphere_sample_t[1] = 0.0
    @test SimulationEngine._drag_state_buffer_current(p_outside, 1, 0.0)
    @test SimulationEngine._all_active_spacecraft_outside_atmosphere(u0_outside.sc, p_outside, 0.0)

    p_outside.shared_buffers.in_atmosphere[1] = true
    p_outside.shared_buffers.in_atmosphere_sample_t[1] = 0.0
    @test SimulationEngine._drag_state_buffer_current(p_outside, 1, 0.0)
    @test !SimulationEngine._all_active_spacecraft_outside_atmosphere(u0_outside.sc, p_outside, 0.0)

    sc_mixed = [
        make_spacecraft(ra_alt_m=220e3, rp_alt_m=100e3, ν_deg=0.0),
        make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=0.0),
    ]
    args_mixed = build_config_multi(
        spacecraft=sc_mixed,
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(AerodynamicCoefficientfM(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    args_mixed.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))
    u0_mixed = build_initial_conditions(args_mixed)
    du_mixed = copy(u0_mixed)
    du_mixed .= 0.0
    p_mixed = ODEParams(n_sats=2, args=args_mixed)
    _initialize_heat_rate_buffers!(p_mixed)
    p_mixed.shared_buffers.in_atmosphere[1] = false
    p_mixed.shared_buffers.in_atmosphere[2] = true
    p_mixed.shared_buffers.in_atmosphere_sample_t .= -1.0
    @test !SimulationEngine._all_active_spacecraft_outside_atmosphere(u0_mixed.sc, p_mixed, 0.0)

    withenv(
        "SPACEAGORA_RHS_EXECUTION_MODE" => "flat",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(max(1, Threads.nthreads()))
    ) do
        spacecraft_dynamics_implicit_atmosphere!(du_mixed, u0_mixed, p_mixed, 0.0)
    end
    @test norm(SVector{3, Float64}(du_mixed.sc[1].vel)) > 0.0
    @test all(==(0.0), du_mixed.sc[2].pos)
    @test all(==(0.0), du_mixed.sc[2].vel)
    @test du_mixed.sc[2].mass == 0.0
    @test all(==(0.0), du_mixed.sc[2].heat_loads)
end

@testset "Split IMEX Allows Zero Implicit Partition" begin
    args_zero_implicit = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    u0_zero_implicit = build_initial_conditions(args_zero_implicit)
    p_zero_implicit = ODEParams(n_sats=1, args=args_zero_implicit)
    _initialize_heat_rate_buffers!(p_zero_implicit)
    split_prob_zero_implicit = withenv("SPACEAGORA_SOLVER_MODE" => "split_imex") do
        _build_typed_solver_problem(u0_zero_implicit, (0.0, 30.0), p_zero_implicit, CallbackSet(), _solver_policy_mode())
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "split_imex", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol_zero_implicit, meta_zero_implicit = _solve_with_solver_policy(split_prob_zero_implicit, _active_solver_config(), args_zero_implicit, 1e-8, 1e-8)
        @test string(sol_zero_implicit.retcode) == "Success"
        @test meta_zero_implicit.solver == "KenCarp4(IMEX)"
    end
end

@testset "Heat Load Derivative Uses Thermal Model" begin
    sc = make_single_link_spacecraft(
        ra_alt_m=220e3,
        rp_alt_m=100e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=0.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=true
    )
    args.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))

    u0 = build_initial_conditions(args)
    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams(n_sats=1, args=args)
    _initialize_heat_rate_buffers!(p)
    spacecraft_dynamics!(du0, u0, p, 0.0)
    drag_force, drag_torque = SimulationModel.calcForceTorque(AerodynamicCoefficientfM(), u0.sc[1], p, 1)

    @test p.shared_buffers.densities[1] > 0.0
    @test isfinite(p.shared_buffers.temperatures[1])
    @test all(isfinite, p.shared_buffers.winds[1])
    @test norm(drag_force) > 0.0
    @test all(isfinite, drag_torque)
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

@testset "Two-Body Kepler Invariants Stay Constant" begin
    ra_alt_m = 1_800e3
    rp_alt_m = 400e3
    sc = make_spacecraft(
        ra_alt_m=ra_alt_m,
        rp_alt_m=rp_alt_m,
        i_deg=28.0,
        ω_deg=15.0,
        Ω_deg=25.0,
        ν_deg=40.0
    )
    ra_m = EARTH.Rp_e + ra_alt_m
    rp_m = EARTH.Rp_e + rp_alt_m
    a0 = 0.5 * (ra_m + rp_m)
    period_s = 2pi * sqrt(a0^3 / EARTH.μ)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=2.5 * period_s,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-9,
            abstol_orbit=1e-9,
            dt_max_orbit=10.0
        )
    )

    df = run_case_silent(args)
    eps = specific_energy(df, EARTH.μ)
    hmag = specific_angular_momentum_magnitude(df)
    a_series = -EARTH.μ ./ (2.0 .* eps)

    eps0 = first(eps)
    h0 = first(hmag)
    a_rel_drift = maximum(abs.((a_series .- a_series[1]) ./ a0))
    eps_rel_drift = maximum(abs.((eps .- eps0) ./ eps0))
    h_rel_drift = maximum(abs.((hmag .- h0) ./ h0))

    @test eps_rel_drift < 1e-5
    @test h_rel_drift < 1e-5
    @test a_rel_drift < 1e-5
end

@testset "Gravity Backbone Improves Long-Horizon Two-Body Energy Drift" begin
    ra_alt_m = 1_400e3
    rp_alt_m = 400e3
    sc = make_spacecraft(
        ra_alt_m=ra_alt_m,
        rp_alt_m=rp_alt_m,
        i_deg=28.0,
        ω_deg=15.0,
        Ω_deg=25.0,
        ν_deg=40.0
    )
    a0 = 0.5 * ((EARTH.Rp_e + ra_alt_m) + (EARTH.Rp_e + rp_alt_m))
    period_s = 2pi * sqrt(a0^3 / EARTH.μ)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=15.0 * period_s,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-5,
            abstol_orbit=1e-5,
            dt_max_orbit=120.0
        )
    )

    df_tsit = DataFrame()
    df_backbone = DataFrame()
    tsit_time = @elapsed begin
        df_tsit = withenv("SPACEAGORA_SOLVER_MODE" => "tsit5") do
            run_case_silent(args)
        end
    end
    backbone_time = @elapsed begin
        df_backbone = withenv(
            "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
            "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "20.0"
        ) do
            run_case_silent(args)
        end
    end

    eps_tsit = specific_energy(df_tsit, EARTH.μ)
    eps_backbone = specific_energy(df_backbone, EARTH.μ)
    drift_tsit = maximum(abs.((eps_tsit .- first(eps_tsit)) ./ first(eps_tsit)))
    drift_backbone = maximum(abs.((eps_backbone .- first(eps_backbone)) ./ first(eps_backbone)))

    @info "gravity_backbone_benchmark" case="two_body_long_horizon" tsit5_seconds=tsit_time gravity_backbone_seconds=backbone_time tsit5_energy_drift=drift_tsit gravity_backbone_energy_drift=drift_backbone
    @test drift_backbone < drift_tsit
end

@testset "J2 Secular Rates Match Standard First-Order Drift" begin
    sc = make_spacecraft(
        ra_alt_m=2_000e3,
        rp_alt_m=400e3,
        i_deg=45.0,
        ω_deg=25.0,
        Ω_deg=40.0,
        ν_deg=20.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=86_400.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-9,
            abstol_orbit=1e-9,
            dt_max_orbit=10.0
        )
    )

    df = run_case_silent(args)
    times = Float64.(df.time)
    Ω_series = Vector{Float64}(undef, nrow(df))
    ω_series = Vector{Float64}(undef, nrow(df))

    @inbounds for idx in 1:nrow(df)
        pos = SVector{3, Float64}(
            Float64(df.sc1_pos_1[idx]),
            Float64(df.sc1_pos_2[idx]),
            Float64(df.sc1_pos_3[idx])
        )
        vel = SVector{3, Float64}(
            Float64(df.sc1_vel_1[idx]),
            Float64(df.sc1_vel_2[idx]),
            Float64(df.sc1_vel_3[idx])
        )
        oe = rvtoorbitalelement(pos, vel, EARTH)
        Ω_series[idx] = Float64(oe[4])
        ω_series[idx] = Float64(oe[5])
    end

    oe0 = rvtoorbitalelement(
        SVector{3, Float64}(Float64(df.sc1_pos_1[1]), Float64(df.sc1_pos_2[1]), Float64(df.sc1_pos_3[1])),
        SVector{3, Float64}(Float64(df.sc1_vel_1[1]), Float64(df.sc1_vel_2[1]), Float64(df.sc1_vel_3[1])),
        EARTH
    )
    Ωdot_expected, ωdot_expected = SimulationModel.GravityEffectors.j2_secular_rates(
        Float64(oe0[1]),
        Float64(oe0[2]),
        Float64(oe0[3]),
        EARTH
    )

    Ωdot_measured = linear_regression_slope(times, unwrap_angle_series(Ω_series))
    ωdot_measured = linear_regression_slope(times, unwrap_angle_series(ω_series))

    # Vallado / Montenbruck-Gill first-order secular rates describe the drift trend.
    # Compare slopes with moderate tolerance because the integrated elements are osculating, not mean.
    @test signbit(Ωdot_measured) == signbit(Ωdot_expected)
    @test signbit(ωdot_measured) == signbit(ωdot_expected)
    @test isapprox(Ωdot_measured, Ωdot_expected; rtol=0.10, atol=0.0)
    @test isapprox(ωdot_measured, ωdot_expected; rtol=0.15, atol=0.0)
end

@testset "Gravity Backbone J2 Preserves Secular Drift Trend" begin
    sc = make_spacecraft(
        ra_alt_m=2_000e3,
        rp_alt_m=400e3,
        i_deg=45.0,
        ω_deg=20.0,
        Ω_deg=15.0,
        ν_deg=40.0
    )
    a0 = 0.5 * ((EARTH.Rp_e + 2_000e3) + (EARTH.Rp_e + 400e3))
    period_s = 2pi * sqrt(a0^3 / EARTH.μ)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=6.0 * period_s,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=30.0
        )
    )

    df = withenv(
        "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
        "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "5.0"
    ) do
        run_case_silent(args)
    end

    times = Float64.(df.time)
    Ω_series = Vector{Float64}(undef, nrow(df))
    ω_series = Vector{Float64}(undef, nrow(df))
    @inbounds for idx in 1:nrow(df)
        pos = SVector{3, Float64}(
            Float64(df.sc1_pos_1[idx]),
            Float64(df.sc1_pos_2[idx]),
            Float64(df.sc1_pos_3[idx])
        )
        vel = SVector{3, Float64}(
            Float64(df.sc1_vel_1[idx]),
            Float64(df.sc1_vel_2[idx]),
            Float64(df.sc1_vel_3[idx])
        )
        oe = rvtoorbitalelement(pos, vel, EARTH)
        Ω_series[idx] = Float64(oe[4])
        ω_series[idx] = Float64(oe[5])
    end

    oe0 = rvtoorbitalelement(
        SVector{3, Float64}(Float64(df.sc1_pos_1[1]), Float64(df.sc1_pos_2[1]), Float64(df.sc1_pos_3[1])),
        SVector{3, Float64}(Float64(df.sc1_vel_1[1]), Float64(df.sc1_vel_2[1]), Float64(df.sc1_vel_3[1])),
        EARTH
    )
    Ωdot_expected, ωdot_expected = SimulationModel.GravityEffectors.j2_secular_rates(
        Float64(oe0[1]),
        Float64(oe0[2]),
        Float64(oe0[3]),
        EARTH
    )

    Ωdot_measured = linear_regression_slope(times, unwrap_angle_series(Ω_series))
    ωdot_measured = linear_regression_slope(times, unwrap_angle_series(ω_series))
    @test signbit(Ωdot_measured) == signbit(Ωdot_expected)
    @test signbit(ωdot_measured) == signbit(ωdot_expected)
    @test isapprox(Ωdot_measured, Ωdot_expected; rtol=0.20, atol=0.0)
    @test isapprox(ωdot_measured, ωdot_expected; rtol=0.25, atol=0.0)
end

@testset "AGORA Earth Regression (Golden)" begin
    golden_path = joinpath(REPO_ROOT, "test", "golden", "agora_earth_regression.csv")
    if !isfile(golden_path)
        @test_skip "Golden regression fixture is not present in this checkout"
        return
    end
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
