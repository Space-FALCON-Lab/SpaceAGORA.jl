using Test
using TOML
using StaticArrays
using LinearAlgebra

# Coverage probes for the R6 routing, RHS calibration and machine-topology
# layers. Included in-process from suite 05 (like the other coverage_*_probes
# files) because the calibration probes need the shared fixtures --
# build_config_multi, make_spacecraft, ODEParams, EARTH, the
# _initialize_*_buffers! chain -- that only exist in that scope.
#
# Every calibration block here points SPACEAGORA_RHS_CALIBRATION_PATH at a temp
# file and, on leaving, empties the in-process cache and clears the loaded flag.
# _rhs_calib_save! rewrites the whole store from the cache, so a cache that
# still holds temp-derived entries when something later saves would truncate
# the real store under output/ to those entries.

const _R6_SE = SimulationEngine
const _R6_SC = SpaceAGORA.SimulationCampaigns
const _R6_PP = SpaceAGORA.ParallelProfiles
const _R6_PC = SimulationModel.ParallelCost
const _R6_CB = SimulationModel.SimulationCallbacks
const _R6_ENV = SimulationModel.EnvironmentModels
const _R6_KIN = SimulationModel.Kinematics

struct _R6CovDensityModel <: SimulationModel.AbstractDensityModel end

function _r6_reset_calib_cache!()
    lock(_R6_SE._rhs_calib_lock) do
        empty!(_R6_SE._rhs_calib_cache)
        _R6_SE._rhs_calib_loaded[] = false
        _R6_SE._rhs_calib_loaded_path[] = ""
        empty!(_R6_SE._rhs_calib_solve_start)
        empty!(_R6_SE._rhs_calib_solve_honoured)
    end
    return nothing
end

# Four satellites, two effectors of different kinds so the plan candidates do
# not collapse to one shape. Harmonics is pre-pass handled and the aero effector
# is flat-queue only, which is the mix the sweep exists to rank.
function _r6_fixture(; n_sats::Int = 4, verbose::Bool = false, mission_time::Float64 = 30.0,
                       density_model = ExponentialAtmosphereModel(EARTH))
    harmonics_file = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    args = build_config_multi(
        spacecraft = [
            make_spacecraft(ra_alt_m = 500e3 + 10e3 * i, rp_alt_m = 480e3 + 5e3 * i, ν_deg = 120.0 + 5.0 * i)
            for i in 1:n_sats
        ],
        density_model = density_model,
        orientation_sim = false,
        mission_time = mission_time,
        EI_km = 120.0,
        dynamic_effectors = (GravitationalHarmonicsModel(4, 4, harmonics_file, EARTH),
                             AerodynamicCoefficientConstant()),
        keplerian = true,
        simulation_settings = SimulationSettings(results = false, verbose = verbose,
                                                 generate_plots = false, normalize = false),
    )
    p = ODEParams(n_sats = n_sats, args = args)
    _initialize_heat_rate_buffers!(p)
    _initialize_harmonics_workspace_buffers!(p)
    _R6_SE._initialize_density_model_instances!(p)
    _R6_SE._initialize_density_cache_buffers!(p)
    _R6_SE._initialize_gram_isolated_pool_buffers!(p)
    _initialize_aero_workspace_buffers!(p)
    _initialize_nbody_workspace_buffers!(p)
    u0 = build_initial_conditions(args)
    return args, p, u0
end

@testset "R6 coverage: link kinematics" begin
    model = make_spacecraft(ra_alt_m = 500e3, rp_alt_m = 480e3)
    root = model.root
    panel = Link(root = false, m = 5.0, ref_area = 0.5, r = MVector{3, Float64}(0.0, 1.0, 0.0))
    rot = _R6_KIN.rot
    @test rotate_to_inertial(model, root, 1) ≈ rot(root.q)'
    @test rotate_to_body(root) == I(3)
    rotate_link(panel, SVector{3, Float64}(0.0, 0.0, 1.0), pi / 2)
    @test panel.q ≈ [0.0, 0.0, sin(pi / 4), cos(pi / 4)]
    @test rotate_to_inertial(model, panel, 1) ≈ rot(root.q)' * rot(panel.q)'
    @test rotate_to_body(panel) ≈ rot(panel.q)'
    # A degenerate axis falls back to the body y axis rather than dividing by zero.
    rotate_link(panel, SVector{3, Float64}(0.0, 0.0, 0.0), 0.3)
    @test panel.q ≈ [0.0, sin(0.15), 0.0, cos(0.15)]
    @test all(isfinite, panel.q)
    rotate_link(panel, SMatrix{3, 3, Float64}(I))
    @test abs(panel.q[4]) ≈ 1.0 atol = 1e-12
    q_new = SVector{4, Float64}(0.0, 0.0, 0.6, 0.8)
    rotate_link(panel, q_new)
    @test panel.q ≈ q_new
    @test_throws AssertionError rotate_link(root, q_new)
end

@testset "R6 coverage: GRAM process-batch gates" begin
    GPB = _R6_CB
    withenv("SPACEAGORA_GRAM_PROCESS_POOL_THRESHOLD" => nothing) do
        @test GPB._gram_process_pool_threshold() == 64
    end
    withenv("SPACEAGORA_GRAM_PROCESS_POOL_THRESHOLD" => "abc") do
        @test GPB._gram_process_pool_threshold() == 64
    end
    withenv("SPACEAGORA_GRAM_PROCESS_POOL_THRESHOLD" => "0") do
        @test GPB._gram_process_pool_threshold() == 1
    end
    withenv("SPACEAGORA_GRAM_PROCESS_POOL_THRESHOLD" => "7") do
        @test GPB._gram_process_pool_threshold() == 7
    end
    withenv("SPACEAGORA_GRAM_PROCESS_POOL" => "off", "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
        @test !GPB._gram_process_pool_enabled(1024)
    end
    # Nested inside an outer process split the pool always declines.
    withenv("SPACEAGORA_GRAM_PROCESS_POOL" => "on", "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1") do
        @test !GPB._gram_process_pool_enabled(1024)
    end
    if GPB.ParallelProcess !== nothing
        withenv("SPACEAGORA_GRAM_PROCESS_POOL" => "on", "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
            @test GPB._gram_process_pool_enabled(1)
        end
        withenv("SPACEAGORA_GRAM_PROCESS_POOL" => "auto", "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
                "SPACEAGORA_GRAM_PROCESS_POOL_THRESHOLD" => "8") do
            @test GPB._gram_process_pool_enabled(8)
            @test !GPB._gram_process_pool_enabled(7)
        end
    end
    # A non-GRAM density model is never a service candidate, whatever the mode.
    _, p_gpb, _ = _r6_fixture()
    withenv("SPACEAGORA_GRAM_PROCESS_POOL" => "on", "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
        @test !GPB._rhs_density_service_candidate(p_gpb, 4)
        @test !GPB._rhs_density_service_candidate(p_gpb, 0)
    end
end

@testset "R6 coverage: RHS calibration buckets, env guards and store parsing" begin
    for (n, bucket) in ((1, "1"), (3, "2_4"), (6, "5_8"), (12, "9_16"), (20, "17_32"),
                        (40, "33_64"), (100, "65_128"), (200, "129_256"), (300, "257p"))
        @test _R6_SE._calib_sat_bucket(n) == bucket
    end

    withenv("SPACEAGORA_RHS_CALIBRATE_OVERRIDE_MARGIN" => "x") do
        @test_throws ArgumentError _R6_SE._rhs_calibrate_override_margin()
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_TIE_MARGIN" => "x") do
        @test_throws ArgumentError _R6_SE._rhs_calibrate_tie_margin()
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_MIN_SOLVE_S" => "x") do
        @test_throws ArgumentError _R6_SE._rhs_calibrate_min_solve_seconds()
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_REVERIFY_SHARE" => "x") do
        @test_throws ArgumentError _R6_SE._rhs_calibrate_reverify_share()
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_REVERIFY_SHARE" => "-0.5") do
        @test _R6_SE._rhs_calibrate_reverify_share() == 0.0
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_REVERIFY_SHARE" => "0.25") do
        @test _R6_SE._rhs_calibrate_reverify_share() == 0.25
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_SCORE" => "best") do
        @test _R6_SE._rhs_calibrate_score_statistic() === :min
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_SCORE" => "bogus") do
        @test_throws ArgumentError _R6_SE._rhs_calibrate_score_statistic()
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_HEURISTIC_VOTES" => "abc") do
        @test _R6_SE._rhs_calibrate_heuristic_votes_needed() == 3
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_HEURISTIC_VOTES" => "0") do
        @test _R6_SE._rhs_calibrate_heuristic_votes_needed() == 1
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_SCHEDULERS" => "dynamic, static") do
        @test _R6_SE._rhs_calibrate_schedulers() == Symbol[:dynamic, :static]
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_SCHEDULERS" => "bogus") do
        @test_throws ArgumentError _R6_SE._rhs_calibrate_schedulers()
    end

    # Machine label: the override is honoured once, then cached.
    saved_label = _R6_SE._CALIB_MACHINE_LABEL[]
    try
        _R6_SE._CALIB_MACHINE_LABEL[] = ""
        withenv("SPACEAGORA_PERF_MACHINE_LABEL" => "r6 cov/box") do
            label = _R6_SE._calib_machine_label()
            @test !isempty(label)
            @test !occursin(" ", label) && !occursin("/", label)
            @test _R6_SE._calib_machine_label() == label
        end
    finally
        _R6_SE._CALIB_MACHINE_LABEL[] = saved_label
    end

    # Store parsing: malformed files and malformed rows degrade to "no entry".
    dir = mktempdir()
    plan_for_save = nothing
    try
        bad = joinpath(dir, "bad.toml")
        write(bad, "calibrations = [ { signature = \n")
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => bad) do
            _r6_reset_calib_cache!()
            @test _R6_SE._rhs_calib_lookup("no_such_signature") === nothing
        end
        notvec = joinpath(dir, "notvec.toml")
        write(notvec, "schema_version = 1\ncalibrations = 5\n")
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => notvec) do
            _r6_reset_calib_cache!()
            @test _R6_SE._rhs_calib_lookup("no_such_signature") === nothing
        end
        mixed = joinpath(dir, "mixed.toml")
        # TOML inline tables must stay on one line each.
        write(mixed, join([
            "schema_version = 1",
            "calibrations = [",
            "  7,",
            "  { mode = \"satellite_batch\" },",
            "  { signature = \"r6cov_row\", mode = \"flat_constellation_effector_queue\", allotment = 2, " *
                "scheduler = \"static\", elapsed_mean_ns = 1234.0, solve_ns = 5.0e9, plan_votes = 2 },",
            "  { signature = \"r6cov_heur\", mode = \"heuristic\", allotment = 1, heuristic_votes = 3 },",
            "]",
        ], "\n") * "\n")
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => mixed) do
            _r6_reset_calib_cache!()
            plan = _R6_SE._rhs_calib_lookup("r6cov_row")
            plan_for_save = plan
            @test plan !== nothing && plan !== :heuristic
            @test plan.mode === :flat_constellation_effector_queue
            @test plan.allotment == 2
            @test _R6_SE._rhs_calib_heuristic_votes("r6cov_heur") == 3
            @test _R6_SE._rhs_calib_heuristic_votes("r6cov_row") == 0
            @test _R6_SE._rhs_calib_lookup("no_such_signature") === nothing
            # A re-pin keeps the shape's solve length and restarts the honoured clock.
            _R6_SE._rhs_calib_store!("r6cov_row", plan, 999.0; sweep_ns = 10.0)
            entry = lock(_R6_SE._rhs_calib_lock) do
                copy(_R6_SE._rhs_calib_cache["r6cov_row"])
            end
            @test entry["solve_ns"] == 5.0e9
            @test entry["honoured_ns"] == 0.0
            @test entry["sweep_ns"] == 10.0
        end
        # A store whose parent "directory" is a regular file cannot be written:
        # the save warns and keeps going rather than failing the solve.
        blocker = joinpath(dir, "blocker")
        write(blocker, "not a directory\n")
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(blocker, "store.toml")) do
            _R6_SE._rhs_calib_store!("r6cov_row", plan_for_save, 1.0)
            @test_logs (:warn, r"failed to save") match_mode = :any _R6_SE._rhs_calib_save!()
        end
    finally
        _r6_reset_calib_cache!()
    end
end

@testset "R6 coverage: RHS sweep verbose path, cached reload and plan precompile" begin
    args_v, p_v, u_v = _r6_fixture(verbose = true)
    effs = args_v.dynamics_model.dynamic_effectors
    dir = mktempdir()
    try
        withenv(
            "SPACEAGORA_RHS_CALIBRATE" => "force",
            "SPACEAGORA_RHS_CALIBRATE_N_WARMUP" => "1",
            "SPACEAGORA_RHS_CALIBRATE_N_TIMED" => "2",
            "SPACEAGORA_RHS_CALIBRATE_INTERLEAVE" => "1",
            "SPACEAGORA_RHS_CALIBRATE_SCORE" => "min",
            "SPACEAGORA_RHS_CALIBRATE_OVERRIDE_MARGIN" => "0.0",
            "SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER" => "1",
            "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(dir, "sweep.toml"),
            "SPACEAGORA_PARALLEL_POLICY_V2" => nothing,
        ) do
            _r6_reset_calib_cache!()
            # The sweep itself, verbose: the ladder print-out and the ranking.
            plan, elapsed, verdict, rival_ns, rival_plan = _R6_SE._run_rhs_sweep!(p_v, u_v, effs, true, args_v)
            @test verdict isa Symbol
            @test p_v.shared_buffers.rhs_plan_override[] === nothing
            if verdict === :pinned
                @test plan !== nothing && elapsed > 0.0
            else
                @test plan === nothing
            end
            # Through the gate: whichever verdict the race produces is stored...
            # The gate itself declines on a one-thread budget (nothing to rank),
            # so the store expectations need a second thread; CI and the
            # suite-05 process have four, the default local entrypoint has one.
            if SimulationModel.ParallelPolicy.effective_inner_thread_budget() > 1
                p_v.shared_buffers.rhs_plan_override[] = nothing
                _R6_SE._calibrate_rhs_plan_if_needed!(p_v, u_v, args_v)
                first_override = p_v.shared_buffers.rhs_plan_override[]
                sig = _R6_SE._rhs_calib_signature(p_v, effs, args_v.environment_model.density_model)
                @test _R6_SE._rhs_calib_lookup(sig) !== nothing
                @test isfile(joinpath(dir, "sweep.toml"))
                # A verdict whose solve length was never recorded (solve_ns 0
                # reads as "sweep regime") is re-verified on the next solve: the
                # second call may sweep again or honour the entry, but the store
                # keeps an entry for the shape either way.
                p_v.shared_buffers.rhs_plan_override[] = nothing
                _R6_SE._calibrate_rhs_plan_if_needed!(p_v, u_v, args_v)
                second_override = p_v.shared_buffers.rhs_plan_override[]
                @test second_override === nothing ||
                      second_override.mode in (:satellite_batch, :flat_constellation_effector_queue)
                @test _R6_SE._rhs_calib_lookup(sig) !== nothing
                # Stamp a short measured solve on the entry: in auto mode a short
                # solve honours whatever verdict is cached, with no sweep (force
                # mode, used above to make the sweep run, never honours one).
                lock(_R6_SE._rhs_calib_lock) do
                    _R6_SE._rhs_calib_cache[sig]["solve_ns"] = 0.5e9
                end
                cached = _R6_SE._rhs_calib_lookup(sig)
                p_v.shared_buffers.rhs_plan_override[] = nothing
                withenv("SPACEAGORA_RHS_CALIBRATE" => "auto") do
                    _R6_SE._calibrate_rhs_plan_if_needed!(p_v, u_v, args_v)
                end
                third_override = p_v.shared_buffers.rhs_plan_override[]
                if cached === nothing || cached === :heuristic
                    @test third_override === nothing ||
                          third_override.mode in (:satellite_batch, :flat_constellation_effector_queue)
                else
                    @test third_override !== nothing
                    @test third_override.mode == cached.mode
                    @test third_override.allotment == cached.allotment
                end
            else
                @test_skip "RHS calibration gate needs an inner thread budget above one"
            end
        end
        # Plan precompilation replaces the sweep when requested.
        withenv("SPACEAGORA_RHS_PRECOMPILE_PLANS" => "1", "SPACEAGORA_RHS_CALIBRATE" => "force",
                "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(dir, "precompile.toml")) do
            _r6_reset_calib_cache!()
            p_v.shared_buffers.rhs_plan_override[] = nothing
            _R6_SE._calibrate_rhs_plan_if_needed!(p_v, u_v, args_v)
            @test p_v.shared_buffers.rhs_plan_override[] === nothing
            @test !isfile(joinpath(dir, "precompile.toml"))
        end
    finally
        p_v.shared_buffers.rhs_plan_override[] = nothing
        _r6_reset_calib_cache!()
    end
end

@testset "R6 coverage: in-run width identification trial" begin
    args_t, p_t, u_t = _r6_fixture()
    effs = args_t.dynamics_model.dynamic_effectors
    dir = mktempdir()
    try
        withenv(
            "SPACEAGORA_RHS_IDENTIFY" => "1",
            "SPACEAGORA_RHS_IDENTIFY_ROUNDS" => "2",
            "SPACEAGORA_RHS_IDENTIFY_MIN_RATIO" => "0.001",
            "SPACEAGORA_RHS_IDENTIFY_PERSIST" => "1",
            "SPACEAGORA_RHS_IDENTIFY_TRACE" => "1",
            "SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER" => "1",
            "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(dir, "identify.toml"),
        ) do
            _r6_reset_calib_cache!()
            @test _R6_SE._rhs_identify_enabled()
            @test _R6_SE._rhs_identify_rounds() == 2
            @test _R6_SE._rhs_estimated_evaluations(p_t) > 0.0
            p_t.shared_buffers.rhs_plan_override[] = nothing
            wt = _R6_SE.build_rhs_width_trial(p_t, effs)
            if wt === nothing
                @test_skip "no width ladder on this thread budget"
            else
                @test length(wt.plans) >= 2
                @test length(wt.widths) == length(wt.plans)
                @test !wt.committed
                du = zero(u_t)
                steps = 0
                while !wt.committed && steps < 200
                    _R6_SE.rhs_width_trial_step!(du, u_t, p_t, 0.0, wt, _R6_SE.spacecraft_dynamics!)
                    steps += 1
                end
                @test wt.committed
                @test steps >= length(wt.plans) * (2 + 1)
                @test p_t.shared_buffers.rhs_width_trial[] === nothing
                verdict = _R6_PC.trial_verdict(wt.trial)
                @test verdict.rounds >= 2
                @test 1 <= verdict.arm <= length(wt.plans)
                speedups = _R6_PC.trial_speedups(wt.trial)
                @test length(speedups) == length(wt.plans)
                @test all(s -> s > 0.0, speedups)
                alpha, beta = _R6_SE._rhs_identify_fit(wt, speedups)
                @test alpha >= 0.0 && beta >= 0.0
                # A committed trial stepped again is a plain dispatch.
                _R6_SE.rhs_width_trial_step!(du, u_t, p_t, 0.0, wt, _R6_SE.spacecraft_dynamics!)
                @test all(isfinite, du)
                # The verdict was persisted, one way or the other.
                @test _R6_SE._rhs_calib_lookup(wt.signature) !== nothing
                # A second build finds the persisted verdict and does not trial again.
                p_t.shared_buffers.rhs_plan_override[] = nothing
                @test _R6_SE.build_rhs_width_trial(p_t, effs) === nothing
            end
        end
    finally
        p_t.shared_buffers.rhs_plan_override[] = nothing
        p_t.shared_buffers.rhs_width_trial[] = nothing
        _r6_reset_calib_cache!()
    end
end

@testset "R6 coverage: density callback width calibration" begin
    args_c, p_c, u_c = _r6_fixture(n_sats = 8, mission_time = 3600.0)
    dir = mktempdir()
    try
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_V2" => "1",
            "SPACEAGORA_CALLBACK_WIDTH_CALIBRATE" => "1",
            "SPACEAGORA_CALLBACK_WIDTH_CALIBRATE_MIN_STEPS" => "0",
            "SPACEAGORA_CALLBACK_WIDTH_CALIBRATE_N_TIMED" => "2",
            "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(dir, "cbw.toml"),
        ) do
            @test _R6_SE._callback_width_calibrate_enabled()
            @test _R6_SE._callback_width_n_timed() == 2
            @test _R6_SE._callback_width_min_steps() == 0.0
            @test _R6_SE._callback_width_candidates(1) == [1]
            @test _R6_SE._callback_width_candidates(6) == [1, 2, 4, 6]
            @test _R6_SE._choose_callback_width(Dict(1 => 10.0, 2 => 4.0), 2, 0.1) == 0
            @test _R6_SE._choose_callback_width(Dict(1 => 4.0, 2 => 10.0), 2, 0.1) == 1
            @test _R6_SE._choose_callback_width(Dict(1 => 4.0), 2, 0.1) == 1
            @test _R6_SE._estimated_accepted_steps(p_c) > 0.0
            _R6_SE._initialize_runtime_env_config!(p_c)
            penv = p_c.shared_buffers.policy_env_config[]
            @test penv !== nothing && penv.policy_v2
            _R6_SE._calibrate_density_callback_width!(p_c, u_c, args_c)
            @test p_c.shared_buffers.density_callback_width[] >= 0
        end
        withenv("SPACEAGORA_CALLBACK_WIDTH_CALIBRATE_N_TIMED" => "abc") do
            @test _R6_SE._callback_width_n_timed() == 4
        end
        withenv("SPACEAGORA_CALLBACK_WIDTH_CALIBRATE_MIN_STEPS" => "-3") do
            @test _R6_SE._callback_width_min_steps() == 200.0
        end
    finally
        p_c.shared_buffers.density_callback_width[] = 0
        _R6_SE._initialize_runtime_env_config!(p_c)
        _r6_reset_calib_cache!()
    end
end

@testset "R6 coverage: lock width cap and identification inside a solve" begin
    args_s, _, _ = _r6_fixture(mission_time = 60.0)
    dir = mktempdir()
    try
        withenv(
            "SPACEAGORA_RHS_CALIBRATE" => "off",
            "SPACEAGORA_RHS_LOCK_WIDTH_CAP" => "1",
            "SPACEAGORA_RHS_IDENTIFY" => "1",
            "SPACEAGORA_RHS_IDENTIFY_ROUNDS" => "2",
            "SPACEAGORA_RHS_IDENTIFY_MIN_RATIO" => "0.001",
            "SPACEAGORA_RHS_IDENTIFY_PERSIST" => "0",
            "SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER" => "1",
            "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(dir, "solve.toml"),
        ) do
            @test _R6_SE._rhs_lock_width_cap_enabled()
            @test_nowarn run_simulation(args_s)
            @test !isfile(joinpath(dir, "solve.toml"))
        end
    finally
        _r6_reset_calib_cache!()
    end
end

@testset "R6 coverage: campaign route state path, persistence and density families" begin
    tn = Base.Threads.nthreads()
    withenv("SPACEAGORA_OUTER_ROUTE_STATE_PATH" => nothing, "SPACEAGORA_PARALLEL_PROFILE" => "R6",
            "SPACEAGORA_PERF_MACHINE_LABEL" => "r6cov") do
        path = _R6_SC.campaign_route_state_path()
        @test isabspath(path)
        @test occursin("outer_route_state_", basename(path))
        @test occursin("r6cov", basename(path))
        @test endswith(basename(path), "_t$(tn).toml")
    end
    withenv("SPACEAGORA_OUTER_ROUTE_STATE_PATH" => joinpath("rel", "state.toml")) do
        path = _R6_SC.campaign_route_state_path()
        @test isabspath(path)
        @test endswith(path, joinpath("rel", "state.toml"))
    end

    dir = mktempdir()
    st = _R6_PP.OuterRouteState()
    try
        bad = joinpath(dir, "bad.toml")
        write(bad, "not = [toml\n")
        withenv("SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "1", "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => bad) do
            _R6_SC.reset_campaign_route_state_persistence!()
            @test _R6_SC.ensure_campaign_route_state_loaded!(st) === nothing   # malformed: cold start
            @test _R6_SC.ensure_campaign_route_state_loaded!(st) === nothing   # loaded once per process
        end
        good = joinpath(dir, "good.toml")
        feat = _R6_SC.campaign_route_features(samples = 8, n_sats = 1, density_family = "none", mission_time_s = 60.0)
        _R6_PP.record_outer_route_feedback!(st, feat; route = :threads, successes = 8, failures = 0,
                                            elapsed_success_s = 1.0, discard_cold_observation = false)
        sig = _R6_PP.outer_route_signature(feat)
        @test haskey(_R6_PP.outer_route_stats_snapshot(st, sig), :threads)
        withenv("SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "1", "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => good) do
            @test _R6_SC.save_campaign_route_state(st) === nothing
            @test isfile(good)
            st2 = _R6_PP.OuterRouteState()
            _R6_SC.reset_campaign_route_state_persistence!()
            @test _R6_SC.ensure_campaign_route_state_loaded!(st2) === nothing
            @test haskey(_R6_PP.outer_route_stats_snapshot(st2, sig), :threads)
        end
        withenv("SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0") do
            @test _R6_SC.save_campaign_route_state(st) === nothing
        end
    finally
        _R6_SC.reset_campaign_route_state_persistence!()
    end

    surrogate = _R6_ENV.GRAMAtmosphereModelSurrogate(ExponentialAtmosphereModel(EARTH), "surrogate.toml", nothing)
    @test _R6_SC._campaign_density_family(surrogate) == "gram_surrogate"
    @test _R6_SC._campaign_density_family(ExponentialAtmosphereModel(EARTH)) == "exponential"
    @test _R6_SC._campaign_density_family(NoAtmosphereModel()) == "none"
    @test _R6_SC._campaign_density_family(_R6CovDensityModel()) == "_r6covdensitymodel"
end

@testset "R6 coverage: machine topology direct reads" begin
    MT = _R6_PP
    core_map = MT._proc_cpuinfo_core_map()
    @test core_map isa Dict{Int, Tuple{Int, Int}}
    q = MT.cgroup_cpu_quota()
    @test q == -1.0 || q > 0.0
    lim = MT.cgroup_memory_limit()
    @test lim == -1 || lim > 0
    allowed = MT.allowed_cpus()
    @test allowed isa Set{Int}
    aff = MT._affinity_physical_cores(MT.physical_core_count())
    @test aff == -1 || aff >= 1
    @test MT.available_memory_bytes() > 0
    @test MT.process_rss_bytes() > 0
    @test MT.memory_local_slot_cap(2) >= 0
    @test MT.memory_local_slot_cap(2; resident = 2) >= 0
    @test MT.memory_local_slot_cap(4; extra_per_worker = 1 << 20, resident = 1) >= 0
    @test MT.memory_worker_cap(extra_per_worker = 1 << 20) >= 0
    @test MT.native_gram_worker_extra_bytes(3) > 0
    withenv("SPACEAGORA_GRAM_SAT_MEMORY_MB" => "1.5") do
        @test MT.native_gram_worker_extra_bytes(2) == 2 * round(Int, 1.5 * (1 << 20))
    end
    withenv("SPACEAGORA_MEMORY_BUDGET_GB" => "0.001") do
        @test MT.memory_worker_cap() == 0
        @test MT.memory_worker_cap(resident = 2) >= 0
        @test MT.memory_local_slot_cap(2; resident = 2) >= 0
    end
    withenv("SPACEAGORA_MEMORY_BUDGET_GB" => "1024", "SPACEAGORA_PERF_WORKER_MEMORY_GB" => "0.001") do
        @test MT.memory_worker_cap(resident = 3) >= 3
    end
    @test MT._positive_int_env("SPACEAGORA_R6COV_UNSET_INT") == -1
    @test MT._positive_float_env("SPACEAGORA_R6COV_UNSET_FLOAT") == -1.0
    withenv("SPACEAGORA_R6COV_INT" => "abc") do
        @test_throws ArgumentError MT._positive_int_env("SPACEAGORA_R6COV_INT")
    end
    withenv("SPACEAGORA_R6COV_FLOAT" => "-2.0") do
        @test_throws ArgumentError MT._positive_float_env("SPACEAGORA_R6COV_FLOAT")
    end
    withenv("SPACEAGORA_R6COV_FLOAT" => "2.5") do
        @test MT._positive_float_env("SPACEAGORA_R6COV_FLOAT") == 2.5
    end
end

@testset "R6 coverage: robust timing sink and campaign dispatcher warm-up" begin
    _R6_PC._TIMING_SINK[] = 0.0
    @test _R6_PC._consume_sink(NaN) === nothing
    @test isnan(_R6_PC._TIMING_SINK[])
    _R6_PC._TIMING_SINK[] = 0.0
    @test _R6_PC._consume_sink(1.0) === nothing
    @test _R6_PC._TIMING_SINK[] == 0.0
    @test _R6_PC._sign_test_two_sided(0, 0) == 1.0
    @test _R6_PC._sign_test_two_sided(5, 5) ≈ 2.0 / 32.0
    @test _R6_PC._sign_test_two_sided(3, 5) == 1.0
    # The body of the package's precompile workload, callable at run time.
    @test _R6_SC._warm_campaign_dispatchers() === nothing
end
