@testset "Run Simulation State Isolation" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=45.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    l_pi_before = args.environment_model.planet.L_PI
    @test_nowarn run_simulation(args)
    @test args.environment_model.planet.L_PI == l_pi_before

    args_no_isolation = deepcopy(args)
    @test_nowarn run_simulation(args_no_isolation; isolate_state=false)

    if Threads.nthreads() >= 2
        args_parallel = deepcopy(args)
        t1 = Threads.@spawn run_simulation(args_parallel)
        t2 = Threads.@spawn run_simulation(args_parallel)
        @test_nowarn fetch(t1)
        @test_nowarn fetch(t2)
    end
end

@testset "RHS Plan Calibration" begin
    # ── Env var parsing ──────────────────────────────────────────────────────
    withenv("SPACEAGORA_RHS_CALIBRATE" => "off") do
        @test SimulationEngine._rhs_calibration_mode() == :off
    end
    withenv("SPACEAGORA_RHS_CALIBRATE" => "auto") do
        @test SimulationEngine._rhs_calibration_mode() == :auto
    end
    withenv("SPACEAGORA_RHS_CALIBRATE" => "force") do
        @test SimulationEngine._rhs_calibration_mode() == :force
    end
    withenv("SPACEAGORA_RHS_CALIBRATE" => "invalid_xyz") do
        @test_throws ArgumentError SimulationEngine._rhs_calibration_mode()
    end

    withenv("SPACEAGORA_RHS_CALIBRATE_N_WARMUP" => "3") do
        @test SimulationEngine._rhs_calibrate_n_warmup() == 3
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_N_WARMUP" => "0") do
        @test SimulationEngine._rhs_calibrate_n_warmup() == 1
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_N_TIMED" => "7") do
        @test SimulationEngine._rhs_calibrate_n_timed() == 7
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_N_TIMED" => "0") do
        @test SimulationEngine._rhs_calibrate_n_timed() == 1
    end

    # ── Machine label stability ───────────────────────────────────────────────
    # The label is cached in a module-level Ref after the first call; test idempotency only.
    label1 = SimulationEngine._calib_machine_label()
    label2 = SimulationEngine._calib_machine_label()
    @test label1 == label2
    @test !isempty(label1)

    # ── Plan constructors return expected field layout ────────────────────────
    sat_batch_plan = SimulationEngine._make_calib_satellite_batch_plan()
    @test sat_batch_plan.mode == :satellite_batch
    @test sat_batch_plan.allotment == 1
    @test sat_batch_plan.dominant_axis == :satellite
    @test sat_batch_plan.policy_applied == true

    flat_plan = SimulationEngine._make_calib_flat_plan(4)
    @test flat_plan.mode == :flat_constellation_effector_queue
    @test flat_plan.allotment == 4
    @test flat_plan.dominant_axis == :flat_effector
    @test flat_plan.policy_applied == true

    flat_plan_floor = SimulationEngine._make_calib_flat_plan(0)
    @test flat_plan_floor.allotment == 1

    # ── In-memory store / lookup round-trip ──────────────────────────────────
    test_sig = "v1|machine=test_gate|budget=8|sats=2_4|effs=1|harm=1"
    SimulationEngine._rhs_calib_store!(test_sig, sat_batch_plan, 1.5e6)
    retrieved = SimulationEngine._rhs_calib_lookup(test_sig)
    @test retrieved !== nothing
    @test retrieved.mode == :satellite_batch

    SimulationEngine._rhs_calib_store!(test_sig, flat_plan, 0.9e6)
    retrieved_flat = SimulationEngine._rhs_calib_lookup(test_sig)
    @test retrieved_flat !== nothing
    @test retrieved_flat.mode == :flat_constellation_effector_queue
    @test retrieved_flat.allotment == 4

    # ── TOML persistence round-trip (temp directory) ─────────────────────────
    mktempdir() do tmp
        calib_path = joinpath(tmp, "calib_test.toml")
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => calib_path) do
            # Store and save to disk
            sig_disk = "v1|machine=disktest|budget=4|sats=2_4|effs=1|harm=1"
            SimulationEngine._rhs_calib_store!(sig_disk, flat_plan, 2.0e6)
            SimulationEngine._rhs_calib_save!()
            @test isfile(calib_path)

            parsed = TOML.parsefile(calib_path)
            @test haskey(parsed, "calibrations")
            @test parsed["schema_version"] == 1
            rows = parsed["calibrations"]
            @test any(r -> get(r, "signature", "") == sig_disk, rows)
        end
    end

    # ── Calibration guard: SPACEAGORA_RHS_CALIBRATE=off skips entirely ───────
    harmonics_file = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    harmonics_model = GravitationalHarmonicsModel(4, 4, harmonics_file, EARTH)

    args_calib = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3 + 10e3 * i, rp_alt_m=480e3 + 5e3 * i, ν_deg=120.0 + 5.0 * i)
            for i in 1:4
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(harmonics_model,),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_calib = ODEParams(n_sats=4, args=args_calib)
    _initialize_heat_rate_buffers!(p_calib)
    _initialize_harmonics_workspace_buffers!(p_calib)
    SimulationEngine._initialize_density_model_instances!(p_calib)
    SimulationEngine._initialize_density_cache_buffers!(p_calib)
    SimulationEngine._initialize_gram_isolated_pool_buffers!(p_calib)
    _initialize_aero_workspace_buffers!(p_calib)
    _initialize_nbody_workspace_buffers!(p_calib)
    u_calib = build_initial_conditions(args_calib)

    withenv("SPACEAGORA_RHS_CALIBRATE" => "off") do
        p_calib.shared_buffers.rhs_plan_override[] = nothing
        SimulationEngine._calibrate_rhs_plan_if_needed!(p_calib, u_calib, args_calib)
        @test p_calib.shared_buffers.rhs_plan_override[] === nothing
    end

    # ── Calibration guard: single satellite skips calibration ────────────────
    args_single_sat = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=480e3, ν_deg=120.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(harmonics_model,),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_single = ODEParams(n_sats=1, args=args_single_sat)
    _initialize_heat_rate_buffers!(p_single)
    _initialize_harmonics_workspace_buffers!(p_single)
    SimulationEngine._initialize_density_model_instances!(p_single)
    SimulationEngine._initialize_density_cache_buffers!(p_single)
    SimulationEngine._initialize_gram_isolated_pool_buffers!(p_single)
    _initialize_aero_workspace_buffers!(p_single)
    _initialize_nbody_workspace_buffers!(p_single)
    u_single = build_initial_conditions(args_single_sat)
    withenv("SPACEAGORA_RHS_CALIBRATE" => "auto") do
        p_single.shared_buffers.rhs_plan_override[] = nothing
        SimulationEngine._calibrate_rhs_plan_if_needed!(p_single, u_single, args_single_sat)
        @test p_single.shared_buffers.rhs_plan_override[] === nothing
    end

    # ── Signature stability ───────────────────────────────────────────────────
    sig_a = SimulationEngine._rhs_calib_signature(p_calib, args_calib.dynamics_model.dynamic_effectors)
    sig_b = SimulationEngine._rhs_calib_signature(p_calib, args_calib.dynamics_model.dynamic_effectors)
    @test sig_a == sig_b
    @test startswith(sig_a, "v1|machine=")
    @test occursin("|sats=", sig_a)
    @test occursin("|effs=1|", sig_a)
    @test occursin("|harm=1", sig_a)

    # ── Multi-sat + harmonics calibration (requires worker threads) ───────────
    if Threads.nthreads() > 1
        # Pin the SIMD batch floor so the 4-sat fixture yields >= 2 viable
        # workers and the sweep produces an override regardless of the
        # default SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER (4).
        # Point the calibration path at a throwaway file so the force-mode
        # save cannot persist this suite's synthetic signatures into the
        # machine-default calibration state (which poisons later suite runs).
        withenv(
            "SPACEAGORA_RHS_CALIBRATE" => "force",
            "SPACEAGORA_RHS_CALIBRATE_N_WARMUP" => "1",
            "SPACEAGORA_RHS_CALIBRATE_N_TIMED" => "2",
            # p_calib only has 4 active satellites; the default SIMD floor of 4
            # sats/worker leaves viable_workers == 1, so no flat-plan candidate
            # is generated to compare against satellite_batch. Loosen the floor
            # so the sweep has more than one candidate to choose from.
            "SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER" => "1",
            "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(mktempdir(), "calib_force_gate.toml")
        ) do
            p_calib.shared_buffers.rhs_plan_override[] = nothing
            SimulationEngine._calibrate_rhs_plan_if_needed!(p_calib, u_calib, args_calib)
            override = p_calib.shared_buffers.rhs_plan_override[]
            @test override !== nothing
            @test override.mode ∈ (:satellite_batch, :flat_constellation_effector_queue)
            @test override.policy_applied == true
        end
    end
end

