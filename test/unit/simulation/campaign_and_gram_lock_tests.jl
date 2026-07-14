@testset "Constellation Ensemble Campaign" begin
    planet_ng = SimulationModel.make_no_gram_planet(:earth)

    function _ensemble_test_sat(id::Int, raan_deg::Float64)
        root = Link{0}(root=true, m=100.0, ref_area=2.0)
        ic = InitialCondition(
            ra=planet_ng.Rp_e + 700e3,
            rp=planet_ng.Rp_e + 650e3,
            i=45.0,
            ω=0.0,
            Ω=raan_deg,
            ν=10.0
        )
        return SpacecraftModel(Joint[], [root], root, true, root.m, 0.0, root.inertia, 0, 0, ic, id)
    end

    sats = [_ensemble_test_sat(11, 0.0), _ensemble_test_sat(22, 40.0), _ensemble_test_sat(33, 80.0)]
    results_dir = mktempdir()
    cfg = build_config_multi(
        spacecraft=sats,
        density_model=SimulationModel.NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=600.0,
        EI_km=120.0,
        dynamic_effectors=(SimulationModel.InverseSquaredJ2GravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(
            results=true,
            verbose=false,
            generate_plots=false,
            normalize=false,
            results_directory=results_dir
        ),
        ephemerides_model=SimulationModel.SimpleEphemeridesModel(),
        planet=planet_ng
    )

    ensemble_threads = min(2, Threads.nthreads())
    res = SimulationCampaigns.run_constellation_ensemble(cfg; threads=ensemble_threads, return_solution=true)
    @test res isa SimulationCampaigns.MonteCarloResult
    @test length(res.samples) == 3
    @test isempty(res.failed)
    @test [s.seed for s in res.samples] == [1, 2, 3]
    @test all(s -> s.value !== nothing, res.samples)

    # Each ensemble member writes into its own per-satellite results directory.
    @test isdir(joinpath(results_dir, "sat_1_id_11"))
    @test isdir(joinpath(results_dir, "sat_2_id_22"))
    @test isdir(joinpath(results_dir, "sat_3_id_33"))

    # An ensemble member must reproduce a direct single-satellite propagation.
    single_cfg = SimulationCampaigns._ensemble_member_configuration(cfg, sats[3], "solo_direct")
    sol_direct = run_simulation(single_cfg; return_solution=true)
    u_direct = collect(sol_direct.u[end])
    u_member = collect(res.samples[3].value.u[end])
    @test maximum(abs.(u_direct .- u_member)) <= 1e-9

    # Uncoupled guard: GNC effectors are rejected unless explicitly allowed.
    cfg_ctrl = build_config_multi(
        spacecraft=sats,
        density_model=SimulationModel.NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=600.0,
        EI_km=120.0,
        dynamic_effectors=(SimulationModel.InverseSquaredJ2GravityModel(),),
        control_effectors=(make_base_thruster_model(thrust=0.0),),
        control_rates=[1.0],
        keplerian=true,
        ephemerides_model=SimulationModel.SimpleEphemeridesModel(),
        planet=planet_ng
    )
    @test_throws ArgumentError SimulationCampaigns.run_constellation_ensemble(cfg_ctrl)
    @test SimulationCampaigns._validate_ensemble_uncoupled(cfg_ctrl, true) === nothing

    # Empty constellation is rejected.
    cfg_empty = build_config_multi(
        spacecraft=SpacecraftModel[],
        density_model=SimulationModel.NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=600.0,
        EI_km=120.0,
        dynamic_effectors=(SimulationModel.InverseSquaredJ2GravityModel(),),
        keplerian=true,
        ephemerides_model=SimulationModel.SimpleEphemeridesModel(),
        planet=planet_ng
    )
    @test_throws ArgumentError SimulationCampaigns.run_constellation_ensemble(cfg_empty)

    # Member splits alias the parent environment (cheap split); parallel isolation
    # therefore relies on the per-worker deepcopy, which must sever mutable model
    # state such as the planet object.
    member_split = SimulationCampaigns._ensemble_member_configuration(cfg, sats[1], "alias_probe")
    @test member_split.environment_model === cfg.environment_model
    member_isolated = deepcopy(member_split)
    @test member_isolated.environment_model !== cfg.environment_model
    @test member_isolated.environment_model.planet !== cfg.environment_model.planet

    # An explicit isolate_state=false must not break the parallel path: the runner
    # overrides it with its own per-worker copy.
    res_noiso = SimulationCampaigns.run_constellation_ensemble(
        cfg; threads=ensemble_threads, return_solution=true, isolate_state=false
    )
    @test isempty(res_noiso.failed)
    u_noiso = collect(res_noiso.samples[3].value.u[end])
    @test maximum(abs.(u_direct .- u_noiso)) <= 1e-9

    # Checkpoint path splitting: any checkpoint interaction gets a per-member path.
    settings_ck_default = SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false,
        results_directory="outdir", checkpoint_enabled=true
    )
    split_ck_default = SimulationCampaigns._ensemble_member_settings(settings_ck_default, "sat_1_id_11")
    @test split_ck_default.results_directory == joinpath("outdir", "sat_1_id_11")

    settings_resume_explicit = SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false,
        results_directory="outdir", checkpoint_directory="ckdir", resume_from_checkpoint=true
    )
    split_resume = SimulationCampaigns._ensemble_member_settings(settings_resume_explicit, "sat_2_id_22")
    @test split_resume.checkpoint_directory == joinpath("ckdir", "sat_2_id_22")
    @test split_resume.results_directory == "outdir"

    settings_plain = SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false
    )
    @test SimulationCampaigns._ensemble_member_settings(settings_plain, "sat_3_id_33") === settings_plain
end

@testset "GRAM Lock Scope" begin
    EMt = SimulationModel.EnvironmentModels

    withenv("SPACEAGORA_GRAM_LOCK_SCOPE" => nothing) do
        @test EMt._gram_lock_scope() === :global
    end
    withenv("SPACEAGORA_GRAM_LOCK_SCOPE" => "global") do
        @test EMt._gram_lock_scope() === :global
    end
    for token in ("model", "per_model", "per-model", "instance", "MODEL")
        withenv("SPACEAGORA_GRAM_LOCK_SCOPE" => token) do
            @test EMt._gram_lock_scope() === :model
        end
    end
    withenv("SPACEAGORA_GRAM_LOCK_SCOPE" => "bogus") do
        @test_throws ArgumentError EMt._gram_lock_scope()
    end

    # Every wrapper construction carries its own instance lock so distinct
    # models never contend under SPACEAGORA_GRAM_LOCK_SCOPE=model.
    m1 = EMt.GRAMAtmosphereModel(:dummy_core_a)
    m2 = EMt.GRAMAtmosphereModel(:dummy_core_b)
    @test m1.instance_lock isa ReentrantLock
    @test m1.instance_lock !== m2.instance_lock
    @test m1.core === :dummy_core_a
    @test :instance_lock in propertynames(m1)
end

