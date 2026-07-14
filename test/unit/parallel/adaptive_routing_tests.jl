@testset "Adaptive Campaign Routing" begin
    # ── Features from an explicit campaign shape ──────────────────────────────
    feats = SimulationCampaigns.campaign_route_features(
        samples=12, n_sats=2, density_family="exponential", mission_time_s=1200.0
    )
    @test feats isa ParallelProfiles.OuterRouteFeatures
    @test feats.category == "montecarlo"
    @test feats.montecarlo_samples == 12
    @test feats.n_sats == 2
    sig = ParallelProfiles.outer_route_signature(feats)
    @test occursin("cat=montecarlo", sig)
    @test occursin("dens=exp", sig)
    @test occursin("mission=short", sig)
    @test_throws ArgumentError SimulationCampaigns.campaign_route_features(samples=-1)
    @test_throws ArgumentError SimulationCampaigns.campaign_route_features(samples=4, n_sats=0)

    # ── Features derived from a SimulationConfiguration ───────────────────────
    planet_ng = SimulationModel.make_no_gram_planet(:earth)

    function _adaptive_test_sat(id::Int, raan_deg::Float64)
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

    adaptive_sats = [_adaptive_test_sat(1, 0.0), _adaptive_test_sat(2, 120.0)]
    cfg = build_config_multi(
        spacecraft=adaptive_sats,
        density_model=SimulationModel.NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=300.0,
        EI_km=120.0,
        dynamic_effectors=(SimulationModel.InverseSquaredJ2GravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        ephemerides_model=SimulationModel.SimpleEphemeridesModel(),
        planet=planet_ng
    )
    cfg_feats = SimulationCampaigns.campaign_route_features(cfg; samples=2, n_sats=1)
    @test cfg_feats.density_family == "none"
    @test cfg_feats.mission_time_s == 300.0
    @test cfg_feats.n_sats == 1
    @test cfg_feats.montecarlo_samples == 2
    @test cfg_feats.dynamic_effector_count == 1
    @test cfg_feats.has_control == false
    @test cfg_feats.orientation_on == false
    full_feats = SimulationCampaigns.campaign_route_features(cfg; samples=2)
    @test full_feats.n_sats == 2

    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => nothing, "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
        # ── threads=:auto Monte Carlo: route selection, env split, feedback ───
        state = ParallelProfiles.OuterRouteState()
        n_seeds = 6
        inner_budget_seen = Vector{String}(undef, n_seeds)
        outer_active_seen = Vector{String}(undef, n_seeds)
        runner = seed -> begin
            inner_budget_seen[seed] = get(ENV, "SPACEAGORA_INNER_THREAD_BUDGET", "")
            outer_active_seen[seed] = get(ENV, "SPACEAGORA_OUTER_PARALLEL_ACTIVE", "")
            seed * 10
        end
        res1 = SimulationCampaigns.run_monte_carlo(runner, 1:n_seeds; threads=:auto, route_state=state)
        @test res1 isa SimulationCampaigns.MonteCarloResult
        @test length(res1.successful) == n_seeds
        @test [s.value for s in res1.samples] == [10, 20, 30, 40, 50, 60]
        # The adaptive env overrides are scoped to the campaign run.
        @test isempty(get(ENV, "SPACEAGORA_INNER_THREAD_BUDGET", ""))

        expected_workers = Threads.nthreads() > 1 ? min(n_seeds, Threads.nthreads()) : 1
        @test res1.threads == expected_workers
        if Threads.nthreads() > 1
            expected_budget = string(max(1, fld(Threads.nthreads(), expected_workers)))
            @test all(==(expected_budget), inner_budget_seen)
            @test all(==("1"), outer_active_seen)
        end

        auto_sig = ParallelProfiles.outer_route_signature(
            SimulationCampaigns.campaign_route_features(samples=n_seeds)
        )
        route1 = Threads.nthreads() > 1 ? :threads : :none
        snap1 = ParallelProfiles.outer_route_stats_snapshot(state, auto_sig)
        @test haskey(snap1, route1)
        @test snap1[route1].samples == n_seeds
        @test snap1[route1].success_rate == 1.0
        # Feedback stores amortized campaign wall time per sample.
        @test isapprox(snap1[route1].mean_s, res1.elapsed_s / n_seeds; rtol=1e-6)

        # A second identical campaign explores the under-sampled serial route,
        # so repeated campaigns accumulate stats for every feasible allocation.
        if Threads.nthreads() > 1
            res2 = SimulationCampaigns.run_monte_carlo(runner, 1:n_seeds; threads=:auto, route_state=state)
            @test res2.threads == 1
            snap2 = ParallelProfiles.outer_route_stats_snapshot(state, auto_sig)
            @test snap2[:none].samples == n_seeds
            @test snap2[:threads].samples == n_seeds
        end

        # Caller-provided route features are patched with the actual sample count.
        if Threads.nthreads() > 1
            override_state = ParallelProfiles.OuterRouteState()
            stale = SimulationCampaigns.campaign_route_features(samples=0, density_family="exponential")
            res_o = SimulationCampaigns.run_monte_carlo(
                x -> x, 1:n_seeds;
                threads=:auto, route_features=stale, route_state=override_state
            )
            @test res_o.threads == expected_workers
            snap_o = ParallelProfiles.outer_route_stats_snapshot(
                override_state,
                ParallelProfiles.outer_route_signature(
                    SimulationCampaigns.campaign_route_features(samples=n_seeds, density_family="exponential")
                )
            )
            @test sum(info.samples for info in values(snap_o)) == n_seeds
        end

        # Failures are recorded and lower the stored success rate.
        fail_state = ParallelProfiles.OuterRouteState()
        res_fail = SimulationCampaigns.run_monte_carlo(1:4; threads=:auto, route_state=fail_state) do seed
            seed == 2 && error("seed two failed")
            seed
        end
        @test length(res_fail.failed) == 1
        snap_fail = ParallelProfiles.outer_route_stats_snapshot(
            fail_state,
            ParallelProfiles.outer_route_signature(SimulationCampaigns.campaign_route_features(samples=4))
        )
        @test sum(info.samples for info in values(snap_fail)) == 4
        @test any(info.success_rate == 0.75 for info in values(snap_fail))

        # Nested adaptive campaigns yield to an enclosing outer split: serial
        # execution, no route feedback (contended timings must not poison the
        # shared statistics).
        nested_state = ParallelProfiles.OuterRouteState()
        res_nested = withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1") do
            SimulationCampaigns.run_monte_carlo(x -> x, 1:n_seeds; threads=:auto, route_state=nested_state)
        end
        @test res_nested.threads == 1
        @test length(res_nested.successful) == n_seeds
        @test isempty(nested_state.history)

        # Caller-provided features with a non-montecarlo category are normalized
        # before routing and recording: other categories' default-route rules can
        # answer :process for shapes whose candidates are [:none, :threads],
        # which the selector would clamp to a serial default on larger machines.
        if Threads.nthreads() > 1
            cat_state = ParallelProfiles.OuterRouteState()
            odd_features = ParallelProfiles.OuterRouteFeatures(
                category="scaling",
                n_sats=2,
                mission_time_s=3600.0,
                has_control=true,
                montecarlo_samples=0
            )
            res_cat = SimulationCampaigns.run_monte_carlo(
                x -> x, 1:n_seeds;
                threads=:auto, route_features=odd_features, route_state=cat_state
            )
            @test res_cat.threads == expected_workers
            normalized_sig = ParallelProfiles.outer_route_signature(
                ParallelProfiles.OuterRouteFeatures(
                    category="montecarlo",
                    n_sats=2,
                    mission_time_s=3600.0,
                    has_control=true,
                    montecarlo_samples=n_seeds
                )
            )
            snap_cat = ParallelProfiles.outer_route_stats_snapshot(cat_state, normalized_sig)
            @test sum(info.samples for info in values(snap_cat)) == n_seeds
        end

        # ── threads=:auto constellation ensemble ──────────────────────────────
        ens_state = ParallelProfiles.OuterRouteState()
        res_ens = SimulationCampaigns.run_constellation_ensemble(
            cfg; threads=:auto, route_state=ens_state, return_solution=true
        )
        @test res_ens isa SimulationCampaigns.MonteCarloResult
        @test length(res_ens.samples) == 2
        @test isempty(res_ens.failed)
        @test all(s -> s.value !== nothing, res_ens.samples)
        ens_sig = ParallelProfiles.outer_route_signature(
            SimulationCampaigns.campaign_route_features(cfg; samples=2, n_sats=1)
        )
        ens_snap = ParallelProfiles.outer_route_stats_snapshot(ens_state, ens_sig)
        @test sum(info.samples for info in values(ens_snap)) == 2
        @test all(info.success_rate == 1.0 for info in values(ens_snap))
    end

    # ── Argument validation ────────────────────────────────────────────────────
    guard_state = ParallelProfiles.OuterRouteState()
    @test_throws ArgumentError SimulationCampaigns.run_monte_carlo(identity, 1:2; threads=:fast)
    @test_throws ArgumentError SimulationCampaigns.run_monte_carlo(identity, 1:2; threads=1, route_state=guard_state)
    @test_throws ArgumentError SimulationCampaigns.run_constellation_ensemble(cfg; threads=:fast)
    @test_throws ArgumentError SimulationCampaigns.run_constellation_ensemble(cfg; threads=1, route_state=guard_state)

    # Empty seed sets return an empty result without consulting the bandit.
    empty_state = ParallelProfiles.OuterRouteState()
    res_empty = SimulationCampaigns.run_monte_carlo(identity, Int[]; threads=:auto, route_state=empty_state)
    @test isempty(res_empty.samples)
    @test isempty(empty_state.history)

    # The process-global state is a stable, inspectable OuterRouteState.
    @test SimulationCampaigns.campaign_outer_route_state() isa ParallelProfiles.OuterRouteState
    @test SimulationCampaigns.campaign_outer_route_state() === SimulationCampaigns.campaign_outer_route_state()
end
