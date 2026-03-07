@testset "Runtime Boundary Probes" begin
    @testset "SaveData Snapshot Boundary Probe" begin
        args = build_config(
            spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=120.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true,
            simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
        )
        save_fields = SimulationModel.default_save_fields(args)
        u = build_initial_conditions(args)
        p = ODEParams{1}(args=args)
        integrator = (p=p,)

        saved_values_probe = @timed SavedValues(Float64, SimulationModel.SaveData)
        saved_values_allocs = saved_values_probe.gcstats.malloc +
            saved_values_probe.gcstats.realloc +
            saved_values_probe.gcstats.poolalloc +
            saved_values_probe.gcstats.bigalloc

        snapshot_probe = @timed begin
            snapshot = nothing
            for step in 1:8
                snapshot = SimulationModel.SimulationCallbacks._save_snapshot(
                    save_fields,
                    u,
                    Float64(step),
                    integrator
                )
            end
            snapshot
        end
        snapshot_allocs = snapshot_probe.gcstats.malloc +
            snapshot_probe.gcstats.realloc +
            snapshot_probe.gcstats.poolalloc +
            snapshot_probe.gcstats.bigalloc

        @test isfinite(saved_values_probe.time)
        @test saved_values_probe.time >= 0.0
        @test saved_values_probe.bytes >= 0
        @test saved_values_allocs >= 0
        @test saved_values_probe.value isa SavedValues

        @test isfinite(snapshot_probe.time)
        @test snapshot_probe.time >= 0.0
        @test snapshot_probe.bytes >= 0
        @test snapshot_allocs >= 0
        @test snapshot_probe.value isa SimulationModel.SaveData
        @test all(haskey(snapshot_probe.value, key) for key in (:position, :velocity, :mass, :drag, :periapsis_altitude, :heat_rate, :heat_load))
    end

    @testset "Parallel Policy Channel Boundary Probe" begin
        policy = SimulationModel.ParallelPolicy
        budget = max(1, Threads.nthreads())
        items = 32

        withenv(
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
            "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "1"
        ) do
            acc = Base.Threads.Atomic{Int}(0)
            probe = @timed begin
                policy.with_policy_context() do
                    policy.threaded_foreach_persistent(:runtime_boundary_probe, items, budget) do idx
                        Base.Threads.atomic_add!(acc, idx)
                    end
                    policy.threaded_foreach_persistent(:runtime_boundary_probe, items, budget) do idx
                        Base.Threads.atomic_add!(acc, idx)
                    end
                end
                nothing
            end
            allocs = probe.gcstats.malloc + probe.gcstats.realloc + probe.gcstats.poolalloc + probe.gcstats.bigalloc

            @test isfinite(probe.time)
            @test probe.time >= 0.0
            @test probe.bytes >= 0
            @test allocs >= 0
            @test acc[] == 2 * sum(1:items)
        end
    end
end
