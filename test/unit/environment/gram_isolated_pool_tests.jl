module GramIsolatedPoolTests
# The isolated per-worker GRAM pool: what decides that it runs, what it builds,
# and -- when a native GRAM is present on the host -- that it returns the same
# density as the single locked instance, to the bit.
#
# Native GRAM is not thread-safe, so every density query in a process funnels
# through one lock. `SPACEAGORA_GRAM_ISOLATED_POOL` replaces the one shared
# model with `workers` independent `deepcopy`ed models, each behind its own
# lock. That is only admissible if a second instance never disagrees with the
# first, which the native-backed testset below asserts on real GRAM values.
#
# The native-backed part is skipped when GRAM is not built on the host; the rest
# runs everywhere.
using Test
using SpaceAGORA
using StaticArrays

const SM = SpaceAGORA.SimulationModel
const EM = SM.EnvironmentModels
const CB = SM.SimulationCallbacks
const PP = SM.ParallelPolicy

const POOL_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const POOL_GRAM_ROOT = joinpath(POOL_REPO, "data", "GRAMSuite.jl")
const POOL_NATIVE_ROOT = get(ENV, "GRAM_ROOT", joinpath(POOL_GRAM_ROOT, "GRAM Suite 2.0"))
const POOL_SPICE_PATH = joinpath(POOL_NATIVE_ROOT, "SPICE")
const POOL_GRAM_LIB = joinpath(POOL_NATIVE_ROOT, "Build", "lib",
    Sys.isapple() ? "libGRAM.dylib" : (Sys.iswindows() ? "libGRAM.dll" : "libGRAM.so"))
const POOL_GRAM_READY = isfile(POOL_GRAM_LIB) && isdir(POOL_SPICE_PATH)

# Load the extension before the native-free cases create persistent workers.
# Existing worker tasks cannot see methods introduced at a later world age.
if POOL_GRAM_READY
    if Base.find_package("GRAMSuite") === nothing
        pushfirst!(LOAD_PATH, POOL_GRAM_ROOT)
    end
    @eval import GRAMSuite
end

@testset "isolated pool gate" begin
    # :off and :on are unconditional; the default (:auto) additionally needs
    # more than one thread and at least `threshold` items, which is what keeps
    # a single-threaded run from paying for clones it cannot use.
    withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "off") do
        @test CB._gram_isolated_pool_mode() === :off
        @test !CB._gram_isolated_pool_enabled(1024)
    end
    withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "on") do
        @test CB._gram_isolated_pool_mode() === :on
        @test CB._gram_isolated_pool_enabled(1)
        @test !CB._gram_isolated_pool_enabled(0)
    end
    # The shipped defaults, each SOURCED from
    # benchmarks/studies/gram_thread_scaling/results/ and argued in the
    # docstrings in density_callbacks/config.jl. Pinned here because a silent
    # drift in any of the three would either switch the pool on where it was
    # measured to lose or off where it was measured to win.
    withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => nothing,
            "SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD" => nothing,
            "SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS" => nothing) do
        @test CB._gram_isolated_pool_mode() === :auto
        @test CB._gram_isolated_pool_threshold() == 1024
        @test CB._gram_isolated_pool_max_workers() == min(4, max(1, Threads.nthreads()))
        @test CB._gram_isolated_pool_enabled(1024) == (Threads.nthreads() > 1)
        @test !CB._gram_isolated_pool_enabled(1023)
    end
    withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "auto",
            "SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD" => "8") do
        @test CB._gram_isolated_pool_threshold() == 8
        @test CB._gram_isolated_pool_enabled(8) == (Threads.nthreads() > 1)
        @test !CB._gram_isolated_pool_enabled(7)
    end
end

# A history-bearing raw core exercises the real pool dispatcher without native
# GRAM. Pool clones own distinct counters, exactly the distinction under test.
mutable struct PoolWindCore
    history_dependent::Bool
    calls::Int
end
EM._gram_core_wind_is_history_dependent(core::PoolWindCore) = core.history_dependent
function EM._gram_core_density_state(core::PoolWindCore, h::Float64, lat::Float64,
        lon::Float64, t::Float64, wind::Bool, lk, temperature::Float64)
    lock(lk) do
        core.calls += 1
        w = wind && core.history_dependent ? Float64(core.calls) : 7.0
        return 1.0e-8, 200.0, SVector{3,Float64}(w, 0.0, 0.0)
    end
end
function wind_guard_params()
    return (args=(environment_model=(EI=600.0, planet=(T_ref=200.0,)),
                  mission_configuration=(keplerian=true,)),
            shared_buffers=(gram_isolated_pool_models=EM.GRAMAtmosphereModel[],
                            gram_isolated_pool_locks=ReentrantLock[],
                            callback_env_config=Ref(CB._snapshot_callback_env_config())))
end

@testset "automatic pool preserves requested wind histories" begin
    @test EM._gram_core_wind_is_history_dependent(Ref(:unknown))
    n = 8
    hs, lats, lons = fill(150_000.0, n), zeros(n), zeros(n)
    alloc() = (fill(-1.0, n), fill(-2.0, n), fill(SVector{3,Float64}(-3, -4, -5), n))
    withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "auto",
            "SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD" => "1",
            "SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS" => "2") do
        model = EM.GRAMAtmosphereModel(PoolWindCore(true, 0))
        p = wind_guard_params()
        for ts in (0.0, collect(1.0:n))
            rhos, temps, winds = alloc()
            @test !CB._gram_isolated_pool_batch_eval!(rhos, temps, winds, model,
                hs, lats, lons, ts, true, p; allotment_hint=2)
            @test rhos == fill(-1.0, n)
            @test temps == fill(-2.0, n)
            @test winds == fill(SVector{3,Float64}(-3, -4, -5), n)
            @test model.core.calls == 0
            @test isempty(p.shared_buffers.gram_isolated_pool_models)
            @test isempty(p.shared_buffers.gram_isolated_pool_locks)
        end
        # Nominal wind queries and wind=false keep their previous eligibility.
        for (history, wind) in ((false, true), (true, false))
            model = EM.GRAMAtmosphereModel(PoolWindCore(history, 0))
            p = wind_guard_params()
            rhos, temps, winds = alloc()
            pooled = CB._gram_isolated_pool_batch_eval!(rhos, temps, winds, model,
                hs, lats, lons, 0.0, wind, p; allotment_hint=2)
            @test pooled == (Threads.nthreads() > 1)
            if pooled
                @test rhos == fill(1.0e-8, n)
                @test temps == fill(200.0, n)
                @test winds == fill(SVector{3,Float64}(7, 0, 0), n)
                @test length(p.shared_buffers.gram_isolated_pool_models) == 2
                @test model.core.calls == 0
            end
        end
    end
    # Explicit on remains an opt-in to separate stochastic histories.
    withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "on",
            "SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS" => "2") do
        model = EM.GRAMAtmosphereModel(PoolWindCore(true, 0))
        p = wind_guard_params()
        rhos, temps, winds = alloc()
        @test CB._gram_isolated_pool_batch_eval!(rhos, temps, winds, model,
            hs, lats, lons, 0.0, true, p; allotment_hint=2) == (Threads.nthreads() > 1)
        @test model.core.calls == 0
    end
end

@testset "isolated pool width is gated by the density callback's minimum budget" begin
    # A regression guard on an interlock that is easy to trip over and silent
    # when tripped. `_gram_isolated_pool_batch_eval!` takes its width from
    # `_density_callback_thread_decision(...; heavy_work=true).allotment`, and
    # `auto_thread_min_budget(:density_callback)` pins that to 1 on any process
    # with fewer than 16 threads (the floor exists because native GRAM is
    # serialized behind the process-wide lock). The pooled call then fails its
    # own `workers > 1` guard and the run takes the locked path with no
    # diagnostic, however `SPACEAGORA_GRAM_ISOLATED_POOL` is set.
    #
    # Asserted on the policy function rather than through a solve so it holds on
    # a host with no GRAM: it is a property of the width decision, not of GRAM.
    #
    # The pool now asks for its width with `lock_free=true`, which routes it to
    # the lock-free source and past this floor -- see the comment above
    # `_density_callback_thread_decision` in density_callbacks/config.jl. The
    # floor itself is unchanged and still governs the locked path, so it is still
    # pinned here.
    @test PP.auto_thread_min_budget(:density_callback) == 16
    withenv("SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET" => "2") do
        @test PP.auto_thread_min_budget(:density_callback) == 2
    end
    # The lock-free family is deliberately not held to the same floor.
    @test PP.auto_thread_min_budget(:density_callback_lockfree) ==
          PP.auto_thread_min_budget()
end

if !POOL_GRAM_READY
    @info "Skipping native GRAM isolated-pool tests: no libGRAM on this host." lib = POOL_GRAM_LIB
else
    const POOL_MODEL = EM.GRAMAtmosphereModel(planet_name="earth")

    function pool_params(n::Int)
        planet = SM.Earth("", POOL_SPICE_PATH)
        sc = SM.SpacecraftModel[]
        for i in 1:n
            root = SM.Link(root=true, m=500.0, ref_area=12.0)
            ic = SM.InitialCondition(
                ra=planet.Rp_e + 420e3, rp=planet.Rp_e + 330e3,
                i=35.0, ω=40.0, Ω=360.0 * (i - 1) / n, ν=120.0
            )
            push!(sc, SM.SpacecraftModel(
                SM.Joint[], SM.Link[root], root, true, root.m, 0.0, root.inertia, 0, 0, ic, i
            ))
        end
        args = SM.SimulationConfiguration(
            simulation_settings=SM.SimulationSettings(
                results=false, verbose=false, generate_plots=false,
                normalize=false, save_csv=false
            ),
            mission_configuration=SM.MissionConfiguration(
                mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
                mission_time=10.0, orientation_sim=false, num_steps_to_save=10,
                data_rate=10.0
            ),
            environment_model=SM.EnvironmentModel(
                planet=planet, EI=600.0, density_model=POOL_MODEL,
                thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
                topography=false, wind=false
            ),
            dynamics_model=SM.DynamicsModel(sc, (SM.InverseSquaredGravityModel(),)),
            guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
            navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
            control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
            initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
            integration_tolerances=SM.IntegrationTolerances(
                reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0
            )
        )
        return args, SM.ODEParams(n_sats=n, args=args)
    end

    @testset "isolated pool build" begin
        _, p = pool_params(4)
        models, locks = CB._ensure_gram_isolated_pool!(p, POOL_MODEL, 3)
        @test length(models) == 3
        @test length(locks) == 3
        # Independent instances behind independent locks is the whole premise:
        # sharing either one back would reintroduce the serialization the pool
        # exists to remove, or the shared native state it exists to avoid.
        @test all(m -> m !== POOL_MODEL, models)
        @test length(unique(objectid.(models))) == 3
        @test length(unique(objectid.(locks))) == 3
        # Same width reuses; a different width rebuilds. Rebuilding every call
        # would pay a native GRAM construction per callback invocation.
        same, _ = CB._ensure_gram_isolated_pool!(p, POOL_MODEL, 3)
        @test all(i -> same[i] === models[i], 1:3)
        wider, _ = CB._ensure_gram_isolated_pool!(p, POOL_MODEL, 5)
        @test length(wider) == 5
        # A non-positive width is a no-op rather than an error: the caller's
        # width comes from a policy decision that can legitimately be zero.
        kept, _ = CB._ensure_gram_isolated_pool!(p, POOL_MODEL, 0)
        @test length(kept) == 5
    end

    @testset "native wind policy and automatic callback fallback" begin
        for mode in (nothing, "auto", "perturbed", "pert", "stochastic",
                     "nominal", "mean", "deterministic", "base", " NOMINAL ")
            withenv("SPACEAGORA_GRAM_WIND_MODE" => mode) do
                @test EM._gram_core_wind_is_history_dependent(POOL_MODEL.core) ==
                      (GRAMSuite._gram_wind_mode() !== :nominal)
            end
        end
        withenv("SPACEAGORA_GRAM_WIND_MODE" => "invalid") do
            @test_throws ArgumentError EM._gram_core_wind_is_history_dependent(POOL_MODEL.core)
        end
        function callback_samples(pool_mode)
            withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => pool_mode,
                    "SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD" => "1",
                    "SPACEAGORA_DENSITY_BATCH_PARALLEL" => "on",
                    "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
                    "SPACEAGORA_GRAM_WIND_MODE" => "perturbed") do
                original, _ = pool_params(4)
                env = SM.EnvironmentModel(planet=original.environment_model.planet,
                    EI=600.0, density_model=deepcopy(POOL_MODEL), wind=true,
                    thermal_model=original.environment_model.thermal_model,
                    ephemerides_model=SM.SimpleEphemeridesModel())
                args = SM.SimConfig._with_configuration(original; environment_model=env)
                p = SM.ODEParams(n_sats=4, args=args)
                p.shared_buffers.et_start[] = SM.ephemerides_time_seconds(args.initial_time, env.ephemerides_model)
                p.shared_buffers.callback_env_config[] = CB._snapshot_callback_env_config()
                u = SpaceAGORA.SimulationEngine.build_initial_conditions(args)
                cb = CB.get_density_callback(4, args.dynamics_model.dynamic_effectors, args)
                samples = []
                for t in (0.0, 1.0, 2.0)
                    cb.affect!((p=p, u=u, t=t))
                    push!(samples, (copy(p.shared_buffers.densities),
                        copy(p.shared_buffers.temperatures), copy(p.shared_buffers.winds)))
                end
                @test isempty(p.shared_buffers.gram_isolated_pool_models)
                @test isempty(p.shared_buffers.gram_isolated_pool_locks)
                return samples
            end
        end
        # Real callback execution must fall through to the same locked instance
        # history, rather than only returning the expected eligibility flag.
        @test isequal(callback_samples("auto"), callback_samples("off"))
    end

    @testset "isolated pool is bit-identical to the locked path with nominal winds" begin
        # Bit identity is a property of deterministic winds only. GRAM's
        # perturbed winds are a random walk over each instance's own call
        # history, so two instances -- or one instance queried twice -- need
        # not agree, and the GRAMSuite revision CI pins defaults to them. The
        # comparison is made with nominal winds, as gram_density_service_probes.jl
        # does for the same reason.
        withenv("SPACEAGORA_GRAM_WIND_MODE" => "nominal") do
            n = 64
            _, p = pool_params(n)
            hs = Vector{Float64}(undef, n)
            lats = Vector{Float64}(undef, n)
            lons = Vector{Float64}(undef, n)
            ts = Vector{Float64}(undef, n)
            for i in 1:n
                x = (i - 1) / (n - 1)
                hs[i] = 150.0e3 + x * 550.0e3
                lats[i] = -0.5pi + pi * mod(x * sqrt(2.0), 1.0)
                lons[i] = -pi + 2pi * mod(x * sqrt(3.0), 1.0)
                ts[i] = x * 100.0
            end

            alloc() = (zeros(Float64, n), zeros(Float64, n),
                       [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n])
            rho_ref, T_ref, w_ref = alloc()
            EM.getDensityBatch!(rho_ref, T_ref, w_ref, POOL_MODEL, hs, lats, lons, ts, true, p)
            @test count(!iszero, rho_ref) == n

            # `==` on Float64 calls -0.0 equal to 0.0 and every NaN unequal to
            # itself, so the comparison is on the bits.
            bits(v) = reinterpret.(UInt64, v)
            wbits(v) = vcat((reinterpret.(UInt64, collect(x)) for x in v)...)

            workers = max(2, min(4, Threads.nthreads()))
            models, locks = CB._ensure_gram_isolated_pool!(p, POOL_MODEL, workers)
            for k in eachindex(models)
                rho_k, T_k, w_k = alloc()
                for i in 1:n
                    rho_k[i], T_k[i], w_k[i] = CB._gram_isolated_pool_density_state(
                        models[k], hs[i], lats[i], lons[i], ts[i], true, p, locks[k]
                    )
                end
                @test bits(rho_k) == bits(rho_ref)
                @test bits(T_k) == bits(T_ref)
                @test wbits(w_k) == wbits(w_ref)
            end

            if Threads.nthreads() > 1
                rho_p, T_p, w_p = alloc()
                pooled = withenv(
                    "SPACEAGORA_GRAM_ISOLATED_POOL" => "on",
                    "SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS" => string(workers),
                    "SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD" => "1"
                ) do
                    CB._gram_isolated_pool_batch_eval!(
                        rho_p, T_p, w_p, POOL_MODEL, hs, lats, lons, ts, true, p;
                        allotment_hint=workers
                    )
                end
                @test pooled
                @test bits(rho_p) == bits(rho_ref)
                @test bits(T_p) == bits(T_ref)
                @test wbits(w_p) == wbits(w_ref)
            end

            # The guard that keeps the pool from building instances nothing will
            # use. Every item above 2000 km is answered as vacuum and never reaches
            # GRAM, so a batch made entirely of those must be declined -- measured,
            # the speculative build costs 1.83x on a 1024-spacecraft run that never
            # enters the atmosphere.
            @test CB._gram_isolated_pool_native_count(hs, p) == n
            vacuum_hs = fill(2_500_000.0, n)
            @test CB._gram_isolated_pool_native_count(vacuum_hs, p) == 0
            mixed_hs = copy(hs)
            mixed_hs[2:end] .= 2_500_000.0
            @test CB._gram_isolated_pool_native_count(mixed_hs, p) == 1
            rho_v, T_v, w_v = alloc()
            @test !withenv(
                "SPACEAGORA_GRAM_ISOLATED_POOL" => "on",
                "SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS" => string(workers),
                "SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD" => "1"
            ) do
                CB._gram_isolated_pool_batch_eval!(
                    rho_v, T_v, w_v, POOL_MODEL, vacuum_hs, lats, lons, ts, true, p;
                    allotment_hint=workers
                )
            end

            # The pooled batch call declines rather than silently running at width
            # one, which is what lets the caller fall through to the locked path.
            rho_d, T_d, w_d = alloc()
            @test !withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "off") do
                CB._gram_isolated_pool_batch_eval!(
                    rho_d, T_d, w_d, POOL_MODEL, hs, lats, lons, ts, true, p;
                    allotment_hint=workers
                )
            end
        end
    end
end

end # module
