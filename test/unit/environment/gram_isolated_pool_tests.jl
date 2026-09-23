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
const POOL_SPICE_PATH = joinpath(POOL_GRAM_ROOT, "GRAM Suite 2.0", "SPICE")
const POOL_GRAM_LIB = joinpath(POOL_GRAM_ROOT, "GRAM Suite 2.0", "Build", "lib", "libGRAM.so")
const POOL_GRAM_READY = isfile(POOL_GRAM_LIB) && isdir(POOL_SPICE_PATH)

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
    # Loaded the way examples/common.jl and the benchmark harness load it: the
    # vendored GRAMSuite is a separate project, not a dependency of the root
    # environment, and the package extension wires itself in on import.
    if Base.find_package("GRAMSuite") === nothing
        pushfirst!(LOAD_PATH, POOL_GRAM_ROOT)
    end
    @eval import GRAMSuite

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

    @testset "isolated pool is bit-identical to the locked path" begin
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

end # module
