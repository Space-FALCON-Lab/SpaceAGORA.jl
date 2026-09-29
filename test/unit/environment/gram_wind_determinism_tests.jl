module GramWindDeterminismTests
# Reproducibility of GRAM's perturbed winds.
#
# GRAM's perturbed winds are a correlated random walk that advances on every
# native update of an atmosphere instance. Two defects made them depend on
# things other than a run's inputs:
#
# - Native Earth-GRAM reads its MAP wind standard-deviation tables before it
#   first loads them, so the first atmosphere updated in a process started from
#   zero wind perturbations (and a NaN north/south component). The first run in
#   a process therefore differed from every later identical run. The GRAMSuite
#   extension now loads the tables once per process before a model is used,
#   and every run reseeds its models right before the solve, so pre-solve
#   probes and earlier runs on the same instance leave no trace.
# - With several threads, per-stage RHS density queries and the threaded
#   density callback advanced the shared walk in scheduling order. With
#   history-dependent winds the run now samples the atmosphere once per
#   accepted step, serially in spacecraft order
#   (SPACEAGORA_DENSITY_FREEZE_PER_STEP=auto), so the result is independent of
#   the thread count.
#
# The native-free testsets run everywhere. The native ones need libGRAM and its
# SPICE kernels; they skip without them unless
# SPACEAGORA_REQUIRE_NATIVE_GRAM_PROBES=1, which turns a skip into a failure.
using Test
using SpaceAGORA
using StaticArrays
using Serialization

const SM = SpaceAGORA.SimulationModel
const EM = SM.EnvironmentModels
const CB = SM.SimulationCallbacks
const SE = SpaceAGORA.SimulationEngine

const WD_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const WD_GRAM_ROOT = joinpath(WD_REPO, "data", "GRAMSuite.jl")
const WD_NATIVE_ROOT = get(ENV, "GRAM_ROOT", joinpath(WD_GRAM_ROOT, "GRAM Suite 2.0"))
const WD_SPICE_PATH = joinpath(WD_NATIVE_ROOT, "SPICE")
const WD_GRAM_LIB = joinpath(WD_NATIVE_ROOT, "Build", "lib",
    Sys.isapple() ? "libGRAM.dylib" : (Sys.iswindows() ? "libGRAM.dll" : "libGRAM.so"))
const WD_GRAM_READY = isfile(WD_GRAM_LIB) && isdir(WD_SPICE_PATH)
const WD_CASE_FILE = joinpath(WD_REPO, "test", "helpers", "gram_wind_determinism_case.jl")

include(joinpath(WD_REPO, "test", "helpers", "native_probe_reporting.jl"))

if WD_GRAM_READY
    if Base.find_package("GRAMSuite") === nothing
        pushfirst!(LOAD_PATH, WD_GRAM_ROOT)
    end
    @eval import GRAMSuite
    include(WD_CASE_FILE)
end

# A raw core whose history dependence the test controls. Unknown cores are
# history-dependent by default (density_models.jl), so `false` must be explicit.
struct WindHistoryCore
    history::Bool
end
EM._gram_core_wind_is_history_dependent(core::WindHistoryCore) = core.history

@testset "history dependence and freeze-per-step resolution" begin
    planet = SM.make_no_gram_planet(:earth)
    @test !EM.density_model_history_dependent(SM.NoAtmosphereModel())
    @test !EM.density_model_history_dependent(SM.ExponentialAtmosphereModel(planet))
    @test EM.density_model_history_dependent(EM.GRAMAtmosphereModel(WindHistoryCore(true)))
    @test !EM.density_model_history_dependent(EM.GRAMAtmosphereModel(WindHistoryCore(false)))
    # Nothing to reset without the extension's native implementation.
    @test !EM.reset_density_model_history!(SM.ExponentialAtmosphereModel(planet))

    for (raw, mode) in ((nothing, :auto), ("auto", :auto), ("", :auto), ("1", :on),
                        ("on", :on), ("0", :off), ("off", :off))
        withenv("SPACEAGORA_DENSITY_FREEZE_PER_STEP" => raw) do
            @test CB._density_freeze_per_step_mode() === mode
            # Without a run, auto keeps the historical default: no freeze.
            @test CB._density_freeze_per_step_enabled() == (mode === :on)
            @test CB._density_freeze_per_step_enabled(false) == (mode === :on)
            @test CB._density_freeze_per_step_enabled(true) == (mode !== :off)
            cfg = CB._snapshot_callback_env_config(; history_dependent=true)
            @test cfg.density_history_dependent
            @test cfg.density_freeze_per_step == (mode !== :off)
            cfg = CB._snapshot_callback_env_config()
            @test !cfg.density_history_dependent
            @test cfg.density_freeze_per_step == (mode === :on)
        end
    end
    withenv("SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "maybe") do
        @test_throws ArgumentError CB._density_freeze_per_step_mode()
    end
end

@testset "a surrogate is history-dependent only when it can reach native GRAM" begin
    planet = SM.make_no_gram_planet(:earth)
    history_base = EM.GRAMAtmosphereModel(WindHistoryCore(true))
    # No point-fallback altitude (the Earth, Mars and Venus default): in-grid
    # queries are table lookups, so the run is not history-dependent.
    grid_only = EM.GRAMAtmosphereModelSurrogate(history_base, "", nothing)
    with_fallback = EM.GRAMAtmosphereModelSurrogate(history_base, "", 60.0e3)
    nominal_fallback = EM.GRAMAtmosphereModelSurrogate(
        EM.GRAMAtmosphereModel(WindHistoryCore(false)), "", 60.0e3)
    @test !EM.density_model_history_dependent(grid_only)
    @test EM.density_model_history_dependent(with_fallback)
    @test !EM.density_model_history_dependent(nominal_fallback)

    env(model) = SM.EnvironmentModel(planet=planet, EI=300.0, density_model=model,
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=true, ephemerides_model=SM.SimpleEphemeridesModel())
    withenv("SPACEAGORA_DENSITY_FREEZE_PER_STEP" => nothing) do
        grid_args = (environment_model=env(grid_only),)
        @test !CB._run_density_history_dependent(grid_args)
        cfg = CB._snapshot_callback_env_config(grid_args)
        @test !cfg.density_history_dependent
        @test !cfg.density_freeze_per_step
        fallback_args = (environment_model=env(with_fallback),)
        @test CB._run_density_history_dependent(fallback_args)
        @test CB._snapshot_callback_env_config(fallback_args).density_freeze_per_step
    end
    # The pure-grid surrogate is lock-free and keeps its threaded callback (the
    # decision itself is checked in the next testset).
    @test CB.density_model_threadsafe(grid_only)
end

# Counts resets; stands in for native GRAM in the pre-solve reset's gating.
mutable struct CountingHistoryModel <: SM.AbstractDensityModel
    history::Bool
    resets::Int
end
EM.density_model_history_dependent(m::CountingHistoryModel) = m.history
EM.reset_density_model_history!(m::CountingHistoryModel) = (m.resets += 1; true)

@testset "the pre-solve reset only touches history-dependent runs" begin
    planet = SM.make_no_gram_planet(:earth)
    function reset_params(main, per_sat, pool; wind=true)
        env = SM.EnvironmentModel(planet=planet, EI=300.0, density_model=main,
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=wind, ephemerides_model=SM.SimpleEphemeridesModel())
        sb = (density_models=Any[per_sat...], gram_isolated_pool_models=Any[pool...],
              gram_density_cache=Any[:stale, :stale], vacuum_gram_caches=Any[:stale])
        return (args=(environment_model=env,), shared_buffers=sb)
    end
    fresh(history=true) = CountingHistoryModel(history, 0)

    main, per_sat, pool = fresh(), [fresh(), fresh()], [fresh(), fresh(false)]
    p = reset_params(main, per_sat, pool)
    SE._reset_density_model_histories!(p)
    @test main.resets == 1
    @test all(m -> m.resets == 1, per_sat)
    @test pool[1].resets == 1
    @test pool[2].resets == 0   # its history does not matter
    @test all(isnothing, p.shared_buffers.gram_density_cache)
    @test all(isnothing, p.shared_buffers.vacuum_gram_caches)

    # Winds off: no query history reaches the dynamics, nothing is reseeded.
    main, per_sat, pool = fresh(), [fresh()], [fresh()]
    p = reset_params(main, per_sat, pool; wind=false)
    SE._reset_density_model_histories!(p)
    @test main.resets == 0 && per_sat[1].resets == 0 && pool[1].resets == 0
    @test p.shared_buffers.gram_density_cache == [:stale, :stale]

    # A history-free configured model: the run is not history-dependent.
    main, per_sat = fresh(false), [fresh()]
    SE._reset_density_model_histories!(reset_params(main, per_sat, CountingHistoryModel[]))
    @test per_sat[1].resets == 0
end

@testset "history-dependent runs evaluate the density callback serially" begin
    planet = SM.make_no_gram_planet(:earth)
    function params(density_model; wind=true)
        env = SM.EnvironmentModel(planet=planet, EI=300.0, density_model=density_model,
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=wind, ephemerides_model=SM.SimpleEphemeridesModel())
        return (args=(environment_model=env,),)
    end
    @test CB._run_density_history_dependent(params(EM.GRAMAtmosphereModel(WindHistoryCore(true))).args)
    @test !CB._run_density_history_dependent(params(EM.GRAMAtmosphereModel(WindHistoryCore(true)); wind=false).args)
    @test !CB._run_density_history_dependent(params(EM.GRAMAtmosphereModel(WindHistoryCore(false))).args)
    @test !CB._run_density_history_dependent(params(SM.ExponentialAtmosphereModel(planet)).args)

    # An explicit `on` and a pinned width both lose to the history guard; the
    # isolated pool's lock-free width query is unaffected by it.
    model = EM.GRAMAtmosphereModel(WindHistoryCore(true))
    root = SM.Link(root=true, m=500.0, ref_area=12.0)
    ic = SM.InitialCondition(ra=planet.Rp_e + 550e3, rp=planet.Rp_e + 550e3,
        i=53.0, ω=0.0, Ω=10.0, ν=0.0)
    spacecraft = SM.SpacecraftModel(SM.Joint[], [root], root, true, 500.0,
        0.0, root.inertia, 0, 0, ic, 1)
    config_for(model) = SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime,
            mission_time=1.0, num_steps_to_save=4, data_rate=2.0),
        environment_model=params(model).args.environment_model,
        dynamics_model=SM.DynamicsModel([spacecraft], (SM.InverseSquaredGravityModel(),)),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
    cfg_args = config_for(model)
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on") do
        for history in (true, false)
            sb = (callback_env_config=Ref(CB._snapshot_callback_env_config(; history_dependent=history)),
                  density_callback_width=Ref(4))
            p = (shared_buffers=sb,)
            decision = CB._density_callback_thread_decision(p, cfg_args, 64)
            if history
                @test !decision.use_threads
                @test decision.allotment == 1
                @test !decision.policy_applied
            else
                @test decision.use_threads == true
                @test decision.allotment == 4
            end
        end
    end

    # A pure-grid surrogate over the same history-dependent core: not frozen,
    # and its callback threads on the lock-free route as before the freeze
    # default existed.
    grid_only = EM.GRAMAtmosphereModelSurrogate(model, "", nothing)
    grid_cfg = config_for(grid_only)
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on",
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => nothing) do
        snapshot = CB._snapshot_callback_env_config(grid_cfg)
        @test !snapshot.density_history_dependent
        @test !snapshot.density_freeze_per_step
        sb = (callback_env_config=Ref(snapshot), density_callback_width=Ref(0))
        decision = CB._density_callback_thread_decision((shared_buffers=sb,), grid_cfg, 64)
        @test decision.use_threads == (Threads.nthreads() > 1)
    end
end

# Runs `script` in a fresh child at each thread count and returns the
# deserialized results. A failed child's stderr is printed, not discarded.
function run_children(script::String, thread_counts; env_overrides=Pair{String, Any}[])
    results = Dict{Int, Any}()
    mktempdir() do dir
        script_path = joinpath(dir, "child.jl")
        write(script_path, script)
        for threads in thread_counts
            out = joinpath(dir, "t$(threads).jls")
            log = joinpath(dir, "t$(threads).log")
            cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(WD_REPO) -t $(threads) $(script_path) $(out)`
            env = copy(ENV)
            delete!(env, "JULIA_NUM_THREADS")
            for (k, v) in env_overrides
                v === nothing ? delete!(env, k) : (env[k] = v)
            end
            ok = success(pipeline(setenv(cmd, env); stdout=log, stderr=log))
            @test ok
            if ok
                results[threads] = deserialize(out)
            else
                @error "Child at $(threads) threads failed" output = read(log, String)
            end
        end
    end
    return results
end

const WD_PROBE_FILE = joinpath(WD_REPO, "test", "helpers", "order_probe_density_case.jl")

# The threaded paths, forced on. At 4 threads and a handful of spacecraft the
# automatic policy keeps the density callback, the RHS atmosphere prefill and
# the per-spacecraft RHS loop serial, so a thread-count comparison under it
# would prove nothing about thread order.
const WD_FORCED_THREADING = [
    "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on",
    "SPACEAGORA_RHS_BATCH_PARALLEL" => "on",
    "SPACEAGORA_RHS_BATCH_THREAD_THRESHOLD" => "2",
    "SPACEAGORA_EFFECTOR_PARALLEL" => "on",
]

@testset "history-dependent results do not depend on the thread count (order probe)" begin
    # A pure-Julia history-dependent model (every query shifts later values),
    # so this runs without native GRAM. Under the default `auto` freeze the
    # result must be identical at 1 and 4 threads with the threaded paths
    # forced on. The negative control (freeze off) must differ, and must have
    # queried from several threads, which shows the comparison can fail.
    script = """
    using SpaceAGORA, StaticArrays, Serialization
    const SM = SpaceAGORA.SimulationModel
    const EM = SM.EnvironmentModels
    const CB = SM.SimulationCallbacks
    include($(repr(WD_PROBE_FILE)))
    args = order_probe_args()
    runs = [order_probe_final(args) for _ in 1:2]
    free = withenv(() -> order_probe_final(args), "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "0")
    serialize(ARGS[1], (runs=runs, free=free, threads=Threads.nthreads()))
    """
    results = run_children(script, (1, 4); env_overrides=[WD_FORCED_THREADING...,
        "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => nothing])
    if length(results) == 2
        one, four = results[1], results[4]
        @test one.threads == 1 && four.threads == 4
        for r in (one, four)
            @test r.runs[1].u == r.runs[2].u
        end
        @test one.runs[1].u == four.runs[1].u
        @test one.runs[1].t == four.runs[1].t
        @test length(four.free.threads) > 1
        @test one.free.u != four.free.u
    end
end

function native_testsets()
    # SpaceAGORA's SPICE kernels, which the GRAM ephemeris hook needs, load
    # with a SPICE-backed planet; runs always build one before querying.
    SM.Earth("", WD_SPICE_PATH)
    @testset "first native atmosphere in a process" begin
        # This process has usually updated GRAM atmospheres already (earlier
        # unit files), so the first-atmosphere defect is checked in a fresh
        # child below. Here: two fresh models agree, and nothing is NaN.
        a = wind_case_first_native_winds(wind_case_model())
        b = wind_case_first_native_winds(wind_case_model())
        @test all(isfinite, a)
        @test isequal(a, b)
    end

    @testset "reseeding reproduces a fresh model" begin
        used = wind_case_model()
        path = [(150.0e3 - 50.0k, 0.3 + 1e-3k, -1.2 + 2e-3k, 0.5k) for k in 0:20]
        q(m) = [EM._gram_point_density(m, x..., true) for x in path]
        reference = q(wind_case_model())
        q(used); q(used)
        @test EM.reset_density_model_history!(used)
        @test isequal(q(used), reference)
        # A raw-core wrapper has no recorded seed and is left alone.
        @test !EM.reset_density_model_history!(EM.GRAMAtmosphereModel(WindHistoryCore(true)))
    end

    withenv("SPACEAGORA_GRAM_WIND_MODE" => "perturbed") do
        @testset "the pre-solve reset leaves every native model initialized" begin
            # Reseeding marks native GRAM uninitialized; the next update would
            # re-run the one-time initialization that reaches CSPICE, which is
            # unsafe on several instances at once. The reset therefore warms
            # every reseeded model again. Observable: each model's next queries
            # equal a fresh model's after the warm query (initialization
            # already taken), not a fresh model's own first queries.
            warmed() = (m = wind_case_model(); CB._warm_gram_pool_model!(m, m); m)
            main = wind_case_model()
            per_sat = [warmed() for _ in 1:2]
            pool = [warmed() for _ in 1:3]
            path = [(150.0e3 - 50.0k, 0.3 + 1e-3k, -1.2 + 2e-3k, 0.5k) for k in 0:10]
            q(m) = [EM._gram_point_density(m, x..., true) for x in path]
            foreach(q, (main, per_sat..., pool...))  # advance every walk
            p = (args=wind_case_args(main),
                 shared_buffers=(density_models=Any[per_sat...], gram_isolated_pool_models=pool,
                                 gram_density_cache=Any[nothing], vacuum_gram_caches=Any[nothing]))
            SE._reset_density_model_histories!(p)
            after_warm = q(warmed())
            first_update = q(wind_case_model())
            @test !isequal(after_warm, first_update)
            for m in (main, per_sat..., pool...)
                @test isequal(q(m), after_warm)
            end
        end
    end

    withenv("SPACEAGORA_GRAM_WIND_MODE" => "perturbed",
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => nothing) do
        @testset "identical runs in one process are identical from the first" begin
            args = wind_case_args(wind_case_model())
            @test CB._snapshot_callback_env_config(args).density_freeze_per_step
            runs = [wind_case_final(args) for _ in 1:3]
            @test all(r -> r.u == runs[1].u && r.t == runs[1].t, runs)
            # A shared, already-advanced instance run without isolation: the
            # pre-solve reseed restores the walk.
            model = args.environment_model.density_model
            [EM._gram_point_density(model, 140.0e3, 0.1k, 0.2k, 1.0k, true) for k in 1:25]
            @test wind_case_final(args; isolate_state=false).u == runs[1].u
            @test wind_case_final(args; isolate_state=false).u == runs[1].u
        end
    end

    withenv("SPACEAGORA_GRAM_WIND_MODE" => "nominal",
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => nothing) do
        @testset "nominal winds keep per-stage sampling" begin
            args = wind_case_args(wind_case_model(); mission_time=20.0)
            cfg = CB._snapshot_callback_env_config(args)
            @test !cfg.density_history_dependent
            @test !cfg.density_freeze_per_step
            default = wind_case_final(args)
            explicit_off = withenv(() -> wind_case_final(args),
                "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "0")
            @test default.u == explicit_off.u
            @test default.t == explicit_off.t
        end
    end

    child_prelude = """
    using SpaceAGORA, Serialization
    const SM = SpaceAGORA.SimulationModel
    const EM = SM.EnvironmentModels
    const CB = SM.SimulationCallbacks
    const SE = SpaceAGORA.SimulationEngine
    const WD_SPICE_PATH = $(repr(WD_SPICE_PATH))
    Base.find_package("GRAMSuite") === nothing && pushfirst!(LOAD_PATH, $(repr(WD_GRAM_ROOT)))
    import GRAMSuite
    include($(repr(WD_CASE_FILE)))
    SM.Earth("", WD_SPICE_PATH)
    """

    @testset "perturbed winds do not depend on the thread count" begin
        # Children at 1 and 4 threads. Each also checks the first atmosphere in
        # its own fresh process, where the static-table defect showed.
        #
        # The threaded paths are forced on (WD_FORCED_THREADING). The negative
        # control lives in the order-probe testset above: per-stage perturbed
        # GRAM with several spacecraft does not finish in test time.
        script = child_prelude * """
        first_winds = wind_case_first_native_winds(wind_case_model())
        second_winds = wind_case_first_native_winds(wind_case_model())
        args = wind_case_args(wind_case_model(); n=8, mission_time=40.0)
        runs = [wind_case_final(args) for _ in 1:2]
        serialize(ARGS[1], (first_winds=first_winds, second_winds=second_winds, runs=runs,
                            threads=Threads.nthreads()))
        """
        results = run_children(script, (1, 4); env_overrides=[WD_FORCED_THREADING...,
            "SPACEAGORA_GRAM_WIND_MODE" => "perturbed",
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => nothing])
        if length(results) == 2
            one, four = results[1], results[4]
            @test one.threads == 1 && four.threads == 4
            for r in (one, four)
                @test all(isfinite, r.first_winds)
                @test isequal(r.first_winds, r.second_winds)
                @test r.runs[1].u == r.runs[2].u
            end
            @test isequal(one.first_winds, four.first_winds)
            @test one.runs[1].u == four.runs[1].u
            @test one.runs[1].t == four.runs[1].t
        end
    end

    @testset "reseeded pool models are safe to update concurrently" begin
        # Eight warmed pool-style instances, reset through the pre-solve reset
        # and then updated concurrently, one thread per instance and each
        # under its own lock, as the isolated pool does. Without the re-warm,
        # every instance's first update after the reseed re-enters GRAM's
        # CSPICE initialization concurrently. The child must survive, and the
        # threaded rounds must equal a serial one.
        script = child_prelude * """
        warmed() = (m = wind_case_model(); CB._warm_gram_pool_model!(m, m); m)
        pool = [warmed() for _ in 1:8]
        p = (args=wind_case_args(pool[1]),
             shared_buffers=(density_models=Any[], gram_isolated_pool_models=pool,
                             gram_density_cache=Any[], vacuum_gram_caches=Any[]))
        point(k, j) = (140.0e3 + 1.0e3 * k, 0.1 * k, 0.2 * k + 0.01 * j, 1.0 * j)
        function one_round(threaded::Bool)
            SE._reset_density_model_histories!(p)
            out = Matrix{Any}(undef, 8, 4)
            function body(k)
                m = pool[k]
                for j in 1:4
                    out[k, j] = EM._gram_core_density_state(m.core, point(k, j)..., true,
                        m.instance_lock, 200.0)
                end
            end
            if threaded
                Threads.@threads :static for k in 1:8
                    body(k)
                end
            else
                foreach(body, 1:8)
            end
            return out
        end
        serial = one_round(false)
        threaded = [one_round(true) for _ in 1:10]
        serialize(ARGS[1], (serial=serial, threaded=threaded, threads=Threads.nthreads()))
        """
        results = run_children(script, (4,);
            env_overrides=["SPACEAGORA_GRAM_WIND_MODE" => "perturbed"])
        if haskey(results, 4)
            r = results[4]
            @test r.threads == 4
            @test all(round -> isequal(round, r.serial), r.threaded)
        end
    end
end

if WD_GRAM_READY
    native_testsets()
else
    @testset "native GRAM wind determinism" begin
        required = NativeProbeReporting.native_probe_required()
        required || @info "Skipping native GRAM wind-determinism tests: no libGRAM on this host." lib = WD_GRAM_LIB
        @test !required
    end
end

end # module
