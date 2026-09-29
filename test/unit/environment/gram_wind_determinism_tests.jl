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
    cfg_args = SM.SimulationConfiguration(
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

    @testset "perturbed winds do not depend on the thread count" begin
        # Children at 1 and 4 threads. Each also checks the first atmosphere in
        # its own fresh process, where the static-table defect showed.
        script = """
        using SpaceAGORA, Serialization
        const SM = SpaceAGORA.SimulationModel
        const EM = SM.EnvironmentModels
        const WD_SPICE_PATH = $(repr(WD_SPICE_PATH))
        Base.find_package("GRAMSuite") === nothing && pushfirst!(LOAD_PATH, $(repr(WD_GRAM_ROOT)))
        import GRAMSuite
        include($(repr(WD_CASE_FILE)))
        SM.Earth("", WD_SPICE_PATH)
        first_winds = wind_case_first_native_winds(wind_case_model())
        second_winds = wind_case_first_native_winds(wind_case_model())
        args = wind_case_args(wind_case_model(); n=4)
        runs = [wind_case_final(args) for _ in 1:2]
        serialize(ARGS[1], (first_winds=first_winds, second_winds=second_winds, runs=runs,
                            threads=Threads.nthreads()))
        """
        mktempdir() do dir
            script_path = joinpath(dir, "child.jl")
            write(script_path, script)
            results = Dict{Int, Any}()
            for threads in (1, 4)
                out = joinpath(dir, "t$(threads).jls")
                cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(WD_REPO) -t $(threads) $(script_path) $(out)`
                env = copy(ENV)
                env["SPACEAGORA_GRAM_WIND_MODE"] = "perturbed"
                delete!(env, "SPACEAGORA_DENSITY_FREEZE_PER_STEP")
                delete!(env, "JULIA_NUM_THREADS")
                ok = success(pipeline(setenv(cmd, env); stdout=devnull, stderr=devnull))
                @test ok
                ok && (results[threads] = deserialize(out))
            end
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
