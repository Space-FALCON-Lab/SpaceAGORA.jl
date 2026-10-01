# One run of the shortened Odyssey MarsGRAM arc under one GRAM density mode.
#
#   julia --project=. benchmarks/studies/telemetry_validation/perturbed_density_modes/run_mode.jl \
#       --tag=A_s11 --mode=step --seed=11 --density-scale=1.0 --orbits=N --out=DIR \
#       [--tight=true] [--maxiters=K] [--pass-dt=1.0] [--log=true]
#
# Base case: manifests/odyssey_marsgram.toml (the record configuration),
# truncated to N orbits through the runner's own builders (_make_orbit_args,
# _with_study_settings, the runner's solver environment). Only
# [scenarios.atmosphere_truth] gram_seed / gram_perturbation_scales change:
#   mode=off  -> the record setting itself (seed 1001, scales 0): the nominal,
#                mean-density baseline;
#   otherwise -> scales (density_scale, 0, 0, 0): GRAM's density perturbation at
#                the given scale, wind perturbation scales off (winds are the
#                nominal mean field regardless: SPACEAGORA_GRAM_WIND_MODE=auto).
# --dt-max-atm=S caps the in-atmosphere step at S seconds (study default 0.2 s
# for this scenario) through SPACEAGORA_TELEMETRY_DT_MAX_ATM.
# --tight=true tightens all four integration tolerances 10x through the study's
# SPACEAGORA_TELEMETRY_{RELTOL,ABSTOL}_{ORBIT,ATM} hooks (which may only tighten).
#
# Writes to DIR/<tag>/: run_summary.toml, extrema.csv (per-orbit apsides, record
# extraction), simulation_results.csv, and perturbation_log.csv (modes != off).

include(joinpath(@__DIR__, "..", "common.jl"))
load_gramsuite!()

using SpaceAGORA
using SpaceAGORA.TelemetryVerification
using DataFrames
using CSV
using TOML
using Dates
using Printf

const TV = SpaceAGORA.TelemetryVerification
const SE = SpaceAGORA.SimulationEngine

const OPTS = parse_kv_args(copy(ARGS))
const TAG = OPTS["tag"]
const MODE = get(OPTS, "mode", "off")
const SEED = parse(Int, get(OPTS, "seed", "1001"))
const DSCALE = parse(Float64, get(OPTS, "density-scale", "1.0"))
const N_ORBITS = parse(Int, OPTS["orbits"])
const TIGHT = parse_bool_flag(get(OPTS, "tight", "false"))
const DT_MAX_ATM = get(OPTS, "dt-max-atm", "")
const MAXITERS = haskey(OPTS, "maxiters") ? parse(Int, OPTS["maxiters"]) : nothing
const PASS_DT = get(OPTS, "pass-dt", "1.0")
const WRITE_LOG = parse_bool_flag(get(OPTS, "log", "true"))
const OUT = abspath(joinpath(OPTS["out"], TAG))
mkpath(OUT)

const RECORD_MANIFEST = joinpath(MANIFEST_DIR, "odyssey_marsgram.toml")

function run_cfg()
    doc = TOML.parsefile(RECORD_MANIFEST)
    truth = doc["scenarios"][1]["atmosphere_truth"]
    if MODE != "off"
        truth["gram_seed"] = SEED
        truth["gram_perturbation_scales"] = [DSCALE, 0.0, 0.0, 0.0]
    end
    path = joinpath(mktempdir(), "manifest.toml")
    open(io -> TOML.print(io, doc), path, "w")
    return only(cd(() -> TV._load_scenarios_from_manifest(path), REPO_ROOT)), truth
end

const TIGHT_ENV = TIGHT ? [
    "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "1e-8",
    "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => "1e-10",
    "SPACEAGORA_TELEMETRY_RELTOL_ATM" => "1e-8",
    "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => "1e-10",
] : Pair{String, String}[]
const DT_ENV = isempty(DT_MAX_ATM) ? Pair{String, String}[] : ["SPACEAGORA_TELEMETRY_DT_MAX_ATM" => DT_MAX_ATM]

function main()
    cfg, truth = run_cfg()
    TV._planet_from_name(cfg.planet_name)       # furnish kernels as a solve does
    args = withenv(TIGHT_ENV..., DT_ENV...) do
        TV._with_study_settings(TV._make_orbit_args(cfg, N_ORBITS); quick=false)
    end
    tmp = mktempdir()
    cfg_run = TV._with_configuration(args;
        simulation_settings=SpaceAGORA.SimulationModel.SimulationSettings(
            results=true, verbose=false, results_directory=tmp, generate_plots=false,
            generate_filenames=false, normalize=false, save_csv=true),
        solver_config=nothing)
    solver_mode = TV._telemetry_solver_mode()
    maxiters = MAXITERS === nothing ? TV._telemetry_solver_maxiters(:full) : MAXITERS
    log_path = (MODE != "off" && WRITE_LOG) ? joinpath(OUT, "perturbation_log.csv") : ""
    env = [
        "SPACEAGORA_WARN_NORMALIZE" => "0",
        "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
        "SPACEAGORA_SOLVER_MODE" => solver_mode,
        "SPACEAGORA_SOLVER_MAXITERS" => string(maxiters),
        "SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => TV._telemetry_solver_save_env("SPACEAGORA_SOLVER_SAVE_EVERYSTEP", solver_mode),
        "SPACEAGORA_SOLVER_SAVE_ON" => TV._telemetry_solver_save_env("SPACEAGORA_SOLVER_SAVE_ON", solver_mode),
        "SPACEAGORA_GRAM_OFFLINE_SURROGATE" => cfg.atmosphere_truth.gram_offline_surrogate,
        "SPACEAGORA_GRAM_STATIC_GRID" => cfg.atmosphere_truth.gram_static_grid ? "on" : "off",
        "SPACEAGORA_GRAM_TRACK_CACHE" => cfg.atmosphere_truth.gram_track_cache ? "on" : "off",
        "SPACEAGORA_GRAM_GLOBAL_LOCK" => cfg.atmosphere_truth.gram_global_lock,
        "SPACEAGORA_GRAM_DENSITY_PERTURBATION" => MODE,
        "SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_DT_S" => PASS_DT,
        "SPACEAGORA_GRAM_DENSITY_PERTURBATION_LOG" => log_path,
    ]
    started = now(UTC)
    result = nothing
    err_text = ""
    wall_s = @elapsed begin
        try
            result = withenv(env...) do
                cd(tmp) do
                    SE.run_simulation(cfg_run; isolate_state=false, save_fields=TV._save_fields_for_study(),
                                      return_solution=true, return_solver_metadata=true,
                                      extra_callbacks=TV._scenario_extra_callbacks(cfg))
                end
            end
        catch err
            err_text = sprint(showerror, err)
            @error "run failed" tag=TAG exception=(err, catch_backtrace())
        end
    end
    finished = now(UTC)

    csv_path = joinpath(tmp, "simulation_results.csv")
    n_peri = 0; n_apo = 0; final_t = NaN
    if isfile(csv_path)
        cp(csv_path, joinpath(OUT, "simulation_results.csv"); force=true)
        df = CSV.read(csv_path, DataFrame)
        final_t = Float64(df.time[end])
        ex = TV._extract_extrema_series(df, args.environment_model.planet, cfg.orbit_altitude_mode)
        n_peri = length(ex.peri.altitude); n_apo = length(ex.apo.altitude)
        off = something(cfg.epoch_orbit_offset, 0.0)
        CSV.write(joinpath(OUT, "extrema.csv"), DataFrame(
            event=vcat(fill("peri", n_peri), fill("apo", n_apo)),
            index=vcat(1:n_peri, 1:n_apo),
            flight_orbit=vcat(off .+ (0:n_peri-1), off .+ (0:n_apo-1)),
            time_s=vcat(ex.peri.time_s, ex.apo.time_s),
            altitude_km=vcat(ex.peri.altitude, ex.apo.altitude)))
    end

    stats = result === nothing ? nothing : result.solution.stats
    tol = args.integration_tolerances
    summary = Dict{String, Any}(
        "tag" => TAG, "mode" => MODE,
        "gram_seed" => truth["gram_seed"], "gram_perturbation_scales" => truth["gram_perturbation_scales"],
        "orbits_requested" => N_ORBITS, "tight" => TIGHT,
        "reltol_orbit" => tol.reltol_orbit, "abstol_orbit" => tol.abstol_orbit,
        "reltol_atmosphere" => tol.reltol_atmosphere, "abstol_atmosphere" => tol.abstol_atmosphere,
        "dt_max_orbit" => tol.dt_max_orbit, "dt_max_atmosphere" => tol.dt_max_atmosphere,
        "solver_mode" => solver_mode, "maxiters" => maxiters, "pass_dt_s" => PASS_DT,
        "wall_s" => wall_s, "started_utc" => string(started), "finished_utc" => string(finished),
        "retcode" => result === nothing ? "ERROR" : string(result.solution.retcode),
        "error" => err_text,
        "naccept" => stats === nothing ? -1 : Int(stats.naccept),
        "nreject" => stats === nothing ? -1 : Int(stats.nreject),
        "nf" => stats === nothing ? -1 : Int(stats.nf),
        "solver_sequence" => result === nothing ? "" : join([string(m.solver) for m in result.solver_trace], "->"),
        "solver_fallbacks" => result === nothing ? -1 : count(m -> m.fallback_used, result.solver_trace),
        "final_t_s" => final_t, "n_peri" => n_peri, "n_apo" => n_apo,
        "commit" => get(ENV, "SPACEAGORA_RUN_COMMIT", "unknown"),
        "host" => gethostname(), "julia_threads" => Threads.nthreads(),
    )
    open(io -> TOML.print(io, summary), joinpath(OUT, "run_summary.toml"), "w")
    @printf("[%s] mode=%s seed=%d wall=%.1f s retcode=%s naccept=%d nreject=%d nf=%d peri=%d apo=%d\n",
            TAG, MODE, summary["gram_seed"], wall_s, summary["retcode"], summary["naccept"],
            summary["nreject"], summary["nf"], n_peri, n_apo)
    return nothing
end

main()
