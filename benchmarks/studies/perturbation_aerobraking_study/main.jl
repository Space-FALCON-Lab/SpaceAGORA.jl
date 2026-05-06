if !isdefined(@__MODULE__, :REPO_ROOT)
    const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
end
const REPO_PROJECT = joinpath(REPO_ROOT, "Project.toml")

if something(Base.active_project(), "") != REPO_PROJECT
    import Pkg
    Pkg.activate(REPO_ROOT; io=devnull)
end

# Load SpaceAGORA and GRAM. `common.jl` uses !isdefined guards internally,
# so this is safe to call when main.jl is included on Distributed workers.
include(joinpath(REPO_ROOT, "examples", "common.jl"))
setup_gram_example!()

include(joinpath(@__DIR__, "study_config.jl"))
include(joinpath(@__DIR__, "corridor_guidance.jl"))
include(joinpath(@__DIR__, "campaign_setup.jl"))
include(joinpath(@__DIR__, "metrics.jl"))

using CSV
using DataFrames
using Dates
using Distributed

# ---------------------------------------------------------------------------
# Backend and mode flags (read once at load time)
# ---------------------------------------------------------------------------

# "process" (default): pmap over Distributed workers — independent GRAM locks.
# "threads": Threads.@threads with inner parallelism disabled.
# "serial": single-threaded, useful for debugging.
if !isdefined(@__MODULE__, :STUDY_BACKEND)
    const STUDY_BACKEND = get(ENV, "SPACEAGORA_STUDY_BACKEND", "process")
end
if !isdefined(@__MODULE__, :SMOKE_MODE)
    const SMOKE_MODE = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
end

# ---------------------------------------------------------------------------
# Per-worker model caches: load once per worker per body
# ---------------------------------------------------------------------------

if !isdefined(@__MODULE__, :_harmonics_cache)
    const _harmonics_cache = Dict{Symbol,Any}()
end
if !isdefined(@__MODULE__, :_gram_cache)
    const _gram_cache = Dict{Symbol,Any}()
end

function cached_harmonics_model(body_name::Symbol, harmonics_file::String, planet)
    return get!(_harmonics_cache, body_name) do
        GravitationalHarmonicsModel(HARMONICS_DEGREE, HARMONICS_ORDER, harmonics_file, planet)
    end
end

function cached_gram_model(body_name::Symbol, gram_planet_name::String)
    return get!(_gram_cache, body_name) do
        Base.invokelatest(GRAMAtmosphereModel; planet_name=gram_planet_name)
    end
end

# ---------------------------------------------------------------------------
# Single-case runner (executes on workers or the main thread)
# ---------------------------------------------------------------------------

function run_single_case(case::NamedTuple)::NamedTuple
    body_name  = case.body
    ic_params  = case.ic
    pert_level = case.pert
    run_id     = case.run_id
    cfg        = BODY_CONFIGS[body_name]
    out_dir    = case.out_dir

    t_start = time()
    try
        planet    = _make_planet(body_name)
        harmonics = cached_harmonics_model(body_name, cfg.harmonics_file, planet)
        gram      = cached_gram_model(body_name, cfg.gram_planet_name)
        args      = make_perturbation_study_config(
            cfg, ic_params, pert_level, planet, harmonics, gram, run_id;
            results_directory=out_dir,
            smoke_mode=SMOKE_MODE,
        )
        sol       = run_simulation(args; return_solution=true)
        runtime_s = time() - t_start
        retcode   = string(sol.retcode)
        success   = retcode in ("Success", "Terminated", "MaxIters")
        return extract_run_metrics(
            sol, args, cfg, ic_params, pert_level, body_name, run_id;
            solver_success=success,
            solve_retcode=retcode,
            runtime_s=runtime_s,
        )
    catch err
        runtime_s = time() - t_start
        @warn "run_id=$run_id ($body_name/$pert_level) failed" exception=(err, catch_backtrace())
        return failed_run_metrics(
            cfg, ic_params, pert_level, body_name, run_id;
            solve_retcode=string(typeof(err)),
            runtime_s=runtime_s,
        )
    end
end

# ---------------------------------------------------------------------------
# Worker initialisation (process backend)
# ---------------------------------------------------------------------------

if !isdefined(@__MODULE__, :_study_workers_initialized)
    const _study_workers_initialized = Ref(false)
end

function ensure_study_workers!()
    _study_workers_initialized[] && return nothing
    ENV["OPENBLAS_NUM_THREADS"] = "1"
    missing_workers = Sys.CPU_THREADS - nworkers()
    if missing_workers > 0
        addprocs(missing_workers; exeflags=`--startup-file=no --project=$(REPO_ROOT)`)
    end
    study_script = abspath(@__FILE__)
    # Workers include main.jl: activates project, loads SpaceAGORA+GRAM, defines all functions.
    @everywhere workers() include($study_script)
    # Warm up GRAM bindings to avoid Julia 1.12 cross-process deserialization issues.
    @everywhere workers() Base.invokelatest(GRAMAtmosphereModel; planet_name="mars")
    _study_workers_initialized[] = true
    return nothing
end

# ---------------------------------------------------------------------------
# Main study loop
# ---------------------------------------------------------------------------

function run_study(; n_grid::Int=N_GRID, bodies=keys(BODY_CONFIGS))
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    out_dir   = joinpath(OUTPUT_DIR, timestamp)
    mkpath(out_dir)
    println("[study] output → $out_dir")
    eff_grid = SMOKE_MODE ? 1 : n_grid
    println("[study] backend=$STUDY_BACKEND | smoke=$SMOKE_MODE | grid=$(eff_grid)³")

    all_cases = NamedTuple[]
    run_id = 0
    active_bodies = SMOKE_MODE ? [:mars] : sort(collect(bodies))
    for body_name in active_bodies
        cfg      = BODY_CONFIGS[body_name]
        n_ic     = SMOKE_MODE ? 1 : n_grid
        ics      = build_ic_grid(cfg; n=n_ic)
        perts    = SMOKE_MODE ? PERTURBATION_LEVELS[1:1] : PERTURBATION_LEVELS
        for ic in ics, pert in perts
            run_id += 1
            push!(all_cases, (; run_id, body=body_name, ic, pert, out_dir))
        end
    end
    n_total = length(all_cases)
    println("[study] total runs = $n_total")

    rows = Vector{NamedTuple}(undef, n_total)

    if STUDY_BACKEND == "process"
        ensure_study_workers!()
        results = pmap(run_single_case, all_cases)
        copyto!(rows, results)

    elseif STUDY_BACKEND == "threads"
        withenv(
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE"     => "0",
            "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "0",
            "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "0",
        ) do
            Threads.@threads for i in eachindex(all_cases)
                rows[i] = run_single_case(all_cases[i])
            end
        end

    else  # serial
        for i in eachindex(all_cases)
            rows[i] = run_single_case(all_cases[i])
            c = all_cases[i]
            println("  [$(i)/$(n_total)] $(c.body)/$(c.pert) run=$(c.run_id) → $(rows[i].solver_success)")
        end
    end

    df = DataFrame(rows)
    out_csv = joinpath(out_dir, "perturbation_study_results.csv")
    CSV.write(out_csv, df)
    println("[study] wrote $(nrow(df)) rows → $out_csv")
    return df
end

# ---------------------------------------------------------------------------
# Entry point: only executes when this file is the main script, not on workers.
# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    SMOKE_MODE && println("[study] smoke mode active — 1 pert level, short mission times")
    run_study()
end
