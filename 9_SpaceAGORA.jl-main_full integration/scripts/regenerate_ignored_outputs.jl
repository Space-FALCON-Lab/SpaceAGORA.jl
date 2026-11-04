#!/usr/bin/env julia

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

function _usage()
    println("""
    Usage:
      julia --project=. scripts/regenerate_ignored_outputs.jl runtime-analysis [quick|full] [outdir]
      julia --project=. scripts/regenerate_ignored_outputs.jl telemetry [quick|full]
      julia --project=. scripts/regenerate_ignored_outputs.jl docs

    Notes:
      - Generated outputs are written under the ignored `output/` or `docs/build/` paths.
      - CI publishes telemetry CSVs as workflow artifacts instead of storing them in git.
    """)
end

function _run(cmd::Cmd)
    println("Running: ", cmd)
    run(cmd)
    return nothing
end

function _runtime_analysis_cmd(profile::String, outdir::String)
    launcher = joinpath(REPO_ROOT, "benchmarks", "studies", "performance_runtime_analysis.jl")
    return `$(Base.julia_cmd()) --startup-file=no --project=$REPO_ROOT $launcher $profile --outdir=$outdir`
end

function _telemetry_cmd(profile::String)
    launcher = joinpath(REPO_ROOT, "benchmarks", "studies", "telemetry_orbit_accuracy_study.jl")
    return `$(Base.julia_cmd()) --startup-file=no --project=$REPO_ROOT $launcher $profile --enforce=true`
end

function _docs_cmd()
    launcher = joinpath(REPO_ROOT, "docs", "make.jl")
    docs_project = joinpath(REPO_ROOT, "docs")
    return `$(Base.julia_cmd()) --startup-file=no --project=$docs_project $launcher`
end

function main(args::Vector{String})
    isempty(args) && return (_usage(); 1)

    mode = args[1]
    if mode == "runtime-analysis"
        profile = length(args) >= 2 ? args[2] : "quick"
        profile in ("quick", "full") || error("runtime-analysis profile must be `quick` or `full`.")
        outdir = length(args) >= 3 ? normpath(args[3]) : joinpath(REPO_ROOT, "output", "performance")
        mkpath(outdir)
        _run(_runtime_analysis_cmd(profile, outdir))
        return 0
    elseif mode == "telemetry"
        profile = length(args) >= 2 ? args[2] : "quick"
        profile in ("quick", "full") || error("telemetry profile must be `quick` or `full`.")
        mkpath(joinpath(REPO_ROOT, "output"))
        _run(_telemetry_cmd(profile))
        return 0
    elseif mode == "docs"
        _run(_docs_cmd())
        return 0
    else
        _usage()
        error("Unknown mode `$mode`.")
    end
end

exit(main(ARGS))
