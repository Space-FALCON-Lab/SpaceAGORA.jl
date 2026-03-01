#!/usr/bin/env julia

using SpaceAGORACalibration

@inline function _parse_bool(raw::AbstractString)::Bool
    token = lowercase(strip(raw))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid boolean token '$raw'."))
end

function parse_cli(args::Vector{String})
    spec_path = joinpath(@__DIR__, "..", "examples", "telemetry_hybrid_spec.toml")
    project_path = abspath(joinpath(@__DIR__, "..", ".."))
    verification_script = abspath(joinpath(project_path, "scripts", "verify_telemetry.jl"))
    manifest_path = ""
    profile = "quick"
    enforce = false
    plots = false

    for arg in args
        if startswith(arg, "--spec=")
            spec_path = abspath(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--project=")
            project_path = abspath(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--verification-script=")
            verification_script = abspath(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--manifest=")
            manifest_path = abspath(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--profile=")
            raw = lowercase(strip(split(arg, "=", limit=2)[2]))
            raw in ("quick", "full") || throw(ArgumentError("profile must be quick|full."))
            profile = raw
        elseif startswith(arg, "--enforce=")
            enforce = _parse_bool(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--plots=")
            plots = _parse_bool(split(arg, "=", limit=2)[2])
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: --spec=..., --project=..., --verification-script=..., --manifest=..., --profile=quick|full, --enforce=0|1, --plots=0|1"
            ))
        end
    end

    return (
        spec_path=spec_path,
        project_path=project_path,
        verification_script=verification_script,
        manifest_path=manifest_path,
        profile=profile,
        enforce=enforce,
        plots=plots
    )
end

function main(args::Vector{String})
    cfg = parse_cli(args)
    spec = load_spec(cfg.spec_path)

    backend = CommandBackend(
        julia_cmd=Base.julia_cmd(),
        project_path=cfg.project_path,
        verification_script=cfg.verification_script,
        manifest_path=cfg.manifest_path,
        profile=cfg.profile,
        enforce=cfg.enforce,
        plots=cfg.plots
    )

    result = run_calibration(spec, backend)
    println("run_id=" * result.run_id)
    println("run_dir=" * result.run_dir)
    println("report=" * result.report_path)
    println("best_candidate_id=" * string(result.best_candidate.id))
    println("best_score=" * string(result.best_score))
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main(copy(ARGS))
end
