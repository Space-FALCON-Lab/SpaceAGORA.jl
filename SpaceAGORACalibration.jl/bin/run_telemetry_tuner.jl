#!/usr/bin/env julia

using SpaceAGORACalibration
using SpaceAGORA

@inline function _parse_bool(raw::AbstractString)::Bool
    token = lowercase(strip(raw))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid boolean token '$raw'."))
end

@inline function _parse_profile(raw::AbstractString)::String
    token = lowercase(strip(raw))
    token in ("quick", "full") || throw(ArgumentError("profile must be quick|full."))
    return token
end

function _run_post_verify_best(
    run_dir::String;
    profile::String,
    plots::Bool,
    enforce::Bool
)
    profile_symbol = Symbol(_parse_profile(profile))
    best_manifest = joinpath(run_dir, "best_manifest.toml")
    isfile(best_manifest) || throw(ArgumentError("Missing best_manifest.toml at $best_manifest"))

    summary_path = joinpath(run_dir, "final_summary_$(profile).csv")
    errors_path = joinpath(run_dir, "final_errors_$(profile).csv")

    req = SpaceAGORA.VerificationRequest(
        profile=profile_symbol,
        out_summary=summary_path,
        out_errors=errors_path,
        manifest_path=best_manifest,
        enforce=enforce,
        generate_plots=plots
    )
    result = SpaceAGORA.run_verification(req)
    println("post_verify_profile=$(String(result.profile))")
    println("post_verify_summary=" * result.summary_path)
    println("post_verify_errors=" * result.errors_path)
    println("post_verify_plots=" * (isempty(result.plots_dir) ? "disabled" : result.plots_dir))
    return nothing
end

function parse_cli(args::Vector{String})
    spec_path = joinpath(@__DIR__, "..", "examples", "telemetry_hybrid_spec.toml")
    project_path = abspath(joinpath(@__DIR__, "..", ".."))
    verification_script = abspath(joinpath(project_path, "scripts", "verify_telemetry.jl"))
    manifest_path = ""
    profile = "quick"
    parallel_profile = nothing
    post_verify_best = false
    post_verify_profile = "full"
    post_verify_plots = true
    post_verify_enforce = false
    enforce = false
    plots = false
    backend_mode = "command"

    for arg in args
        if startswith(arg, "--spec=")
            spec_path = abspath(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--project=")
            project_path = abspath(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--backend=")
            raw = lowercase(strip(split(arg, "=", limit=2)[2]))
            raw in ("command", "inprocess") || throw(ArgumentError("backend must be command|inprocess."))
            backend_mode = raw
        elseif startswith(arg, "--verification-script=")
            verification_script = abspath(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--manifest=")
            manifest_path = abspath(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--profile=")
            profile = _parse_profile(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--parallel-profile=")
            raw = strip(split(arg, "=", limit=2)[2])
            SpaceAGORA.parse_parallel_profile(raw) # validate token early
            parallel_profile = raw
        elseif startswith(arg, "--post-verify-best=")
            post_verify_best = _parse_bool(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--post-verify-profile=")
            post_verify_profile = _parse_profile(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--post-verify-plots=")
            post_verify_plots = _parse_bool(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--post-verify-enforce=")
            post_verify_enforce = _parse_bool(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--enforce=")
            enforce = _parse_bool(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--plots=")
            plots = _parse_bool(split(arg, "=", limit=2)[2])
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: --spec=..., --project=..., --backend=command|inprocess, --verification-script=..., --manifest=..., --profile=quick|full, --parallel-profile=R0|R1_a|R1_b|R2|R3|R4|R4_full_auto, --post-verify-best=0|1, --post-verify-profile=quick|full, --post-verify-plots=0|1, --post-verify-enforce=0|1, --enforce=0|1, --plots=0|1"
            ))
        end
    end

    return (
        spec_path=spec_path,
        project_path=project_path,
        backend_mode=backend_mode,
        verification_script=verification_script,
        manifest_path=manifest_path,
        profile=profile,
        parallel_profile=parallel_profile,
        post_verify_best=post_verify_best,
        post_verify_profile=post_verify_profile,
        post_verify_plots=post_verify_plots,
        post_verify_enforce=post_verify_enforce,
        enforce=enforce,
        plots=plots
    )
end

function main(args::Vector{String})
    cfg = parse_cli(args)
    spec = load_spec(cfg.spec_path)

    backend = if cfg.backend_mode == "inprocess"
        InProcessBackend(
            manifest_path=cfg.manifest_path,
            profile=cfg.profile,
            parallel_profile=cfg.parallel_profile,
            enforce=cfg.enforce,
            plots=cfg.plots
        )
    else
        CommandBackend(
            julia_cmd=Base.julia_cmd(),
            project_path=cfg.project_path,
            verification_script=cfg.verification_script,
            manifest_path=cfg.manifest_path,
            profile=cfg.profile,
            parallel_profile=cfg.parallel_profile,
            enforce=cfg.enforce,
            plots=cfg.plots
        )
    end

    result = run_calibration(spec, backend)
    println("run_id=" * result.run_id)
    println("run_dir=" * result.run_dir)
    println("report=" * result.report_path)
    println("best_candidate_id=" * string(result.best_candidate.id))
    println("best_score=" * string(result.best_score))
    if cfg.post_verify_best
        _run_post_verify_best(
            result.run_dir;
            profile=cfg.post_verify_profile,
            plots=cfg.post_verify_plots,
            enforce=cfg.post_verify_enforce
        )
    end
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main(copy(ARGS))
end
