# Odyssey and Venus Express aerobraking reconstructions — benchmark of record.
#
# Runs the validation scenarios documented in
# docs/spaceagora_aerobraking_reconstruction_record.md through
# SpaceAGORA.TelemetryVerification, one frozen manifest per run
# (manifests/*.toml), writing per-run summary/error CSVs under
# results/<hostname>/<run>/.
#
# Usage:
#   julia --project=. benchmarks/studies/telemetry_validation/run_validation.jl \
#       [quick|full] [--runs=core|all|envelope|<name,...>] [--enforce=true|false] \
#       [--plots=true|false] [--out-root=DIR]
#
# The benchmark-of-record command is:
#   GRAM_ROOT="<rig>/GRAM Suite 2.0" julia --startup-file=no --depwarn=error \
#       --project=. benchmarks/studies/telemetry_validation/run_validation.jl \
#       full --runs=all --enforce=true
#
# Defaults: profile=quick, runs=core (odyssey_tolson + vex_venusgram),
# enforce=false, plots=false. `--enforce=true` gates only the record scenarios
# (odyssey_tolson, vex_venusgram, odyssey_marsgram); the +-1 sigma envelope
# runs are always informational, matching the record's caveat that the
# +1 sigma case impacts before the end of the campaign.

include(joinpath(@__DIR__, "common.jl"))

# Argument parsing happens before the (expensive) package loads so a bad
# invocation fails immediately, and so GRAMSuite is only pulled in when a
# selected run actually needs it. GRAM loading must stay at top level: the
# package-extension methods it triggers would not be visible to an
# already-running function body (world age).
const CLI_OPTS = parse_kv_args(copy(ARGS))
const PROFILE = Symbol(get(CLI_OPTS, "profile", "quick"))
PROFILE in (:quick, :full) || error("Unsupported profile '$PROFILE'. Use quick or full.")
const RUN_NAMES = resolve_runs(get(CLI_OPTS, "runs", "core"))
const ENFORCE = parse_bool_flag(get(CLI_OPTS, "enforce", "false"))
const GENERATE_PLOTS = parse_bool_flag(get(CLI_OPTS, "plots", "false"))
const OUT_ROOT = results_root(get(CLI_OPTS, "out-root", ""))
const SPECS = [run_spec(name) for name in RUN_NAMES]

if any(spec -> spec.needs_gram, SPECS)
    println("Loading GRAMSuite (required by: " *
            join([spec.name for spec in SPECS if spec.needs_gram], ", ") * ") ...")
    load_gramsuite!()
end

using SpaceAGORA
using SpaceAGORA.TelemetryVerification
using DataFrames
using Printf

const ENFORCED_RUNS = ("odyssey_tolson", "vex_venusgram", "odyssey_marsgram")

function main()
    println("Telemetry validation study — Odyssey / Venus Express reconstructions of record")
    println("record   = $RECORD_DOC")
    println("profile  = $PROFILE")
    println("runs     = $(join(RUN_NAMES, ", "))")
    println("enforce  = $ENFORCE (envelope runs are never enforced)")
    println("out root = $OUT_ROOT")

    failures = String[]
    summaries = Any[]
    for spec in SPECS
        run_dir = joinpath(OUT_ROOT, spec.name)
        mkpath(run_dir)
        run_enforce = ENFORCE && spec.name in ENFORCED_RUNS
        println("\n=== $(spec.name): $(spec.description)")
        request = VerificationRequest(
            profile=PROFILE,
            out_summary=joinpath(run_dir, "telemetry_orbit_accuracy_summary.csv"),
            out_errors=joinpath(run_dir, "telemetry_orbit_accuracy_errors.csv"),
            manifest_path=spec.manifest,
            enforce=run_enforce,
            generate_plots=GENERATE_PLOTS
        )
        try
            result = run_verification(request)
            push!(summaries, (name=spec.name, summary=result.summary))
        catch err
            # Keep going so one gated failure still produces the other runs'
            # artifacts; the process exits nonzero below.
            @error "Run $(spec.name) failed" exception=(err, catch_backtrace())
            push!(failures, spec.name)
        end
    end

    if !isempty(summaries)
        println("\n================ Record comparison ================")
        println("Recorded full-profile values ($RECORD_DOC):")
        println("  odyssey_tolson   apo MAE 162.60 km nMAE 0.00667 | peri MAE 4.16 km | decay ratio 0.984")
        println("  odyssey_marsgram apo MAE 1399.52 km nMAE 0.05745 | peri MAE 2.53 km | decay ratio 0.840")
        println("  vex_venusgram    apo MAE 254.34 km nMAE 0.0701  | peri MAE 0.422 km | decay 1.262 median")
        println("This run ($PROFILE profile):")
        for entry in summaries
            for row in eachrow(entry.summary)
                decay = row.event == "apo" ?
                    @sprintf(" | decay ratio %.3f median / %.3f total",
                             row.drag_decay_ratio_median, row.drag_decay_ratio_total) : ""
                @printf(
                    "  %-28s %-4s MAE %10.3f km  RMSE %10.3f km  max %10.3f km  nMAE %.5f  bias %+9.3f km  Cr %.2f%s\n",
                    entry.name, row.event, row.mae_km, row.rmse_km, row.max_abs_km,
                    row.nmae, row.calibrated_bias_km, row.calibrated_cr, decay
                )
            end
        end
        PROFILE == :quick && println(
            "NOTE: quick profile propagates only the first orbits; " *
            "record values are only reproduced by the full profile."
        )
    end

    if !isempty(failures)
        error("Telemetry validation run(s) failed: $(join(failures, ", "))")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
