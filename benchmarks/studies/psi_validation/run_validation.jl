# Plume-surface interaction validation study.
#
# Evaluates SpaceAGORA's plume-surface interaction model at the conditions of
# every reference case under `manifests/` and writes a model-versus-reference
# comparison to `results/<hostname>/`. Each manifest names its source precisely
# and declares its own tolerance and whether it may gate; this script only
# drives the model's public API and scores the result.
#
# Usage:
#   julia --project=. benchmarks/studies/psi_validation/run_validation.jl \
#       [--cases=all|gates|<id,...>] [--enforce=true|false] [--out-root=DIR] \
#       [--field=analytic|table|<path to a plume field json>]
#
# Defaults: cases=all, enforce=false, field=analytic. The results land under
# `results/<hostname>/<analytic|table>/`, so the two fields can be compared. `--enforce=true` exits nonzero if any
# gate-eligible case falls outside its manifest's tolerance; cases that are
# calibration targets, cross-regime, unsourced or not produced by the model are
# never enforced, whatever the flag says.
#
# The run is cheap: it evaluates closed-form functions, loads no kernels and
# runs no simulation.

include(joinpath(@__DIR__, "common.jl"))

const CLI_OPTS = parse_kv_args(copy(ARGS))
const CASE_SELECTOR = get(CLI_OPTS, "cases", "all")
const ENFORCE = parse_bool_flag(get(CLI_OPTS, "enforce", "false"))
const FIELD_NAME = get(CLI_OPTS, "field", "analytic")
const OUT_ROOT = joinpath(results_root(get(CLI_OPTS, "out-root", "")), lowercase(FIELD_NAME) == "analytic" ? "analytic" : "table")

using SpaceAGORA

"""
    evaluate(case, cfg) -> Vector{NamedTuple}

One comparison row per reference row of the case.
"""
function evaluate(case, cfg)
    quantity = String(case.model["quantity"])
    thrust = Float64(case.conditions["thrust_n"])
    want_mean_shear = Bool(get(case.model, "mean_shear_over_reference_radius", false))
    rows = NamedTuple[]

    if quantity == "crater_depth_m"
        depth, peak_r, edge_r, _ = model_profile_crater(cfg, thrust, case.profile)
        for (k, ref) in enumerate(case.references)
            refval = Float64(ref["value"])
            ratio = refval > 0.0 ? depth / refval : NaN
            push!(rows, (
                height_m=nothing, reference=refval, model=depth, ratio=ratio,
                context_only=Bool(get(ref, "context_only", false)),
                primary=(k == 1), mean_shear=nothing,
                note=String(get(ref, "note", "")),
            ))
        end
        return rows, [(t_s=NaN, height_m=peak_r, rate_kg_s=edge_r)]
    end

    if quantity == "eroded_mass_kg"
        total, profile_rows = model_profile_mass(cfg, thrust, case.profile)
        for (k, ref) in enumerate(case.references)
            refval = Float64(ref["value"])
            ratio = refval > 0.0 ? total / refval : NaN
            push!(rows, (
                height_m=nothing, reference=refval, model=total, ratio=ratio,
                context_only=Bool(get(ref, "context_only", false)),
                primary=(k == 1), mean_shear=nothing,
                note=String(get(ref, "note", "")),
            ))
        end
        return rows, profile_rows
    end

    for (k, ref) in enumerate(case.references)
        h = haskey(ref, "height_m") ? Float64(ref["height_m"]) : 0.0
        modelval = model_quantity(cfg, quantity, thrust, h)
        refval = Float64(ref["value"])
        unsourced = Bool(get(ref, "unsourced", false))
        ratio = (modelval === nothing || unsourced || refval == 0.0) ? nothing : modelval / refval
        mean_shear = nothing
        if want_mean_shear && haskey(ref, "mean_radius_m")
            mean_shear = model_mean_shear(cfg, thrust, h, Float64(ref["mean_radius_m"]))
        end
        push!(rows, (
            height_m=haskey(ref, "height_m") ? h : nothing,
            reference=(unsourced ? nothing : refval), model=modelval, ratio=ratio,
            context_only=Bool(get(ref, "context_only", false)),
            primary=(k == 1), mean_shear=mean_shear,
            note=String(get(ref, "note", "")),
        ))
    end
    return rows, NamedTuple[]
end

"""
    case_verdict(case, rows) -> String

A case passes when every non-context reference row is inside the tolerance.
Context rows are printed for perspective and never decide anything.
"""
function case_verdict(case, rows)
    case.status == GATEABLE_STATUS || return case.status
    scored = [r for r in rows if !r.context_only]
    isempty(scored) && return "no-scored-rows"
    for r in scored
        verdict(case, r.ratio) == "pass" || return "fail"
    end
    return "pass"
end


"""
    print_findings(evaluated)

The run's own findings, computed from the rows it just produced rather than
written down in advance: where each scored case sits relative to its reference,
whether the model's trend with height has the same sign as the reference's, and
what is on record as unmeasurable. The prose discussion of what these mean is in
this directory's README, which is written from an actual run of this script.
"""
function print_findings(evaluated)
    println("\n---------------- Findings from this run ----------------")
    for (case, rows, v) in evaluated
        scored = [r for r in rows if !r.context_only && r.ratio !== nothing && isfinite(r.ratio)]
        if isempty(scored)
            reason = case.status == "not_modeled" ? "the model produces no comparable output" :
                     case.status == "different_quantity" ? "the model's output is not a definition of the same quantity" :
                     case.status == "unsourced" ? "no published reference value exists" :
                     "no scored row produced a ratio"
            println("  $(case.id): no ratio -- $reason.")
            continue
        end
        ratios = [r.ratio for r in scored]
        lo, hi = minimum(ratios), maximum(ratios)
        direction = hi < 1.0 ? "below" : (lo > 1.0 ? "above" : "straddling")
        @printf("  %-32s model/reference %.3f to %.3f (%s the reference), %s\n",
                case.id * ":", lo, hi, direction,
                case.status == GATEABLE_STATUS ? v : "not gated ($(case.status))")
        # A trend check the tolerance band cannot make: over a height series, do
        # the model and the reference move the same way as the vehicle descends?
        heights = [r.height_m for r in scored]
        if length(scored) >= 3 && all(h -> h !== nothing, heights)
            order = sortperm(Float64[h for h in heights])
            mv = [scored[i].model for i in order]
            rv = [scored[i].reference for i in order]
            if all(x -> x !== nothing, mv) && all(x -> x !== nothing, rv)
                model_falls_with_height = last(mv) < first(mv)
                ref_falls_with_height = last(rv) < first(rv)
                agree = model_falls_with_height == ref_falls_with_height
                @printf("  %-32s trend with height: model %s, reference %s -- %s\n", "",
                        model_falls_with_height ? "falls" : "rises",
                        ref_falls_with_height ? "falls" : "rises",
                        agree ? "same sign" : "OPPOSITE SIGN")
                # Where each series peaks. A model whose extremum sits at a
                # different height than the reference's has a shape error the
                # ratio band cannot see.
                mh = Float64[heights[i] for i in order]
                mi, ri = argmax(Float64[x for x in mv]), argmax(Float64[x for x in rv])
                if mi != ri
                    @printf("  %-32s peaks at a different height: model at %.2f m, reference at %.2f m\n",
                            "", mh[mi], mh[ri])
                end
            end
        end
    end
    missing_outputs = [c.id for (c, _, _) in evaluated if c.status == "not_modeled"]
    unsourced = [c.id for (c, _, _) in evaluated if c.status == "unsourced"]
    isempty(missing_outputs) ||
        println("  published quantities the model cannot produce: " * join(missing_outputs, ", "))
    isempty(unsourced) ||
        println("  model outputs with no published reference at all: " * join(unsourced, ", "))
    println("  the discussion of what this means is in benchmarks/studies/psi_validation/README.md")
    return nothing
end

function main()
    cases = select_cases(load_cases(), CASE_SELECTOR)
    cfg, field_label = study_config(FIELD_NAME)
    mkpath(OUT_ROOT)

    println("Plume-surface interaction validation study")
    println("cases    = $(join([c.id for c in cases], ", "))")
    println("enforce  = $ENFORCE (only 'sourced' cases can be enforced)")
    println("out root = $OUT_ROOT")
    println("model    = PlumeSurfaceConfig() defaults (Apollo LM DPS over lunar mare regolith)")
    println("field    = $field_label")
    println("erosion  = $(cfg.erosion_model), regimes $(map(typeof, cfg.regimes))")

    comparison_rows = Vector{String}[]
    case_rows = Vector{String}[]
    evaluated = Tuple{Any, Vector{NamedTuple}, String}[]
    failures = String[]
    gate_pass = 0
    gate_total = 0

    for case in cases
        rows, profile_rows = evaluate(case, cfg)
        v = case_verdict(case, rows)
        push!(evaluated, (case, rows, v))
        unit = String(case.model["unit"])
        println("\n=== $(case.id) [$(case.status)$(case.gate_eligible ? ", gate-eligible" : "")] -> $v")
        println("    $(case.title)")
        println("    source: $(source_line(case))")
        if case.status == GATEABLE_STATUS
            gate_total += 1
            v == "pass" && (gate_pass += 1)
        end

        println("    unit: $unit")
        @printf("    %10s  %12s  %12s  %9s  %s\n",
                "height[m]", "reference", "model", "model/ref", "row")
        for r in rows
            hcol = r.height_m === nothing ? "-" : @sprintf("%.2f", r.height_m)
            kind = r.context_only ? "context" : "scored"
            extra = r.mean_shear === nothing ? "" :
                @sprintf("  [model area-avg over a0: %.4g %s]", r.mean_shear, unit)
            refcol = r.reference === nothing ? "(none)" : _fmt(r.reference)
            modcol = r.model === nothing ? "(none)" : _fmt(r.model)
            @printf("    %10s  %12s  %12s  %9s  %s%s\n",
                    hcol, refcol, modcol,
                    r.ratio === nothing ? "-" : @sprintf("%.3f", r.ratio), kind, extra)
            push!(comparison_rows, [
                case.id, case.status, string(case.gate_eligible), String(case.model["quantity"]),
                unit, hcol, _fmt(r.reference), _fmt(r.model),
                r.ratio === nothing ? "" : @sprintf("%.6g", r.ratio),
                r.mean_shear === nothing ? "" : @sprintf("%.6g", r.mean_shear),
                r.context_only ? "context" : "scored",
                r.context_only ? "" : verdict(case, r.ratio),
                r.note,
            ])
        end
        if !isempty(profile_rows)
            if String(case.model["quantity"]) == "crater_depth_m"
                @printf("    crater shape: deepest point at r = %.2f m, edge (a tenth of the peak) at r = %.2f m\n",
                        profile_rows[1].height_m, profile_rows[1].rate_kg_s)
            else
                println("    descent profile the model was integrated over:")
                for p in profile_rows
                    @printf("      t %+7.1f s  h %6.2f m  model rate %8.3f kg/s\n",
                            p.t_s, p.height_m, p.rate_kg_s)
                end
            end
        end
        if case.status == GATEABLE_STATUS && v != "pass"
            push!(failures, case.id)
            println("    tolerance: $(case.tolerance["low"]) to $(case.tolerance["high"]) (model/reference)")
        end

        push!(case_rows, [
            case.id, case.title, case.status, string(case.gate_eligible), v,
            String(case.model["quantity"]), unit,
            haskey(case.conditions, "thrust_n") ? string(case.conditions["thrust_n"]) : "",
            string(get(case.tolerance, "low", "")), string(get(case.tolerance, "high", "")),
            source_line(case), String(get(case.source, "url", "")),
        ])
    end

    comparison_path = joinpath(OUT_ROOT, "psi_validation_comparison.csv")
    cases_path = joinpath(OUT_ROOT, "psi_validation_cases.csv")
    write_csv(comparison_path,
              ["case", "status", "gate_eligible", "quantity", "unit", "height_m", "reference",
               "model", "ratio", "model_area_avg", "row_kind", "verdict", "note"],
              comparison_rows)
    write_csv(cases_path,
              ["case", "title", "status", "gate_eligible", "verdict", "quantity", "unit",
               "thrust_n", "tolerance_low", "tolerance_high", "source", "url"],
              case_rows)

    println("\n================ Summary ================")
    @printf("gate-eligible cases: %d of %d inside tolerance\n", gate_pass, gate_total)
    for case in cases
        case.status == GATEABLE_STATUS && continue
        println("  not gated ($(case.status)): $(case.id)")
    end
    print_findings(evaluated)
    println("wrote $comparison_path")
    println("wrote $cases_path")

    if !isempty(failures)
        msg = "gate-eligible case(s) outside tolerance: $(join(failures, ", "))"
        if ENFORCE
            error(msg)
        else
            println("NOTE: $msg (not enforced; pass --enforce=true to fail the run)")
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
