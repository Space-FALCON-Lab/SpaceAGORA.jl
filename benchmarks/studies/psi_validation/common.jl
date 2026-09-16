# Shared helpers for the plume-surface interaction (PSI) validation study.
# See README.md in this directory for the reference cases and their sources.
#
# Every number this study compares against comes from a manifest under
# `manifests/`, one per reference case, and every manifest names its source
# precisely (authors, title, venue, year, and the table/figure/equation the
# number is on). Nothing is hard-coded here: this file only knows how to read a
# manifest, drive the model's public API at the case's conditions, and score the
# result against the tolerance the manifest declares.

using TOML
using Printf

const STUDY_DIR = normpath(@__DIR__)
const REPO_ROOT = normpath(joinpath(STUDY_DIR, "..", "..", ".."))
const MANIFEST_DIR = joinpath(STUDY_DIR, "manifests")

# Status of a case, which decides whether it can gate anything:
#   sourced            - the reference value is published and the model produces
#                        the same quantity, so the case can pass or fail.
#   cross_regime       - the reference value is published and the model produces
#                        a comparable quantity, but the measurement comes from a
#                        different body, atmosphere or flow regime and its own
#                        source cautions against transferring it. Reported with
#                        its ratio; never gates.
#   calibration_target - the reference was used to fit a model parameter, so
#                        agreement is circular and proves nothing. Never gates.
#   not_modeled        - the reference is published but the model has no output
#                        of that kind. Recorded so the gap is visible. Never gates.
#   unsourced          - no published value was found for the model's behavior.
#                        Never gates; the model's value is an assumption.
const GATEABLE_STATUS = "sourced"

results_root(out_root::AbstractString="") = isempty(out_root) ?
    joinpath(STUDY_DIR, "results", gethostname()) : abspath(out_root)

"""
    manifest_paths() -> Vector{String}

Every case manifest in `manifests/`, in file-name order, so a run's output is
deterministic.
"""
manifest_paths() = sort(filter(p -> endswith(p, ".toml"), readdir(MANIFEST_DIR; join=true)))

"""
    load_case(path) -> NamedTuple

Parse one case manifest. The manifest is the single source of truth for the
reference values, their provenance and the tolerance; this function only checks
that the parts the harness needs are present.
"""
function load_case(path::AbstractString)
    raw = TOML.parsefile(path)
    get(raw, "version", 0) == 1 || error("$(basename(path)): unsupported manifest version")
    case = raw["case"]
    source = raw["source"]
    model = raw["model"]
    status = String(case["status"])
    status in ("sourced", "cross_regime", "calibration_target", "not_modeled", "unsourced") ||
        error("$(basename(path)): unknown case status '$status'")
    gate = Bool(case["gate_eligible"])
    if gate && status != GATEABLE_STATUS
        error("$(basename(path)): only '$GATEABLE_STATUS' cases may be gate-eligible")
    end
    return (
        id=String(case["id"]),
        title=String(case["title"]),
        status=status,
        gate_eligible=gate,
        summary=String(get(case, "summary", "")),
        source=source,
        conditions=get(raw, "conditions", Dict{String, Any}()),
        model=model,
        tolerance=get(raw, "tolerance", Dict{String, Any}()),
        references=Vector{Any}(get(raw, "reference", Any[])),
        profile=Vector{Any}(get(raw, "profile", Any[])),
        path=path,
    )
end

load_cases() = [load_case(p) for p in manifest_paths()]

function select_cases(cases::Vector, selector::AbstractString)
    sel = strip(selector)
    lowercase(sel) == "all" && return cases
    if lowercase(sel) == "gates"
        return [c for c in cases if c.gate_eligible]
    end
    wanted = [String(strip(tok)) for tok in split(sel, ",") if !isempty(strip(tok))]
    isempty(wanted) && error("No cases selected (got '$selector').")
    out = similar(cases, 0)
    for name in wanted
        idx = findfirst(c -> c.id == name, cases)
        idx === nothing && error(
            "Unknown case '$name'. Available: " * join([c.id for c in cases], ", ") *
            " (or the aliases all, gates)."
        )
        push!(out, cases[idx])
    end
    return out
end

# ---- driving the model ---------------------------------------------------------------
#
# Everything below goes through the model's public API only
# (`PlumeSurfaceConfig`, `plume_surface_footprint`, `plume_quantities`,
# `plume_erosion_onset_height`, `plume_ground_effect_force`), so the harness keeps
# working when the effector's internals change.

"""
    model_quantity(cfg, quantity, thrust_n, height_m) -> Union{Float64, Nothing}

One model output named by a manifest's `model.quantity`. `nothing` means the
model has no output of that kind, which is how a `not_modeled` case reports.
"""
function model_quantity(cfg, quantity::AbstractString, thrust_n::Float64, height_m::Float64)
    if quantity == "erosion_kg_s"
        return plume_quantities(cfg, thrust_n, height_m).erosion_kg_s
    elseif quantity == "scour_radius_m"
        return plume_quantities(cfg, thrust_n, height_m).outer_m
    elseif quantity == "peak_shear_pa"
        return plume_quantities(cfg, thrust_n, height_m).shear_pa
    elseif quantity == "peak_pressure_pa"
        return plume_quantities(cfg, thrust_n, height_m).pressure_pa
    elseif quantity == "ejecta_speed_mps"
        return plume_quantities(cfg, thrust_n, height_m).ejecta_mps
    elseif quantity == "footprint_radius_m"
        return plume_surface_footprint(cfg, thrust_n, height_m)[2]
    elseif quantity == "onset_height_m"
        return plume_erosion_onset_height(cfg, thrust_n)
    elseif quantity == "threshold_shear_pa"
        # The shear stress at which the model starts eroding, read off the model
        # itself rather than off a config field: the peak wall shear at the
        # erosion onset height is the threshold by the model's own definition,
        # whatever the internals are.
        onset = plume_erosion_onset_height(cfg, thrust_n)
        onset > 0.0 || return nothing
        return plume_quantities(cfg, thrust_n, onset).shear_pa
    elseif quantity == "ground_effect_fraction"
        thrust_n > 0.0 || return nothing
        return plume_ground_effect_force(cfg, thrust_n, height_m) / thrust_n
    elseif quantity == "eroded_mass_kg"
        return nothing            # handled by `model_profile_mass`, needs the whole descent
    elseif quantity == "none"
        return nothing
    end
    error("Unknown model quantity '$quantity'.")
end

"""
    model_profile_mass(cfg, thrust_n, profile) -> (total_kg, rows)

Total eroded mass the model predicts over a descent given as `(t_s, height_m)`
samples, by trapezoidal integration of the model's erosion rate in time. The
rows are returned so the run's output can show the rate at every sample.

`profile` rows are ordered by time; `t_s` is negative before touchdown, matching
the reference profile's own convention.
"""
function model_profile_mass(cfg, thrust_n::Float64, profile::Vector{Any})
    length(profile) >= 2 || error("A profile case needs at least two samples.")
    ts = [Float64(row["t_s"]) for row in profile]
    hs = [Float64(row["height_m"]) for row in profile]
    order = sortperm(ts)
    ts, hs = ts[order], hs[order]
    rates = [plume_quantities(cfg, thrust_n, h).erosion_kg_s for h in hs]
    total = 0.0
    for k in 1:(length(ts) - 1)
        total += 0.5 * (rates[k] + rates[k + 1]) * (ts[k + 1] - ts[k])
    end
    return total, [(t_s=ts[k], height_m=hs[k], rate_kg_s=rates[k]) for k in eachindex(ts)]
end

"""
    model_mean_shear(cfg, thrust_n, height_m, radius_m) -> Union{Float64, Nothing}

Supplementary diagnostic: the model's wall shear stress averaged over a disc of
`radius_m`, `2/R² ∫₀^R τ(r) r dr`, which is the definition the Apollo CFD
reference averages with. It reconstructs the model's shear profile from
`plume_surface_footprint` and the config's documented `friction_coefficient`,
so it returns `nothing` if a future field model no longer exposes that shape.
Reported alongside the peak shear; never used for a pass or fail.
"""
function model_mean_shear(cfg, thrust_n::Float64, height_m::Float64, radius_m::Float64)
    hasproperty(cfg, :friction_coefficient) || return nothing
    (isfinite(radius_m) && radius_m > 0.0) || return nothing
    p0, R = plume_surface_footprint(cfg, thrust_n, height_m)
    (isfinite(p0) && isfinite(R) && R > 0.0) || return nothing
    cf = Float64(getproperty(cfg, :friction_coefficient))
    n = 2048
    dr = radius_m / n
    acc = 0.0
    for k in 1:n
        r = (k - 0.5) * dr
        x = r / R
        tau = cf * p0 * 2.0 * x * exp(-x * x)
        acc += tau * r * dr
    end
    return 2.0 * acc / (radius_m * radius_m)
end

# ---- scoring -------------------------------------------------------------------------

"""
    verdict(case, ratio) -> String

`pass`/`fail` for a sourced case whose ratio falls inside (outside) the
manifest's tolerance band; otherwise the case's own status, which never gates.
"""
function verdict(case, ratio::Union{Float64, Nothing})
    case.status == GATEABLE_STATUS || return case.status
    ratio === nothing && return "no-model-output"
    isfinite(ratio) || return "fail"
    lo = Float64(case.tolerance["low"])
    hi = Float64(case.tolerance["high"])
    return (lo <= ratio <= hi) ? "pass" : "fail"
end

_fmt(x::Nothing) = ""
_fmt(x::Float64) = isfinite(x) ? @sprintf("%.4g", x) : string(x)
_fmt(x) = string(x)

"""
    csv_escape(s) -> String

Minimal RFC 4180 quoting, so source strings with commas survive the CSV.
"""
function csv_escape(s)
    str = string(s)
    if occursin(',', str) || occursin('"', str) || occursin('\n', str)
        return '"' * replace(str, '"' => "\"\"") * '"'
    end
    return str
end

write_csv(path::AbstractString, header::Vector{String}, rows::Vector{Vector{String}}) = open(path, "w") do io
    println(io, join(header, ","))
    for row in rows
        println(io, join((csv_escape(c) for c in row), ","))
    end
end

function parse_kv_args(argv::Vector{String})
    opts = Dict{String, String}()
    for arg in argv
        startswith(arg, "--") && occursin("=", arg) ||
            error("Unsupported argument '$arg'. Use --key=value flags.")
        k, v = split(arg[3:end], "="; limit=2)
        opts[String(k)] = String(v)
    end
    return opts
end

parse_bool_flag(raw::AbstractString) = lowercase(strip(raw)) in ("1", "true", "yes", "on")

"""
    source_line(case) -> String

The case's source as one citation line: authors, title, venue, year, and the
table/figure/equation the number is on.
"""
function source_line(case)
    s = case.source
    parts = [String(s["authors"]), "\"" * String(s["title"]) * "\"", String(s["venue"]),
             string(s["year"]), String(s["locator"])]
    line = join(parts, ", ")
    doi = String(get(s, "doi", ""))
    isempty(doi) || (line *= ", doi:" * doi)
    return line
end
