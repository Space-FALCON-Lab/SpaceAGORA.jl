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
#   different_quantity - the reference is published and the model produces
#                        something close to it, but the two are not definitions
#                        of the same thing, so no tolerance can be argued in
#                        either direction. Reported with its ratio; never gates.
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
    status in ("sourced", "cross_regime", "calibration_target", "not_modeled", "unsourced",
               "different_quantity") ||
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
    elseif quantity == "ejection_angle_deg"
        # The mass-weighted ejection angle of the grain-transport distribution.
        # It is bounded by the configuration's own `ejection_angle_min_deg` and
        # `ejection_angle_max_deg`, which come from this very reference, so the
        # comparison is circular by construction and the case never gates; it is
        # evaluated so the circularity is visible rather than implied.
        model = _angle_probe_model(cfg)
        model === nothing && return nothing
        summary = plume_ejecta_summary(model, 1)
        return summary === nothing ? nothing : summary.mean_angle_deg
    elseif quantity == "erosion_radius_m"
        return plume_quantities(cfg, thrust_n, height_m).outer_m
    elseif quantity == "eroded_mass_kg"
        return nothing            # handled by `model_profile_mass`, needs the whole descent
    elseif quantity == "crater_depth_m"
        return nothing            # handled by `model_profile_crater`, needs the whole descent
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
reference averages with. It goes through `plume_mean_shear`, so it averages
whatever field the configuration carries rather than assuming a shape.
Reported alongside the peak shear; never used for a pass or fail.
"""
function model_mean_shear(cfg, thrust_n::Float64, height_m::Float64, radius_m::Float64)
    (isfinite(radius_m) && radius_m > 0.0) || return nothing
    return plume_mean_shear(cfg.field, cfg, thrust_n, height_m, radius_m)
end

"""
    study_config(field_name) -> (config, label)

The model configuration a run evaluates. `"analytic"` is `PlumeSurfaceConfig()`
with its default `PlumeAnalyticField`; `"table"` loads
`data/psi/apollo_lmde.json`; any other value is taken as a path to a plume field
table. Everything else about the configuration is the shipped default, so a run
reports the model as it is, not as it could be tuned.
"""
function study_config(field_name::AbstractString)
    name = strip(String(field_name))
    lowercase(name) == "analytic" && return PlumeSurfaceConfig(), "analytic (PlumeAnalyticField)"
    path = lowercase(name) == "table" ? joinpath(REPO_ROOT, "data", "psi", "apollo_lmde.json") : abspath(name)
    isfile(path) || error("No plume field table at $path (use --field=analytic, --field=table or a path).")
    return PlumeSurfaceConfig(field=load_plume_field(path)), "table ($(basename(path)))"
end

"""
    model_profile_crater(cfg, thrust_n, profile; bins=64, min_radius_m=0.05, max_radius_m=40.0)
        -> (depth_m, peak_radius_m, edge_radius_m, rows)

Depth of the crater the model digs over a descent given as `(t_s, height_m)`
samples: the erosion regimes' local mass flux at each radius, integrated in time
by the trapezoidal rule and divided by the soil's in-situ bulk density, on the
same logarithmic radial grid the effector carries. `depth_m` is the deepest
point, `peak_radius_m` the radius it sits at, and `edge_radius_m` the outermost
radius still at a tenth of it and at least `min_depth_m` deep — the effector's
own edge rule.

This is the same arithmetic `PlumeSurfaceInteractionModel` does at accepted
steps; it is repeated here so a case can be scored without running a
simulation, exactly as `model_profile_mass` repeats the mass integral. It goes
through the model's public API only.
"""
function model_profile_crater(cfg, thrust_n::Float64, profile::Vector{Any};
                              bins::Int=64, min_radius_m::Float64=0.05,
                              max_radius_m::Float64=40.0, gravity_m_s2::Float64=1.625,
                              min_depth_m::Float64=1.0e-3)
    length(profile) >= 2 || error("A crater case needs at least two samples.")
    ts = [Float64(row["t_s"]) for row in profile]
    hs = [Float64(row["height_m"]) for row in profile]
    order = sortperm(ts)
    ts, hs = ts[order], hs[order]
    radii = exp.(range(log(min_radius_m), log(max_radius_m); length=bins))
    local_rate(h) = begin
        _, R = plume_surface_footprint(cfg, thrust_n, h)
        env = erosion_environment(footprint_radius_m=max(R, cfg.nozzle_exit_radius_m),
                                  residence_time_s=cfg.gas_residence_time_s,
                                  bearing_width_m=2.0 * max(R, cfg.nozzle_exit_radius_m))
        [regolith_erosion_rate(cfg.regimes, plume_gas_state(cfg.field, cfg, thrust_n, h, r),
                               cfg.soil, gravity_m_s2, env).rate_kg_m2_s for r in radii]
    end
    depth = zeros(length(radii))
    previous = local_rate(hs[1])
    for k in 1:(length(ts) - 1)
        now = local_rate(hs[k + 1])
        depth .+= 0.5 .* (previous .+ now) .* (ts[k + 1] - ts[k]) ./ cfg.bulk_density_kg_m3
        previous = now
    end
    dmax, imax = findmax(depth)
    # The same edge rule the effector uses: a tenth of the deepest point, but
    # never shallower than the reporting floor.
    edge_depth = max(0.1 * dmax, min_depth_m)
    edge = 0.0
    for k in eachindex(radii)
        depth[k] >= edge_depth && (edge = radii[k])
    end
    rows = [(radius_m=radii[k], depth_m=depth[k]) for k in eachindex(radii)]
    return dmax, radii[imax], edge, rows
end

"""
    _angle_probe_model(cfg) -> Union{Nothing, PlumeSurfaceInteractionModel}

A bare `PlumeSurfaceInteractionModel` primed at one thrust and height, so the
ejecta distribution can be evaluated through the effector's own public entry
point without running a simulation. Returns `nothing` when the state cannot be
primed, which is how a condition below the erosion onset reports.
"""
_ANGLE_PROBE_CONTROL = (actuators = (thrust_n = [0.0],),)

function _angle_probe_model(cfg; thrust_n::Float64=12172.0, height_m::Float64=5.0)
    model = PlumeSurfaceInteractionModel(_ANGLE_PROBE_CONTROL; config=cfg, num_sats=1)
    q = plume_quantities(cfg, thrust_n, height_m)
    q.erosion_kg_s > 0.0 || return nothing
    st = model.state
    st.erosion_kg_s[1] = q.erosion_kg_s
    st.last_thrust_n[1] = thrust_n
    st.last_query_height_m[1] = height_m
    st.last_inner_m[1] = q.inner_m
    st.last_outer_m[1] = q.outer_m
    st.last_gravity_m_s2[1] = 1.625
    st.last_body_radius_m[1] = 1.7374e6
    return model
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
