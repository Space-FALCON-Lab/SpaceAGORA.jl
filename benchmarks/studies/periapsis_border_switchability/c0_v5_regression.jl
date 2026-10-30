# Estimator v5 Earth regression (the c0 pattern of c0_v4_regression.jl): the same Earth cells predicted by
#   v3 (drag on / drag off)   predict_history3
#   v5 (drag on / drag off)   predict_history5, default options (the carry rule is false at Earth: uncarried corrected stepper)
#   v5 legacy options         predict_history5 with window = :kepler, post_drag_period = false, kick_post_drag = false (must be identical to v3)
# R4 of estimator_v5_plan.md: where the timing fix does not matter (no drag) v5 == v3 to 1 s / 1 m over 60 d; the first 3 days of the drag runs agree to 5 m.
#   julia --project=. c0_v5_regression.jl --out=DIR --results-root=R
include(joinpath(@__DIR__, "estimator_v5.jl"))
using Statistics
const CLI0 = Dict(m.captures[1] => m.captures[2] for a in ARGS for m in [match(r"^--([^=]+)=(.*)$", a)] if m !== nothing)
const OUT0 = abspath(CLI0["out"]); const ROOT0 = abspath(CLI0["results-root"])
const STUDY = load_study(); const CFG2 = TOML.parsefile(joinpath(STUDY_DIR, "estimator_v2_config.toml"))
const SETT = settings_from(CFG2)
const OPTS = (v = CFG2["v2"]; V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])))
const TAB = SPTable(OPTS.n_sp)
cell_dir(tag, inc) = inc == 50.0 ? joinpath(ROOT0, "border", "grid", tag[7:end]) : joinpath(ROOT0, "border2", "grid", tag)
same_pred(p, q) = all(isequal(getfield(p, f), getfield(q, f)) for f in fieldnames(Prediction))

function maxdiff(p::Prediction, q::Prediction; tmax=Inf, kmax=typemax(Int))
    n = min(length(p.t_s), length(q.t_s)); n == 0 && return (NaN, NaN, 0)
    k = [i for i in 1:min(n, kmax) if p.t_s[i] <= tmax]; isempty(k) && return (NaN, NaN, 0)
    return maximum(abs.(p.t_s[k] .- q.t_s[k])), maximum(abs.(p.alt_m[k] .- q.alt_m[k])), abs(length(p.t_s) - length(q.t_s))
end

function main0()
    nbody_planet()
    rec = let r = earth_record(STUDY), cfg = load_corridor_config(CORRIDOR_CONFIG_PATH)
        PlanetRecord(; (f => getfield(r, f) for f in fieldnames(PlanetRecord) if f != :harm)...,
            harm=GravitationalHarmonicsModel(20, 20, joinpath(REPO_ROOT, String(cfg.values["gravity"]["coefficient_file"])), r.planet))
    end
    ctx = make_context(rec, SETT); zf = zonal_field(rec.harm); tf = tesseral_field(rec.harm); zo = ZonalOpts()
    Snd = with_settings(SETT; drag=false)
    legacy = TimingOpts(window=:kepler, post_drag_period=false, kick_post_drag=false, zonal_sp=false); to = TimingOpts(); tt = TesseralOpts()
    labels = CSV.read(joinpath(ROOT0, "oracle", "oracle_cell_labels.csv"), DataFrame)
    att = TOML.parsefile(joinpath(STUDY_DIR, "attribution_config.toml"))["cells"]["list"]
    atags = ["I$(lpad(round(Int, 10 * c[1]), 4, '0'))_R$(lpad(round(Int, c[2]), 3, '0'))_W$(lpad(round(Int, c[3]), 3, '0'))_E$(Int(c[4]))" for c in att]
    tags = unique(vcat(atags, labels.tag[1:28:end])); rows = NamedTuple[]
    for tag in tags
        c = labels[labels.tag .== tag, :][1, :]
        inc, raan, om = deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg); et0 = epoch_et(STUDY, Int(c.epoch))
        seed = CSV.read(joinpath(cell_dir(tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end]
        p3, _ = predict_history3(ctx, SETT, OPTS, TAB, zf, zo, seed, inc, raan, om, et0); p3n, _ = predict_history3(ctx, Snd, OPTS, TAB, zf, zo, seed, inc, raan, om, et0)
        p5, i5 = predict_history5(ctx, SETT, OPTS, TAB, zf, zo, tf, tt, to, seed, inc, raan, om, et0); p5n, _ = predict_history5(ctx, Snd, OPTS, TAB, zf, zo, tf, tt, to, seed, inc, raan, om, et0)
        p5l, _ = predict_history5(ctx, SETT, OPTS, TAB, zf, zo, tf, tt, legacy, seed, inc, raan, om, et0)
        dtn, dan, cn = maxdiff(p5n, p3n); dt3, da3, c3 = maxdiff(p5, p3; tmax=3 * 86400.0); dtp, dap, _ = maxdiff(p5, p3; kmax=3); dtf, daf, cf = maxdiff(p5, p3)
        push!(rows, (tag=tag, carried=i5[4], legacy_identical=same_pred(p3, p5l), nodrag_max_dt_s=dtn, nodrag_max_dalt_m=dan, nodrag_dpasses=cn,
            drag_first3passes_max_dt_s=dtp, drag_first3passes_max_dalt_m=dap, drag_first3d_max_dt_s=dt3, drag_first3d_max_dalt_m=da3, drag_full_max_dt_s=dtf, drag_full_max_dalt_m=daf, drag_full_dpasses=cf))
    end
    D = DataFrame(rows); mkpath(OUT0); CSV.write(joinpath(OUT0, "c0_v5_regression.csv"), D)
    println("cells: ", nrow(D), "; legacy options identical to v3 in all: ", all(D.legacy_identical), "; carried: ", sum(D.carried))
    println("no-drag run, v5 vs v3 over 60 d: max pass-time difference ", maximum(D.nodrag_max_dt_s), " s, max altitude difference ", maximum(D.nodrag_max_dalt_m), " m, pass-count differences ", sum(D.nodrag_dpasses))
    println("drag run, first 3 passes: max pass-time difference ", maximum(skipmissing(D.drag_first3passes_max_dt_s)), " s, max altitude difference ", maximum(skipmissing(D.drag_first3passes_max_dalt_m)), " m")
    println("drag run, first 3 d: max pass-time difference ", maximum(skipmissing(D.drag_first3d_max_dt_s)), " s, max altitude difference ", maximum(skipmissing(D.drag_first3d_max_dalt_m)), " m")
    println("drag run, full: median max pass-time difference ", median(D.drag_full_max_dt_s), " s (max ", maximum(D.drag_full_max_dt_s), "), median max altitude difference ", median(D.drag_full_max_dalt_m), " m (max ", maximum(D.drag_full_max_dalt_m), ")")
end
main0()
