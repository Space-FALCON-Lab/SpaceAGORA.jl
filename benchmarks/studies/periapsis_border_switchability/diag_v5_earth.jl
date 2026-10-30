# Pass-timing error budget at Earth by truth-state injection (DIAGNOSTIC for estimator v5; never enters a verdict). The Earth analog of diag_v5_timing.jl: the reference is a Cowell
# leg of the estimator's own model class (pole = J2000 z, zonal gravity 2..20 of the truth file, Sun and Moon from SPICE, US76 fit drag), so the stepper is tested against its own
# physics; the Cowell leg itself is compared with the engine truth (periapsis time of pass k+1 from the truth state of pass k) to show the size of what the model class omits
# (non-zonal field, Earth orientation). Usage as diag_v5_timing.jl, plus --root and --study (default study_config.toml).
using Statistics
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "estimator_v5.jl"))
include(joinpath(@__DIR__, "diag_v5_lib.jl"))
function parse_cli_e()
    kv = Dict{String, String}()
    for a in ARGS
        m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2])
    end
    return kv
end
const CLIE = parse_cli_e()
const STUDYE = load_study(abspath(get(CLIE, "study", STUDY_CONFIG_PATH)))
const ROOTE = abspath(CLIE["root"])
const CFG2E = TOML.parsefile(joinpath(STUDY_DIR, "estimator_v2_config.toml"))
const SETTE = settings_from(CFG2E)
const OPTSE = (v = CFG2E["v2"]; V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])))
const TABE = SPTable(OPTSE.n_sp)
cell_dir_e(tag, inc) = inc == 50.0 ? joinpath(ROOTE, "border", "grid", tag[7:end]) : joinpath(ROOTE, "border2", "grid", tag)
function earth_record_deg(study, degree::Int, order::Int)
    rec = earth_record(study); cfg = load_corridor_config(CORRIDOR_CONFIG_PATH)
    harm = GravitationalHarmonicsModel(degree, order, joinpath(REPO_ROOT, String(cfg.values["gravity"]["coefficient_file"])), rec.planet)
    return PlanetRecord(; (f => getfield(rec, f) for f in fieldnames(PlanetRecord) if f != :harm)..., harm=harm)
end
with_J2(c::EstContext, J2) = EstContext(c.mu, c.Re, J2, c.Rp_e, c.planet, c.omega_e, c.fit, c.beta, c.q_low, c.q_high, c.q_hard, c.w2m, c.apo_alt_m, c.body_names, c.body_mu, c.gl_x, c.gl_w, c.center, c.rot)

const VARS_E = Dict(
    "full"    => (j2=true,  sun=true,  drag=true,  zonal=:rule, grav=:zonal),
    "nodrag"  => (j2=true,  sun=true,  drag=false, zonal=:rule, grav=:zonal),
    "zonaldrag" => (j2=true, sun=false, drag=true,  zonal=:rule, grav=:zonal),
    "j2sundrag" => (j2=true, sun=true,  drag=true,  zonal=:none, grav=:j2),
    "j2drag"  => (j2=true,  sun=false, drag=true,  zonal=:none, grav=:j2),
    "j2only"  => (j2=true,  sun=false, drag=false, zonal=:none, grav=:j2),
    "drag2b"  => (j2=false, sun=false, drag=true,  zonal=:none, grav=:point),
)
function stepper_e(name, to_kw)
    tt = TesseralOpts(force=false)
    name == "v3" && return (ctx, S, zf, zo, tf, apo_el, t_apo, et0, np) -> predict_history3(ctx, S, OPTSE, TABE, zf, zo, NaN, 0.0, 0.0, 0.0, et0 + t_apo; max_passes=np, init_override=apo_el)
    to = name == "v5" ? TimingOpts(; to_kw...) : throw(ArgumentError("stepper $name"))
    return (ctx, S, zf, zo, tf, apo_el, t_apo, et0, np) -> predict_history5(ctx, S, OPTSE, TABE, zf, zo, tf, tt, to, NaN, 0.0, 0.0, 0.0, et0 + t_apo; max_passes=np, init_override=apo_el)
end

function main_e()
    nbody_planet(); rec = earth_record_deg(STUDYE, 20, 20); zf = zonal_field(rec.harm); tf = tesseral_field(rec.harm); ctx = make_context(rec, SETTE); cm = CowellModel(ctx, zf, tf)
    labels = CSV.read(joinpath(ROOTE, "oracle", "oracle_cell_labels.csv"), DataFrame)
    tags = split(CLIE["cells"], ","); ks = parse.(Int, split(get(CLIE, "ks", "1,2,3,5,8,12,16,20,25,30,40,50"), ","))
    vns = String.(split(get(CLIE, "variants", "full,nodrag,j2drag,drag2b,zonaldrag,j2sundrag"), ",")); npass = parse(Int, get(CLIE, "npass", "30")); mode = get(CLIE, "mode", "inject")
    stp = stepper_e(get(CLIE, "stepper", "v5"), (zonal_sp=false,))
    rows = NamedTuple[]
    for tag in tags
        c = labels[labels.tag .== tag, :][1, :]; et0 = epoch_et(STUDYE, Int(c.epoch))
        tr = CSV.read(joinpath(cell_dir_e(tag, c.inc_deg), "passes_3b.csv"), DataFrame)
        for k in (mode == "chain" ? [1] : ks)
            k + 2 <= nrow(tr) || continue
            el = (a=tr.a_m[k], e=tr.e[k], i=tr.i_rad[k], raan=tr.raan_rad[k], omega=tr.omega_rad[k], nu=0.0)
            r0, v0 = elements_state(el, ctx.mu); t0 = tr.time_s[k]
            span = mode == "chain" ? (npass + 2) * 90000.0 : 2.5 * 90000.0
            bodies = sun_track(ctx, et0, t0, t0 + span; dt=3600.0)
            ev0 = cowell_leg(cm, CowellForces(grav=:zonal), r0, v0, t0, 0.0, 0.0, bodies, mode == "chain" ? 2 : 3; start=:peri)
            # context for the Cowell-vs-truth comparison (non-zonal field and Earth orientation omitted)
            valid_dt = ev0[2].t - tr.time_s[k + 1]; valid_dr = norm(ev0[2].r) - tr.radius_m[k + 1]
            t_apo = ev0[1].t; apo_el = state_elements(ev0[1].r, ev0[1].v, ctx.mu)
            for vn in vns
                V = VARS_E[vn]; ctxv = V.j2 ? ctx : with_J2(ctx, 0.0); Sv = with_settings(SETTE; tides=V.sun, drag=V.drag, j2=V.j2)
                zo = V.zonal === :none ? ZonalOpts(force=Int[]) : ZonalOpts(); cmv = CowellModel(ctxv, zf, tf, cm.ns, cm.Js); r1, v1 = elements_state(apo_el, ctx.mu)
                F = CowellForces(grav=V.grav, sun=V.sun, drag=V.drag)
                if mode == "chain"
                    ev = cowell_leg(cmv, F, r1, v1, t_apo, 0.0, 0.0, bodies, 2 * npass; start=:apo, reltol=1e-11)
                    p, _ = stp(ctxv, Sv, zf, zo, tf, apo_el, t_apo, et0, npass); perv = [e for e in ev if e.kind === :peri]
                    for j in 1:min(length(p.t_s), length(perv))
                        push!(rows, (; tag, k0=k, variant=vn, j, dt=p.t_s[j] - (perv[j].t - t_apo), drp=p.radius_m[j] - norm(perv[j].r), T_C=j > 1 ? perv[j].t - perv[j - 1].t : NaN))
                    end
                else
                    ev = cowell_leg(cmv, F, r1, v1, t_apo, 0.0, 0.0, bodies, 3; start=:apo)
                    p, _ = stp(ctxv, Sv, zf, zo, tf, apo_el, t_apo, et0, 2)
                    (length(p.t_s) >= 2 && length(ev) >= 3) || continue
                    tC1 = ev[1].t - t_apo; tC2 = ev[3].t - t_apo
                    push!(rows, (; tag, k, variant=vn, T_C=tC2 - tC1, E_half=p.t_s[1] - tC1, E_full=(p.t_s[2] - p.t_s[1]) - (tC2 - tC1), dr1_m=p.radius_m[1] - norm(ev[1].r), dr2_m=p.radius_m[2] - norm(ev[3].r),
                        valid_dt_s=valid_dt, valid_dr_m=valid_dr))
                end
            end
        end
    end
    CSV.write(abspath(get(CLIE, "out", "diag_v5_earth.csv")), DataFrame(rows))
end
main_e()
