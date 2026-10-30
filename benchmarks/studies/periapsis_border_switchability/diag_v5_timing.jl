# Pass-timing error budget of the stepper by truth-state injection (DIAGNOSTIC for estimator v5; never enters a verdict).
# For a truth cell and a pass k: the truth periapsis state of pass k is propagated by the Cowell emulator (diag_v5_lib.jl, validated against the engine truth to ms)
# to the apoapsis A_k; from A_k the stepper (v4, or v5 variants) predicts two passes (k+1, k+2) and the Cowell emulator, with the SAME force subset, produces the reference.
# Errors: E_half = t_hat(k+1) - t(k+1) (apoapsis -> periapsis), E_full = (t_hat(k+2) - t_hat(k+1)) - (t(k+2) - t(k+1)) (one fresh periapsis-to-periapsis period).
# The force subsets (variants) separate: gravity/period model (no drag), drag quadrature (two-body + drag), J2 map interplay (J2 + drag), everything.
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. diag_v5_timing.jl --planet=mars|earth --root=DIR [--cells=tag,tag,..] [--ks=1,2,3,..] [--variants=a,b,..] [--stepper=v4|v5] [--out=FILE]
using Statistics
include(joinpath(@__DIR__, "mars_common.jl"))
include(joinpath(@__DIR__, "diag_v5_lib.jl"))
isfile(joinpath(@__DIR__, "estimator_v5.jl")) && include(joinpath(@__DIR__, "estimator_v5.jl"))

function parse_cli_d()
    kv = Dict{String, String}()
    for a in ARGS
        m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2])
    end
    return kv
end
const CLID = parse_cli_d()


"""Variants: name => (J2 on, sun, drag, zonal rule (:rule or :none), tesseral (:rule, :off), Cowell gravity)."""
const VARIANTS = Dict(
    "full"         => (j2=true,  sun=true,  drag=true,  zonal=:rule, tess=:rule, grav=:full),
    "nodrag"       => (j2=true,  sun=true,  drag=false, zonal=:rule, tess=:rule, grav=:full),
    "grav_sun_nodrag_notess" => (j2=true, sun=true, drag=false, zonal=:rule, tess=:off, grav=:zonal),
    "j2only"       => (j2=true,  sun=false, drag=false, zonal=:none, tess=:off,  grav=:j2),
    "zonalonly"    => (j2=true,  sun=false, drag=false, zonal=:rule, tess=:off,  grav=:zonal),
    "j2drag"       => (j2=true,  sun=false, drag=true,  zonal=:none, tess=:off,  grav=:j2),
    "drag2b"       => (j2=false, sun=false, drag=true,  zonal=:none, tess=:off,  grav=:point),
    "point"        => (j2=false, sun=false, drag=false, zonal=:none, tess=:off,  grav=:point),
    "nodrag_notess"=> (j2=true,  sun=true,  drag=false, zonal=:rule, tess=:off,  grav=:zonal),
    "zonaldrag"    => (j2=true,  sun=false, drag=true,  zonal=:rule, tess=:off,  grav=:zonal),
    "j2sundrag"    => (j2=true,  sun=true,  drag=true,  zonal=:none, tess=:off,  grav=:j2),
    "j2sun_nodrag" => (j2=true,  sun=true,  drag=false, zonal=:none, tess=:off,  grav=:j2),
    "tessdrag"     => (j2=true,  sun=false, drag=true,  zonal=:rule, tess=:rule, grav=:full),
    "drag_notess"  => (j2=true,  sun=true,  drag=true,  zonal=:rule, tess=:off,  grav=:zonal),
)

function stepper_fn(name)
    name == "v4" && return predict_history4
    name == "v5" && return (a...; kw...) -> predict_history5(a[1:8]..., TimingOpts(), a[9:end]...; kw...)
    name == "v5kep" && return (a...; kw...) -> predict_history5(a[1:8]..., TimingOpts(window=:kepler), a[9:end]...; kw...)   # window = v4 drag, everything else v5
    name == "v5zsp" && return (a...; kw...) -> predict_history5(a[1:8]..., TimingOpts(zonal_sp=true), a[9:end]...; kw...)   # v5 plus the zonal short period (tested and rejected: development note)
    name == "v5nokp" && return (a...; kw...) -> predict_history5(a[1:8]..., TimingOpts(kick_post_drag=false), a[9:end]...; kw...)   # integrated window, v4 kick (pre-drag outbound half)
    throw(ArgumentError("unknown stepper $name"))
end

function run_variant(stepper, V, cm::CowellModel, ctx0, S0, O, TAB, zf, tf, apo_el, t_apo, et0, th_apo, w, rbsun, extra...)
    ctx = V.j2 ? ctx0 : with_J2(ctx0, 0.0)
    S = with_settings(S0; tides=V.sun, drag=V.drag, j2=V.j2)
    zo = V.zonal === :none ? ZonalOpts(force=Int[]) : ZonalOpts()
    tt = TesseralOpts(force=true)      # always the carried (timing-corrected) path; a 2-pass history cannot evaluate the carry rule
    tfv = V.tess === :off ? zero_field(tf) : tf
    cm2 = CowellModel(ctx, zf, tf, cm.ns, cm.Js)
    r0, v0 = elements_state(apo_el, ctx.mu)
    F = CowellForces(grav=V.grav, sun=V.sun, drag=V.drag)
    ev = cowell_leg(cm2, F, r0, v0, t_apo, th_apo, w, rbsun, 3; start=:apo)       # peri(k+1), apo(k+1), peri(k+2)
    p, info = stepper(ctx, S, O, TAB, zf, zo, tfv, tt, NaN, 0.0, 0.0, 0.0, et0 + t_apo; max_passes=2, init_override=apo_el, extra...)
    return ev, p, info
end

function diag_cell(stepper, ms, mc, tag, ks, varnames, S, O, TAB, out, ei_of)
    ei = ei_of; zf = zonal_field(ms.harm); tf = tesseral_field(ms.harm)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); etm(i) = SimulationModel.ephemerides_time_seconds(ep(i), SpiceEphemeridesModel())
    ctx = make_context(mars_record(ms, ep(ei)), S); cm = CowellModel(ctx, zf, tf); et0 = etm(ei)
    tr = CSV.read(joinpath(out, "truth", tag, "passes_3b.csv"), DataFrame); rotm = RotationModel(ctx, et0)
    rows = NamedTuple[]
    for k in ks
        k + 2 <= nrow(tr) || continue
        el = (a=tr.a_m[k], e=tr.e[k], i=tr.i_rad[k], raan=tr.raan_rad[k], omega=tr.omega_rad[k], nu=0.0)
        r0, v0 = elements_state(el, ctx.mu); t0 = tr.time_s[k]; th, w = rotation_at(rotm, t0)
        rbsun = body_positions(ctx, et0 + t0 + 0.5 * (tr.time_s[k + 1] - t0))[1]
        ev0 = cowell_leg(cm, CowellForces(), r0, v0, t0, th, w, rbsun, 1; start=:peri)
        t_apo = ev0[1].t; apo_el = state_elements(ev0[1].r, ev0[1].v, ctx.mu)
        th_a = th + w * (t_apo - t0)
        for vn in varnames
            V = VARIANTS[vn]
            ev, p, info = run_variant(stepper, V, cm, ctx, S, O, TAB, zf, tf, apo_el, t_apo, et0, th_a, w, rbsun)
            length(p.t_s) >= 2 && length(ev) >= 3 || continue
            e1 = state_elements(ev[1].r, ev[1].v, ctx.mu); e3 = state_elements(ev[3].r, ev[3].v, ctx.mu)
            tC1 = ev[1].t - t_apo; tC2 = ev[3].t - t_apo
            push!(rows, (; tag, k, variant=vn, T_C=tC2 - tC1, E_half=p.t_s[1] - tC1, E_full=(p.t_s[2] - p.t_s[1]) - (tC2 - tC1),
                dr1_m=p.radius_m[1] - norm(ev[1].r), dr2_m=p.radius_m[2] - norm(ev[3].r), da1_m=p.a_m[1] - e1.a, da2_m=p.a_m[2] - e3.a,
                dec_a=(e1.a - apo_el.a)))
        end
    end
    return rows
end


"""Multi-pass consistency by force subset: from the Cowell apoapsis state after the truth pass `k0`, both the Cowell chain (variant forces) and the stepper (same subset) run `npass` passes;
returns the pass-time and radius errors of the stepper pass by pass (accumulated error budget by physics subset)."""
function chain_cell(stepper, ms, mc, tag, k0s, varnames, S, O, TAB, out, ei, npass)
    zf = zonal_field(ms.harm); tf = tesseral_field(ms.harm)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); etm(i) = SimulationModel.ephemerides_time_seconds(ep(i), SpiceEphemeridesModel())
    ctx = make_context(mars_record(ms, ep(ei)), S); cm = CowellModel(ctx, zf, tf); et0 = etm(ei)
    tr = CSV.read(joinpath(out, "truth", tag, "passes_3b.csv"), DataFrame); rotm = RotationModel(ctx, et0)
    rows = NamedTuple[]
    for k0 in k0s
        k0 + 2 <= nrow(tr) || continue
        el = (a=tr.a_m[k0], e=tr.e[k0], i=tr.i_rad[k0], raan=tr.raan_rad[k0], omega=tr.omega_rad[k0], nu=0.0)
        r0, v0 = elements_state(el, ctx.mu); t0 = tr.time_s[k0]; th, w = rotation_at(rotm, t0)
        rbsun = sun_track(ctx, et0, t0, t0 + (npass + 2) * 70000.0)
        ev0 = cowell_leg(cm, CowellForces(), r0, v0, t0, th, w, rbsun, 1; start=:peri)
        t_apo = ev0[1].t; apo_el = state_elements(ev0[1].r, ev0[1].v, ctx.mu); th_a = th + w * (t_apo - t0)
        for vn in varnames
            V = VARIANTS[vn]
            ctxv = V.j2 ? ctx : with_J2(ctx, 0.0); Sv = with_settings(S; tides=V.sun, drag=V.drag, j2=V.j2)
            zo = V.zonal === :none ? ZonalOpts(force=Int[]) : ZonalOpts(); tt = TesseralOpts(force=true); tfv = V.tess === :off ? zero_field(tf) : tf
            cmv = CowellModel(ctxv, zf, tf, cm.ns, cm.Js); r1, v1 = elements_state(apo_el, ctx.mu)
            ev = cowell_leg(cmv, CowellForces(grav=V.grav, sun=V.sun, drag=V.drag), r1, v1, t_apo, th_a, w, rbsun, 2 * npass; start=:apo, reltol=1e-11)
            perv = [e for e in ev if e.kind === :peri]
            tth = haskey(CLID, "oracle") ? [e.t - t_apo for e in perv] : nothing        # UNFROZEN diagnostic: the rotation angle of every kick at the REFERENCE pass time
            p, info = stepper(ctxv, Sv, O, TAB, zf, zo, tfv, tt, NaN, 0.0, 0.0, 0.0, et0 + t_apo; max_passes=npass, init_override=apo_el, theta_times=tth)
            for j in 1:min(length(p.t_s), length(perv))
                push!(rows, (; tag, k0, variant=vn, j, dt=p.t_s[j] - (perv[j].t - t_apo), drp=p.radius_m[j] - norm(perv[j].r), T_C=j > 1 ? perv[j].t - perv[j - 1].t : NaN))
            end
        end
    end
    return rows
end

function main_d()
    planet = get(CLID, "planet", "mars"); root = abspath(CLID["root"])
    mc = load_mars(abspath(get(CLID, "config", joinpath(@__DIR__, "mars_generality_config_g20.toml")))); ms = mars_setup(mc); out = joinpath(root, String(mc["meta"]["subdir"]))
    cfg = TOML.parsefile(joinpath(@__DIR__, "estimator_v2_config.toml")); v = cfg["v2"]
    O = V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])); TAB = SPTable(O.n_sp)
    S = settings_from(cfg); nf = joinpath(out, "estimator_numerics.toml"); isfile(nf) && (S = with_settings(S; n_drag_half=Int(TOML.parsefile(nf)["n_drag_half"])))
    cells = CSV.read(joinpath(out, get(CLID, "cellfile", "cells.csv")), DataFrame)
    tags = haskey(CLID, "cells") ? split(CLID["cells"], ",") : String.(cells.tag)
    ks = haskey(CLID, "ks") ? parse.(Int, split(CLID["ks"], ",")) : [1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 60, 80]
    vns = haskey(CLID, "variants") ? String.(split(CLID["variants"], ",")) : ["full", "nodrag", "j2only", "j2drag", "drag2b"]
    stepper = stepper_fn(get(CLID, "stepper", "v4"))
    rows = NamedTuple[]
    for tag in tags
        ei = Int(cells[cells.tag .== tag, :].epoch_idx[1])
        t0 = time()
        if get(CLID, "mode", "inject") == "chain"
            k0s = haskey(CLID, "k0s") ? parse.(Int, split(CLID["k0s"], ",")) : [1]
            append!(rows, chain_cell(stepper, ms, mc, tag, k0s, vns, S, O, TAB, out, ei, parse(Int, get(CLID, "npass", "40"))))
        else
            append!(rows, diag_cell(stepper, ms, mc, tag, ks, vns, S, O, TAB, out, ei))
        end
        println(tag, " done ", round(time() - t0; digits=1), " s")
    end
    CSV.write(abspath(get(CLID, "out", "diag_v5_timing.csv")), DataFrame(rows))
end
abspath(PROGRAM_FILE) == @__FILE__() && main_d()
