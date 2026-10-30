# R2 of estimator_v5_plan.md: the drag increment of one pass against the Cowell reference (DIAGNOSTIC). For 3 Mars development cells and the truth passes k = 1, 5, 10, 20, 30, 40, the Cowell
# propagator (J2 + drag, same fit) is run from the Cowell apoapsis state after truth pass k with and without drag; the reference is the change of the osculating a at the next apoapsis
# (phase-aligned, so the J2 short period cancels). The stepper's total change of the MEAN a through `shift_mean` (the conversion the stepper uses) is compared for the Kepler quadrature
# of v1 to v4 and for the integrated window at several resolutions. Output: relative errors (percent) per state and the summary.
#   julia --project=. diag_v5_window_check.jl --root=R [--config=MARS_TOML] [--out=FILE]
using Statistics
include(joinpath(@__DIR__, "mars_common.jl")); include(joinpath(@__DIR__, "estimator_v5.jl")); include(joinpath(@__DIR__, "diag_v5_lib.jl"))
const CLIW = Dict(m.captures[1] => m.captures[2] for a in ARGS for m in [match(r"^--([^=]+)=(.*)$", a)] if m !== nothing)
const ROOTW = abspath(CLIW["root"])
mc = load_mars(abspath(get(CLIW, "config", joinpath(@__DIR__, "mars_generality_config_g20.toml")))); ms = mars_setup(mc); out = mars_results_dir(ROOTW, mc)
S = settings_from(TOML.parsefile(joinpath(@__DIR__, "estimator_v2_config.toml")))
zf = zonal_field(ms.harm); tf = tesseral_field(ms.harm); cells = CSV.read(joinpath(out, "cells.csv"), DataFrame)
ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i])
combos = [(24, 6, 1.0), (32, 8, 1.0), (48, 8, 1.0), (96, 12, 1.0), (32, 8, 2.0), (48, 8, 4.0)]
rows = NamedTuple[]
for ci in (1, 40, 90)
    c = cells[ci, :]; ei = Int(c.epoch_idx); ctx = make_context(mars_record(ms, ep(ei)), S); cm = CowellModel(ctx, zf, tf)
    tr = CSV.read(joinpath(out, "truth", c.tag, "passes_3b.csv"), DataFrame); rho = fit_rho(ctx)
    for k in (1, 5, 10, 20, 30, 40)
        k + 1 <= nrow(tr) || continue
        el = (a=tr.a_m[k], e=tr.e[k], i=tr.i_rad[k], raan=tr.raan_rad[k], omega=tr.omega_rad[k], nu=0.0); r0, v0 = elements_state(el, ctx.mu)
        ev0 = cowell_leg(cm, CowellForces(sun=false, drag=true), r0, v0, tr.time_s[k], 0.0, 0.0, SV3(1.0, 0, 0), 1; start=:peri)
        A = state_elements(ev0[1].r, ev0[1].v, ctx.mu); ra, va = elements_state(A, ctx.mu)
        evn = cowell_leg(cm, CowellForces(grav=:j2, sun=false, drag=false), ra, va, 0.0, 0.0, 0.0, SV3(1.0, 0, 0), 2; start=:apo)
        evd = cowell_leg(cm, CowellForces(grav=:j2, sun=false, drag=true), ra, va, 0.0, 0.0, 0.0, SV3(1.0, 0, 0), 2; start=:apo)
        Pn = state_elements(evn[1].r, evn[1].v, ctx.mu); An = state_elements(evn[2].r, evn[2].v, ctx.mu); Ad = state_elements(evd[2].r, evd[2].v, ctx.mu)
        ref = Ad.a - An.a; osc0 = (a=Pn.a, e=Pn.e, i=Pn.i, raan=Pn.raan, omega=Pn.omega)
        m_pre = mean_state(osculating_to_mean((a=Pn.a, e=Pn.e, i=Pn.i, raan=Pn.raan, omega=Pn.omega, nu=0.0); mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements, ctx); ov = osc_vec(m_pre, ctx)
        din, dout = drag_deltas(osc0, ctx, S, rho); mq = shift_mean(shift_mean(m_pre, rp_form(ov, din), ctx), rp_form(ov, dout), ctx)
        push!(rows, (; tag=c.tag, k, method="kepler_quadrature", err_pct=100 * ((mq[1] - m_pre[1]) / ref - 1), ref_da_m=ref))
        for cb in combos
            fi, fo = drag_window(osc0, ctx, S, rho, TimingOpts(n_window=cb[1], n_back=cb[2], cluster=cb[3]))
            mq = shift_mean(shift_mean(m_pre, fi, ctx), fo, ctx)
            push!(rows, (; tag=c.tag, k, method="window_$(cb[1])_$(cb[2])_$(cb[3])", err_pct=100 * ((mq[1] - m_pre[1]) / ref - 1), ref_da_m=ref))
        end
    end
end
D = DataFrame(rows); CSV.write(abspath(get(CLIW, "out", "diag_v5_window_check.csv")), D)
for (m, g) in pairs(groupby(D, :method))
    println(rpad(m.method, 24), " states ", nrow(g), "  mean ", round(mean(g.err_pct); digits=4), " %  median |err| ", round(median(abs.(g.err_pct)); digits=4), " %  max |err| ", round(maximum(abs.(g.err_pct)); digits=4), " %")
end
ctx1 = make_context(mars_record(ms, ep(1)), S); rho1 = fit_rho(ctx1); osc1 = (a=1.6e7, e=0.78, i=0.44, raan=0.0, omega=0.1); drag_window(osc1, ctx1, S, rho1, TimingOpts())
t0 = time(); for _ in 1:200; drag_window(osc1, ctx1, S, rho1, TimingOpts()); end
println("cost of one drag_window call at the defaults: ", round((time() - t0) / 200 * 1e3; digits=3), " ms (single call, warm)")
