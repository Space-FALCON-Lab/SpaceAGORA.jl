# Diagnostic library for the pass-timing study (estimator v5): a self-contained Cowell propagator of ONE orbit leg with switchable forces, used to
# separate the pass-timing error budget of the stepper (v4) into drag quadrature, mean-motion/period formula, osculating-to-mean map and their interplay.
# This is a DIAGNOSTIC (never enters a verdict). Forces mirror the truth run: point-mass + zonal 2..L + non-zonal field (rotating planet) + Sun tide + ballistic
# drag with the piecewise/exponential fit at geodetic altitude, zero wind, relative velocity v - w x r. Frame: the body-equator inertial frame B of the estimator.
# No src/ or examples/ edit. Included by diag_v5_timing.jl.
isdefined(@__MODULE__, :predict_history4) || include(joinpath(@__DIR__, "estimator_v4.jl"))
using OrdinaryDiffEq

"""Forces of a Cowell leg: `grav` in (:point, :j2, :zonal, :full) (point mass; + J2; + zonal 2..L; + non-zonal field), `sun`, `drag`."""
Base.@kwdef struct CowellForces
    grav::Symbol = :full
    sun::Bool = true
    drag::Bool = true
end

"""Everything a Cowell leg needs besides the state: context (mu, planet, drag fit, beta), zonal and tesseral fields."""
struct CowellModel
    ctx::EstContext
    zf::ZonalField
    tf::TesseralField
    ns::Vector{Int}
    Js::Vector{Float64}
end
function CowellModel(ctx::EstContext, zf::ZonalField, tf::TesseralField)
    ns = collect(2:zf.L); return CowellModel(ctx, zf, tf, ns, [zf.J[n] for n in ns])
end

"""Integrate from the state (r0, v0) at absolute time t0 (s from the run epoch), starting AT an apse (`start` = :peri or :apo), until `nev` radial-velocity
sign changes after the start have been recorded. Rotation angle theta(t) = th0 + w (t - t0) (uniform over a leg). Sun held at its position at t0 + `sun_dt`.
Returns the events `(kind, t, r, v)`."""
function cowell_leg(cm::CowellModel, F::CowellForces, r0::SV3, v0::SV3, t0::Float64, th0::Float64, w::Float64, rb_sun, nev::Int; start::Symbol=:peri,
                    reltol=1e-12, abstol=1e-7, saveat=nothing)
    ctx = cm.ctx; mu = ctx.mu; kw = kick_work(cm.tf.L); Re = ctx.Re
    rbf = rb_sun isa SV3 ? (t -> rb_sun) : rb_sun       # a fixed position or a function of the absolute time
    Lt = cm.tf.L; ns = cm.ns; Js = cm.Js; nzmax = F.grav === :j2 ? [2] : ns
    Jmax = F.grav === :j2 ? [cm.zf.J[2]] : Js
    function f(u, p, t)
        x, y, z = u[1], u[2], u[3]; r = SV3(x, y, z); vv = SV3(u[4], u[5], u[6]); rn = sqrt(x * x + y * y + z * z)
        a = -mu / rn^3 * r
        if F.grav !== :point; a += zonal_acc(r, nzmax, Jmax, mu, Re); end
        if F.grav === :full
            th = th0 + w * (t - t0); c, s = cos(th), sin(th)
            abx, aby, abz = nz_acc(cm.tf, kw, Lt, Lt, c * x + s * y, -s * x + c * y, z)
            a += SV3(c * abx - s * aby, s * abx + c * aby, abz)
        end
        if F.sun
            bs = rbf(t); bs isa SV3 && (bs = (bs,))
            for (b, mub) in zip(bs, ctx.body_mu); a += tide_acc(r, b, mub); end
        end
        if F.drag
            alt = geodetic_altitude_latitude(r, ctx.planet)[1]; rho = corridor_density(ctx.fit, alt)
            vrel = vv - cross(ctx.omega_e, r); V = norm(vrel)
            a += (-0.5 * rho * V / ctx.beta) * vrel
        end
        return SVector(u[4], u[5], u[6], a[1], a[2], a[3])
    end
    T0 = 2π * sqrt((norm(r0) / (2 - norm(r0) * dot(v0, v0) / mu))^3 / mu)
    ev = NamedTuple[]; sgn0 = start === :peri ? 1.0 : -1.0
    cond(u, t, integ) = u[1] * u[4] + u[2] * u[5] + u[3] * u[6]
    function hit!(integ, kind)
        integ.t - t0 > 0.02 * T0 || return
        push!(ev, (kind=kind, t=integ.t, r=SV3(integ.u[1], integ.u[2], integ.u[3]), v=SV3(integ.u[4], integ.u[5], integ.u[6])))
        length(ev) >= nev && terminate!(integ)
    end
    cb = ContinuousCallback(cond, integ -> hit!(integ, :peri), integ -> hit!(integ, :apo); save_positions=(false, false), abstol=1e-9, reltol=0.0)
    u0 = SVector(r0..., v0...)
    sol = saveat === nothing ? solve(ODEProblem(f, u0, (t0, t0 + 1.2 * T0 * (nev + 1) / 2 + 3600.0)), Vern9(); abstol=abstol, reltol=reltol, callback=cb, save_everystep=false) :
        solve(ODEProblem(f, u0, (t0, t0 + 1.2 * T0 * (nev + 1) / 2 + 3600.0)), Vern9(); abstol=abstol, reltol=reltol, callback=cb, save_everystep=false, saveat=saveat)
    saveat === nothing || return ev, [(t=sol.t[j], r=SV3(sol.u[j][1], sol.u[j][2], sol.u[j][3]), v=SV3(sol.u[j][4], sol.u[j][5], sol.u[j][6])) for j in eachindex(sol.t)]
    return ev
end

"""Osculating elements `(a, e, i, Omega, omega, nu)` of a Cartesian state (B frame)."""
function state_elements(r::SV3, v::SV3, mu)
    rn = norm(r); a = 1 / (2 / rn - dot(v, v) / mu); hv = cross(r, v); hn = norm(hv); W = hv / hn
    ev = cross(v, hv) / mu - r / rn; e = norm(ev)
    i = acos(clamp(W[3], -1.0, 1.0)); Om = atan(W[1], -W[2]); node = SV3(cos(Om), sin(Om), 0.0)
    om = atan(dot(ev, cross(W, node)), dot(ev, node)); nu = atan(dot(cross(ev, r), W) / (e * rn), dot(ev, r) / (e * rn))
    return (a=a, e=e, i=i, raan=mod(Om, 2π), omega=mod(om, 2π), nu=mod(nu, 2π))
end

"""Cartesian state (B frame) of osculating elements."""
function elements_state(el, mu)
    r, v = elements_to_cartesian(el; mu=mu); return SV3(r...), SV3(v...)
end

"""The context with another J2 (J2 = 0 gives the two-body limit of the mean-element machinery)."""
with_J2(c::EstContext, J2) = EstContext(c.mu, c.Re, J2, c.Rp_e, c.planet, c.omega_e, c.fit, c.beta, c.q_low, c.q_high, c.q_hard, c.w2m, c.apo_alt_m, c.body_names, c.body_mu, c.gl_x, c.gl_w, c.center, c.rot)

"""Third-body position(s) (B frame; one SV3 when the context has a single body, else a vector) as a function of the absolute time `t` (s from the run epoch): linear interpolation of the SPICE positions every `dt` seconds over [t0, t1]."""
function sun_track(ctx::EstContext, et0::Float64, t0::Float64, t1::Float64; dt::Float64=21600.0)
    ts = collect(t0:dt:(t1 + dt)); P = [body_positions(ctx, et0 + t) for t in ts]
    one = length(ctx.body_mu) == 1
    return t -> (j = clamp(floor(Int, (t - t0) / dt) + 1, 1, length(ts) - 1); u = (t - ts[j]) / dt;
                 one ? (1 - u) * P[j][1] + u * P[j + 1][1] : [(1 - u) * P[j][b] + u * P[j + 1][b] for b in eachindex(P[j])])
end
