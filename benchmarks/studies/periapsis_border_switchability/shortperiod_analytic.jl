# Closed-form (Pi-form) predictor of the periapsis-radius motion from one or two perturbers, and of its oscillation amplitude, period and phase.
#
# Everything here is derived from the orbit-averaged quadrupole and octupole disturbing functions (Lagrange's planetary equation for e) and
# the first-order short-period quadrature of estimator_v2.jl; no coefficient is fitted. Dimensionless groups:
#   Pi_3B = (mu_b / mu)(a / r_b)^3     (quadrupole strength),     a / r_b  (octupole relative strength),     nu = lambda_dot / n  (perturber-to-spacecraft rate)
# Slow (single-averaged) rate:   d r_p / dt = n a Pi_3B [ D2(e, bx, by) + (a / r_b) D3(e, bx, by) ],      `rp_rate_nd` in estimator_v2.jl
# with (bx, by) the direction cosines of the perturber along the apse line and the transverse in-plane direction. For a perturber on a circular orbit
# in a plane at inclination i_p to the orbit, bx, by are linear in cos(lambda), sin(lambda), so the rate has harmonics 0, 1, 2, 3 of the perturber
# longitude lambda (measured from the node of the orbit on the perturber plane); integrating the harmonic k gives an oscillation of amplitude
#     A_k = n a Pi_3B f_k R_k / (k |lambda_dot|),   R_k = sqrt(c_k^2 + s_k^2),  f_2 = 1 (quadrupole), f_1 = f_3 = a / r_b (octupole).
include(joinpath(@__DIR__, "estimator_v2.jl"))

"""Longitude `lam` of the perturber from the node line of the orbit on the perturber plane, the orbit's argument of periapsis `omega_p` and inclination `i_p` in that plane."""
function pert_frame(rb::SV3, vb::SV3, W::SV3, P::SV3)
    n = normalize(cross(rb, vb)); k = normalize(cross(n, W))
    lam = atan(dot(rb, cross(n, k)), dot(rb, k))
    return lam, atan(dot(W, cross(k, P)), dot(k, P)), acos(clamp(dot(W, n), -1.0, 1.0))
end

"""Harmonic coefficients (c_k, s_k) of  d r_p / dt / (n a Pi_3B)  in cos(k lambda), sin(k lambda): k = 0, 2 quadrupole, k = 1, 3 octupole (to be multiplied by a / r_b)."""
function harmonic_coeffs(e, i, w)
    eta = sqrt(1 - e^2); si2 = sin(i)^2; sw = sin(w); cw = cos(w); ci = cos(i)
    c0 = -15 * e * eta * si2 * sw * cw / 4
    c2 = 15 * e * eta * (si2 - 2) * sw * cw / 4;   s2 = 15 * e * eta * ci * cos(2w) / 4
    c1 = -15 * eta * (105 * e^2 * si2 * sw^2 - 75 * e^2 * si2 - 3e^2 + 5si2 - 4) * sw / 64
    s1 = 15 * eta * (105 * e^2 * si2 * sw^2 - 15 * e^2 * si2 - 3e^2 + 15si2 - 4) * ci * cw / 64
    c3 = 75 * eta * (21 * e^2 * si2 * sw^2 - 15 * e^2 * si2 - 28 * e^2 * sw^2 + 21e^2 + si2) * sw / 64
    s3 = -75 * eta * (7 * e^2 * si2 * sw^2 - e^2 * si2 - 28 * e^2 * sw^2 + 7e^2 + si2) * ci * cw / 64
    return (c0=c0, c1=c1, s1=s1, c2=c2, s2=s2, c3=c3, s3=s3)
end

"""Rate of r_p / (n a Pi) synthesized from the harmonics (consistency check of `harmonic_coeffs` against `rp_rate_nd`)."""
function rate_from_harmonics(e, i, w, lam, a_over_b)
    h = harmonic_coeffs(e, i, w)
    return h.c0 + h.c2 * cos(2lam) + h.s2 * sin(2lam) + a_over_b * (h.c1 * cos(lam) + h.s1 * sin(lam) + h.c3 * cos(3lam) + h.s3 * sin(3lam))
end

bodystate(b::String, et; center::String="EARTH", rot=nothing) = (st = spkezr(uppercase(b), et, "J2000", "NONE", center)[1];
    rot === nothing ? (SV3(st[1:3]) .* 1e3, SV3(st[4:6]) .* 1e3) : (SV3(rot' * (SV3(st[1:3]) .* 1e3)), SV3(rot' * (SV3(st[4:6]) .* 1e3))))

"""
Closed-form mean periapsis-radius change and short-period term at the pass times `tp` (s from `et0`), starting from the mean elements `m0 = (a, e, i, Omega, omega)`.
Only J2 secular rates of Omega and omega move the orbit (a, e fixed: no secular a, small e change), the bodies are at their SPICE positions.
Returns `(rbar, sp, rate_t, ...)`: `rbar[k]` mean r_p change from the first pass, `sp[k]` short-period r_p offset (periapsis) per pass, both summed over `bodies`.
"""
function analytic_series(ctx::EstContext, m0::SV5, et0, tp::Vector{Float64}, bodies::Vector{String}; dt_s::Float64=0.05 * SECONDS_PER_DAY, mu_of=(b -> MU_BODY[b]))
    a, k0, h0, inc, Om0 = m0; e = hypot(k0, h0); w0 = atan(h0, k0)
    F0 = OrbitFrame(a, k0, h0, inc, Om0, ctx.mu); dOm, dw, _ = j2_rates(F0, ctx.mu, ctx.J2, ctx.Re); n = F0.n
    tgrid = collect(tp[1]:dt_s:tp[end] + dt_s)
    rate = zeros(length(tgrid))
    spk = zeros(length(tp))
    function frame_at(t)
        w = w0 + dw * t; Om = Om0 + dOm * t
        return OrbitFrame(a, e * cos(w), e * sin(w), inc, Om, ctx.mu)
    end
    for b in bodies
        mub = mu_of(b)
        for (j, t) in enumerate(tgrid)
            rb, _ = bodystate(b, et0 + t; center=ctx.center, rot=ctx.rot); F = frame_at(t); rn = norm(rb); u = rb / rn
            Pi = mub / ctx.mu * (a / rn)^3
            rate[j] += n * a * Pi * rp_rate_nd(e, dot(u, F.P), dot(u, F.Q), a / rn)
        end
        for (kk, t) in enumerate(tp)
            rb, _ = bodystate(b, et0 + t; center=ctx.center, rot=ctx.rot); F = frame_at(t); rn = norm(rb); u = rb / rn
            spk[kk] += a * mub / ctx.mu * (a / rn)^3 * quad_drp_coeff(e, dot(u, F.P), dot(u, F.Q))
        end
    end
    cum = [0.0; cumsum(0.5 .* (rate[1:end-1] .+ rate[2:end]) .* diff(tgrid))]
    rbar = [interp1(tgrid, cum, t) for t in tp]
    return (rbar=rbar .- rbar[1], sp=spk, tgrid=tgrid, rate=rate)
end

interp1(x, y, xv) = (j = clamp(searchsortedlast(x, xv), 1, length(x) - 1); y[j] + (y[j+1] - y[j]) * (xv - x[j]) / (x[j+1] - x[j]))

"""
Closed-form oscillation parameters of one perturber over the window `[t_a, t_b]` (s): the harmonics are evaluated at the window midpoint (e, omega', i' there), the perturber
longitude rate `lam_dot` by a linear fit of the unwrapped longitude over the window (it includes the nodal motion of the orbit on the perturber plane).
"""
function oscillation_params(ctx::EstContext, m0::SV5, et0, t_a, t_b, body::String; npts::Int=121)
    a, k0, h0, inc, Om0 = m0; e = hypot(k0, h0); w0 = atan(h0, k0); mub = MU_BODY[body]
    F0 = OrbitFrame(a, k0, h0, inc, Om0, ctx.mu); dOm, dw, _ = j2_rates(F0, ctx.mu, ctx.J2, ctx.Re); n = F0.n
    ts = collect(range(t_a, t_b, length=npts)); lam = zeros(npts); omp = zeros(npts); ip = zeros(npts); rbs = zeros(npts)
    for (j, t) in enumerate(ts)
        w = w0 + dw * t; F = OrbitFrame(a, e * cos(w), e * sin(w), inc, Om0 + dOm * t, ctx.mu)
        rb, vb = bodystate(body, et0 + t; center=ctx.center, rot=ctx.rot); lam[j], omp[j], ip[j] = pert_frame(rb, vb, F.W, F.P); rbs[j] = norm(rb)
    end
    lam = unwrap(lam)
    X = [ones(npts) ts .- ts[1]]; coef = X \ lam; lamdot = coef[2]
    jm = (npts + 1) ÷ 2; tm = ts[jm]
    h = harmonic_coeffs(e, ip[jm], omp[jm]); rbm = mean(rbs); Pi = mub / ctx.mu * (a / rbm)^3; aob = a / rbm
    lam_m = coef[1] + lamdot * (tm - ts[1])
    out = Dict{Symbol, Any}(:a => a, :e => e, :Pi => Pi, :aob => aob, :nu => lamdot / n, :lamdot => lamdot, :omega_p => omp[jm], :i_p => ip[jm], :lam_mid => lam_m,
        :t_mid => tm, :n => n, :c0 => h.c0)
    hs = [harmonic_coeffs(e, ip[j], omp[j]) for j in 1:npts]
    for (k, f) in ((1, aob), (2, 1.0), (3, aob))
        c = getfield(h, Symbol(:c, k)); s_ = getfield(h, Symbol(:s, k)); R = hypot(c, s_)
        # phase of the harmonic, psi_k(t) = atan(s_k, c_k), drifts as omega' and i' change: oscillation  R cos(k lam - psi_k)  has frequency k lam_dot - psi_dot
        psis = unwrap([atan(getfield(x, Symbol(:s, k)), getfield(x, Symbol(:c, k))) for x in hs]); pc = X \ psis; psidot = pc[2]; psi_m = pc[1] + psidot * (tm - ts[1])
        for (suffix, psi, weff) in (("", psi_m, k * lamdot - psidot), ("f", atan(s_, c), k * lamdot))     # "f": phase frozen at the window midpoint
            amp = n * a * Pi * f * R / weff                                    # signed: osc(t) = amp * sin(k lam(t) - psi(t))
            out[Symbol(:A, k, suffix)] = abs(amp); out[Symbol(:P, k, suffix)] = 2π / abs(weff); out[Symbol(:psi, k, suffix)] = psi
            target = amp >= 0 ? π / 2 : -π / 2
            phase_now = mod(k * lam_m - psi - target + π, 2π) - π
            out[Symbol(:tmax, k, suffix)] = tm - phase_now / weff
        end
        out[Symbol(:c, k)] = c; out[Symbol(:s, k)] = s_; out[Symbol(:R, k)] = R; out[Symbol(:psidot, k)] = psidot
    end
    return out
end

function unwrap(x)
    y = copy(x)
    for j in 2:length(y); y[j] = y[j-1] + (mod(x[j] - x[j-1] + π, 2π) - π); end
    return y
end
