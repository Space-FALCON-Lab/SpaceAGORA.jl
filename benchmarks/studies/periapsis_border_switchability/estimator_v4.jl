# Estimator v4: estimator v3 plus a per-pass TESSERAL/SECTORIAL kick (non-zonal gravity, orders m >= 1) inside the pass-to-pass stepper.
#
# Derivation and the closed forms it rests on: tasks/2026-10-15_periapsis-region-border-targeting/tesseral_derivation.md (sections 3 to 5);
# reference implementation (Python): tesseral_lib.py (`gauss_kick`, `term_kicks`), tesseral_propagator.py. Nothing is fitted.
#
# What is added: for every pass k, the first-order change of the osculating elements (a, e, i, Omega, omega) over the apoapsis-to-apoapsis Kepler arc of the
# current (mean) elements, by Gauss quadrature of the non-zonal acceleration of the planet's gravity file (degree 2..Lt, order 1..Mt, Cunningham recursion)
# with the Gauss variational equations, the planet turning by theta(t) = theta_p + w_rot tau (tau = time from periapsis). The arc is split at periapsis:
#   inbound half  (apoapsis -> periapsis): added to the mean state BEFORE the drag step of the pass, so that the physical periapsis radius is the osculating
#                  radius at the previous apoapsis plus the inbound half kick (the derivation's T5a "pass-1 offset"), and
#   outbound half (periapsis -> apoapsis): added after the pass.
# Timing: the kick also moves the times of the periapsis and of the next apoapsis (through the change of a along the arc and the Gauss mean-anomaly
# equation): dt = 1.5 (int Delta a dt)/a - (int dM dt)/N, carried into t_apo; the planet's rotation angle at each periapsis time comes from the engine's own
# frame routine (`SimulationModel.planet_frame_lpi`, SPICE IAU_<BODY>, pck00011 constants) through `RotationModel`.
#
# Carry rule (no planet-name branching): D_T > budget (1 km), D_T = the horizon rms of the tesseral periapsis-radius change from the parameter Pi_T of
# tesseral_param.py (sigma_n from the gravity file, unit kicks at the first-pass elements, accumulation gain from the v3 history's commensurability sweep).
# With the rule false the v3 prediction is returned UNCHANGED (same object, byte-identical output).
# Usage: included by run_estimator_v4.jl. No main().
isdefined(@__MODULE__, :predict_history3) || include(joinpath(@__DIR__, "estimator_v3.jl"))

# ---- tesseral data of a planet record --------------------------------------------------------------------------------------

"""Non-zonal gravity data of a harmonics model: unnormalized C, S (index [n+1, m+1], degree 2..L, order 1..n, zero elsewhere), normalized Cbar, Sbar, factors `Nf`
(unnormalized = Nf * normalized), `sigma[n+1]` (rms of the normalized non-zonal coefficients of degree n: the Pi_T input), GM, reference radius."""
struct TesseralField
    L::Int
    mu::Float64
    R::Float64
    C::Matrix{Float64}
    S::Matrix{Float64}
    Cbar::Matrix{Float64}
    Sbar::Matrix{Float64}
    Nf::Matrix{Float64}
    sigma::Vector{Float64}
end

"""Fully normalized -> unnormalized factor sqrt((2 - delta_m0)(2n+1)(n-m)!/(n+m)!), (n-m)!/(n+m)! as 1/prod_{k=n-m+1}^{n+m} k."""
norm_factor(n::Int, m::Int) = sqrt((m == 0 ? 1.0 : 2.0) * (2n + 1) / prod(float, (n - m + 1):(n + m); init=1.0))

function tesseral_field(harm; L::Int=Int(harm.L))
    L = min(L, Int(harm.L)); C = zeros(L + 1, L + 1); S = zeros(L + 1, L + 1); Cb = zeros(L + 1, L + 1); Sb = zeros(L + 1, L + 1); Nf = zeros(L + 1, L + 1)
    sig = zeros(L + 1)
    for n in 2:L
        for m in 1:n
            f = norm_factor(n, m); Nf[n + 1, m + 1] = f
            Cb[n + 1, m + 1] = harm.C[n + 1, m + 1]; Sb[n + 1, m + 1] = harm.S[n + 1, m + 1]
            C[n + 1, m + 1] = f * Cb[n + 1, m + 1]; S[n + 1, m + 1] = f * Sb[n + 1, m + 1]
        end
        sig[n + 1] = sqrt(sum(Cb[n + 1, m + 1]^2 + Sb[n + 1, m + 1]^2 for m in 1:n) / (2n))
    end
    return TesseralField(L, harm.gm_m3s2, harm.reference_radius_m, C, S, Cb, Sb, Nf, sig)
end

"""Same field with all coefficients zero (the no-tesseral limit of the stepper, selftest f)."""
zero_field(t::TesseralField) = TesseralField(t.L, t.mu, t.R, zero(t.C), zero(t.S), t.Cbar, t.Sbar, t.Nf, t.sigma)

mutable struct KickWork
    V::Matrix{Float64}
    W::Matrix{Float64}
end
KickWork(L::Int) = KickWork(zeros(L + 3, L + 3), zeros(L + 3, L + 3))
const _KICK_WORK = Dict{Int, KickWork}()
kick_work(L::Int) = get!(() -> KickWork(L), _KICK_WORK, L)

const _GL01 = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
gl01(n::Int) = get!(() -> gauss_legendre01(n), _GL01, n)

# ---- Cunningham recursion and the non-zonal acceleration ---------------------------------------------------------------------

"""Fill V, W (index [n+1, m+1]) with V_nm + i W_nm = (R/r)^(n+1) P_nm(sin phi) e^(i m lambda), n <= Lt + 1, m <= min(n, Mt + 1), at the body-fixed point (x, y, z)."""
@inline function cunningham!(kw::KickWork, R, x, y, z, Lt::Int, Mt::Int)
    V, W = kw.V, kw.W
    r2 = x * x + y * y + z * z; rho = R * R / r2; x0 = R * x / r2; y0 = R * y / r2; z0 = R * z / r2
    V[1, 1] = R / sqrt(r2); W[1, 1] = 0.0
    @inbounds for n in 1:Lt + 1
        V[n + 1, n + 1] = (2n - 1) * (x0 * V[n, n] - y0 * W[n, n]); W[n + 1, n + 1] = (2n - 1) * (x0 * W[n, n] + y0 * V[n, n])
        for m in 0:min(n - 1, Mt + 1)
            c1 = (2n - 1) / (n - m); v = c1 * z0 * V[n, m + 1]; w = c1 * z0 * W[n, m + 1]
            if n >= 2
                c2 = (n + m - 1) / (n - m) * rho; v -= c2 * V[n - 1, m + 1]; w -= c2 * W[n - 1, m + 1]
            end
            V[n + 1, m + 1] = v; W[n + 1, m + 1] = w
        end
    end
    return nothing
end

"""Body-fixed acceleration of the non-zonal terms (degree 2..Lt, order 1..min(n, Mt)) at (x, y, z) (Montenbruck and Gill 2000, sec. 3.2; same arrangement as
tesseral_lib.accel_body with the m = 0 column removed)."""
@inline function nz_acc(t::TesseralField, kw::KickWork, Lt::Int, Mt::Int, x, y, z)
    cunningham!(kw, t.R, x, y, z, Lt, Mt)
    V, W = kw.V, kw.W; ax = 0.0; ay = 0.0; az = 0.0
    @inbounds for n in 2:Lt
        for m in 1:min(n, Mt)
            C = t.C[n + 1, m + 1]; S = t.S[n + 1, m + 1]
            Vp = V[n + 2, m + 2]; Wp = W[n + 2, m + 2]; Vm = V[n + 2, m]; Wm = W[n + 2, m]; V0 = V[n + 2, m + 1]; W0 = W[n + 2, m + 1]
            fac = (n - m + 2) * (n - m + 1)
            ax += 0.5 * (-C * Vp - S * Wp) + 0.5 * fac * (C * Vm + S * Wm)
            ay += 0.5 * (-C * Wp + S * Vp) + 0.5 * fac * (-C * Wm + S * Vm)
            az += (n - m + 1) * (-C * V0 - S * W0)
        end
    end
    k = t.mu / (t.R * t.R)
    return k * ax, k * ay, k * az
end

# ---- planet rotation -----------------------------------------------------------------------------------------------------

"""Rotation angle of the planet's prime meridian about the pole of the estimator's body-equator frame B (theta = 0 when the body-fixed x axis points along the x axis
of B), from the ENGINE's own frame routine: `SimulationModel.planet_frame_lpi(planet, et, SpiceEphemeridesModel())` (SPICE IAU_<BODY> with the pck00011 constants,
the routine the truth's gravity-harmonics effector and relative-velocity use). `B` is the J2000 -> B rotation of the planet record (`ctx.rot`, v_J2000 = B v_B),
`nothing` for the J2000 z pole. `kind = :uniform` uses theta0 + w t (diagnostic)."""
struct RotationModel
    kind::Symbol
    planet::Any
    B::Any
    et0::Float64
    model::Any
    theta0::Float64
    w::Float64
end

function body_angle(planet, B, model, et)
    l = SimulationModel.planet_frame_lpi(planet, et, model)
    xj = SV3(l[1, 1], l[1, 2], l[1, 3]); xb = B === nothing ? xj : SV3(B' * xj)
    return atan(xb[2], xb[1])
end

function RotationModel(ctx::EstContext, et0::Float64; kind::Symbol=:spice)
    model = SpiceEphemeridesModel(); th0 = body_angle(ctx.planet, ctx.rot, model, et0)
    w = wrap_pi(body_angle(ctx.planet, ctx.rot, model, et0 + 100.0) - th0) / 100.0
    return RotationModel(kind, ctx.planet, ctx.rot, et0, model, th0, w)
end

"""Uniform rotation (test model): theta(t) = theta0 + w t."""
uniform_rotation(theta0::Float64, w::Float64) = RotationModel(:uniform, nothing, nothing, 0.0, nothing, theta0, w)

"""(theta, w_rot) at `t_s` seconds after the run epoch."""
function rotation_at(rm::RotationModel, t_s::Float64)
    rm.kind === :uniform && return rm.theta0 + rm.w * t_s, rm.w
    th = body_angle(rm.planet, rm.B, rm.model, rm.et0 + t_s)
    w = wrap_pi(body_angle(rm.planet, rm.B, rm.model, rm.et0 + t_s + 100.0) - th) / 100.0
    return th, w
end

# ---- the per-pass kick ---------------------------------------------------------------------------------------------------

"""Options of the tesseral kick. `nodes`: Gauss-Legendre nodes in true anomaly for a FULL arc (half per half arc; numerical resolution, selftests a to c).
`Lt`, `Mt`: degree and order carried by the kick (0 = the whole file field); the carry rule always uses the whole file field. `budget_km`: carry iff D_T > budget_km
(ASSUMED, the 1 km oracle margin of the Kaula and Pi_T pre-registrations). `e_min`: below this eccentricity the 1/e terms (omega, mean anomaly) of the kick are dropped
(ASSUMED numerical guard; near-circular passes are inside the atmosphere). `force`: `true` carry, `false` do not (diagnostic), `nothing` = the rule.
`rotation`: `:spice` (engine routine at each pass time) or `:uniform` (diagnostic). `dt_terms`: carry the kick's timing terms. `post_drag_period`: the outbound leg lasts half the POST-drag period (addendum 1)."""
Base.@kwdef struct TesseralOpts
    nodes::Int = 64
    Lt::Int = 0
    Mt::Int = 0
    budget_km::Float64 = 1.0
    e_min::Float64 = 1.0e-3
    force::Union{Nothing, Bool} = nothing
    rotation::Symbol = :spice
    dt_terms::Bool = true
    post_drag_period::Bool = true   # outbound leg on the post-drag ellipse (addendum 1)
    base_period_without_kick::Bool = true   # that period from the state without the transient inbound half kick (addendum 2)
    nodes_dt::Int = 96             # nodes of the unit-kick quadrature of the carry rule (tesseral_param.py --nf default 96)
end

with_force(tt::TesseralOpts, f) = TesseralOpts(; (k => getfield(tt, k) for k in fieldnames(TesseralOpts) if k != :force)..., force=f)

"""Sum over the Gauss-Legendre nodes of the true-anomaly interval [f_lo, f_hi] of the per-time Gauss rates (da, de, di, dOmega, domega, dM) times dt, the first moment
sum(da_j tau_j) and sum(dM_j) needed for the timing, for the (n <= Lt, m <= Mt) non-zonal field. Elements frozen on the Kepler arc; theta = theta_p + w tau."""
function kick_arc(t::TesseralField, kw::KickWork, Lt::Int, Mt::Int, a, e, inc, raan, om, theta_p, wrot, f_lo, f_hi, nn::Int, e_min)
    xs, ws = gl01(nn); mu = t.mu
    F = OrbitFrame(a, e * cos(om), e * sin(om), inc, raan, mu); P, Q, Wv = F.P, F.Q, F.W
    p = a * (1 - e * e); h = sqrt(mu * p); Nm = sqrt(mu / a^3); b = a * sqrt(1 - e * e); sini = sin(inc); cosi = cos(inc)
    useew = e > e_min; usei = sini > 1e-8
    sda = 0.0; sde = 0.0; sdi = 0.0; sdO = 0.0; sdw = 0.0; sdM = 0.0; sdatau = 0.0
    span = f_hi - f_lo; vc = sqrt(mu / p)
    for j in 1:nn
        f = f_lo + span * xs[j]; wj = span * ws[j]
        sf, cf = sincos(f); rr = p / (1 + e * cf)
        Ef = 2 * atan(sqrt(1 - e) * sin(f / 2), sqrt(1 + e) * cos(f / 2))     # continuous in (-2 pi, 2 pi): no wrap at f = -pi
        tau = (Ef - e * sin(Ef)) / Nm
        th = theta_p + wrot * tau; ct, st = cos(th), sin(th)
        rhat = cf * P + sf * Q; shat = -sf * P + cf * Q
        pos = rr * rhat
        xb = ct * pos[1] + st * pos[2]; yb = -st * pos[1] + ct * pos[2]
        abx, aby, abz = nz_acc(t, kw, Lt, Mt, xb, yb, pos[3])
        ax = ct * abx - st * aby; ay = st * abx + ct * aby; az = abz
        aR = ax * rhat[1] + ay * rhat[2] + az * rhat[3]; aS = ax * shat[1] + ay * shat[2] + az * shat[3]; aW = ax * Wv[1] + ay * Wv[2] + az * Wv[3]
        wt = wj * rr * rr / h; u = om + f; su, cu = sincos(u)
        da = 2 * a * a / h * (e * sf * aR + p / rr * aS); de = (p * sf * aR + ((p + rr) * cf + rr * e) * aS) / h
        di = rr * cu * aW / h
        dO = usei ? rr * su * aW / (h * sini) : 0.0
        dw = useew ? ((-p * cf * aR + (p + rr) * sf * aS) / (h * e) - (usei ? rr * su * cosi * aW / (h * sini) : 0.0)) : 0.0
        dM = useew ? b / (a * h * e) * ((p * cf - 2 * e * rr) * aR - (p + rr) * sf * aS) : 0.0
        sda += da * wt; sde += de * wt; sdi += di * wt; sdO += dO * wt; sdw += dw * wt; sdM += dM * wt; sdatau += da * wt * tau
    end
    return (da=sda, de=sde, di=sdi, dO=sdO, dw=sdw, dM=sdM, datau=sdatau, Nm=Nm, a=a, e=e)
end

"""
    tesseral_pass_kick(el, t_peri, grav, rot; nodes=64, Lt=grav.L, Mt=Lt, half=true, f_start=-pi, e_min=1e-3)

First-order change of the osculating elements over the apoapsis-to-apoapsis Kepler arc centered on the periapsis of the pass at `t_peri` (s from the run epoch),
elements `el = (a, e, i, Omega, omega)` in the planet's body-equator frame, from the non-zonal gravity (degree 2..Lt, order 1..Mt) of `grav` and the rotation `rot`
(theta(t_peri) and w_rot). Returns `(da, de, di, dOmega, domega, drp, drp_half)` as in the derivation's sketch (drp = (1 - e) da - a de, drp_half its inbound half) and,
for the stepper, the inbound-half elements `(dah, deh, dih, dOh, dwh)`, the timing shifts `dt_in` (periapsis time) and `dt_arc` (next apoapsis time) in seconds. `f_start`:
true anomaly at which the first arc starts (the stepper's first pass starts at the initial anomaly); `half=false` leaves the inbound-half fields NaN.
"""
function tesseral_pass_kick(el::NTuple{5, Float64}, t_peri::Float64, grav::TesseralField, rot::RotationModel; nodes::Int=64, Lt::Int=grav.L, Mt::Int=Lt,
                            half::Bool=true, f_start::Float64=-float(π), e_min::Float64=1.0e-3)
    Lt = min(Lt, grav.L); Mt = min(Mt, Lt); a, e, inc, raan, om = el
    theta, wrot = rotation_at(rot, t_peri); kw = kick_work(grav.L)
    nin = max(8, round(Int, nodes * (-f_start) / (2π))); nout = max(8, nodes ÷ 2)
    kin = kick_arc(grav, kw, Lt, Mt, a, e, inc, raan, om, theta, wrot, f_start, 0.0, nin, e_min)
    kout = kick_arc(grav, kw, Lt, Mt, a, e, inc, raan, om, theta, wrot, 0.0, float(π), nout, e_min)
    Nm = kin.Nm; Tend = π / Nm      # tau at the arc end (apoapsis)
    da = kin.da + kout.da; de = kin.de + kout.de
    # timing: int Delta a dt over the arc = sum da_j (tau_end - tau_j); the inbound part alone ends at tau = 0
    int_da_full = Tend * da - (kin.datau + kout.datau); int_da_in = -kin.datau
    dt_in = 1.5 * int_da_in / a - kin.dM / Nm
    dt_arc = 1.5 * int_da_full / a - (kin.dM + kout.dM) / Nm
    drp = (1 - e) * da - a * de; drp_half = (1 - e) * kin.da - a * kin.de
    return (da=da, de=de, di=kin.di + kout.di, dOmega=kin.dO + kout.dO, domega=kin.dw + kout.dw, drp=drp, drp_half=half ? drp_half : NaN,
            dah=half ? kin.da : NaN, deh=half ? kin.de : NaN, dih=half ? kin.di : NaN, dOh=half ? kin.dO : NaN, dwh=half ? kin.dw : NaN,
            dt_in=dt_in, dt_arc=dt_arc, theta=theta, wrot=wrot)
end

# ---- the carry rule: Pi_T and D_T ------------------------------------------------------------------------------------------------

"""Per-term unit kicks Y[n+1, m+1] (complex): the first-order periapsis-radius kick over one apoapsis-to-apoapsis arc of a unit FULLY NORMALIZED coefficient,
kick_rp(n, m) = Cbar Re Y + Sbar Im Y (tesseral_lib.term_kicks, one Cunningham recursion per node for all (n, m)); GL nodes in true anomaly over [-pi, pi]."""
function unit_kicks_rp(t::TesseralField, a, e, inc, om, wrot, nn::Int)
    L = t.L; kw = kick_work(L); xs, ws = gl01(nn); mu = t.mu; R = t.R
    F = OrbitFrame(a, e * cos(om), e * sin(om), inc, 0.0, mu); P, Q, Wv = F.P, F.Q, F.W
    p = a * (1 - e * e); h = sqrt(mu * p); Nm = sqrt(mu / a^3); sini = sin(inc)
    Y = zeros(ComplexF64, L + 1, L + 1); g0 = mu / (R * R)
    for j in 1:nn
        f = -π + 2π * xs[j]; wj = 2π * ws[j]
        sf, cf = sincos(f); rr = p / (1 + e * cf)
        Ef = 2 * atan(sqrt(1 - e) * sin(f / 2), sqrt(1 + e) * cos(f / 2)); tau = (Ef - e * sin(Ef)) / Nm
        th = wrot * tau; ct, st = cos(th), sin(th)
        rhat = cf * P + sf * Q; shat = -sf * P + cf * Q; pos = rr * rhat
        xb = ct * pos[1] + st * pos[2]; yb = -st * pos[1] + ct * pos[2]
        cunningham!(kw, R, xb, yb, pos[3], L, L); V, W = kw.V, kw.W
        wt = wj * rr * rr / h
        @inbounds for n in 2:L, m in 1:n
            Nf = t.Nf[n + 1, m + 1]; g = g0 * Nf; fac = (n - m + 2) * (n - m + 1)
            Vp = V[n + 2, m + 2]; Wp = W[n + 2, m + 2]; Vm = V[n + 2, m]; Wm = W[n + 2, m]; V0 = V[n + 2, m + 1]; W0 = W[n + 2, m + 1]
            # body-frame accelerations of a unit C (index 1) and a unit S (index 2)
            cx = g * (-0.5 * Vp + 0.5 * fac * Vm); cy = g * (-0.5 * Wp - 0.5 * fac * Wm); cz = g * (-(n - m + 1) * V0)
            sx = g * (-0.5 * Wp + 0.5 * fac * Wm); sy = g * (0.5 * Vp + 0.5 * fac * Vm); sz = g * (-(n - m + 1) * W0)
            for (bx, by, bz, isS) in ((cx, cy, cz, false), (sx, sy, sz, true))
                ax = ct * bx - st * by; ay = st * bx + ct * by; az = bz
                aR = ax * rhat[1] + ay * rhat[2] + az * rhat[3]; aS = ax * shat[1] + ay * shat[2] + az * shat[3]
                da = 2 * a * a / h * (e * sf * aR + p / rr * aS)
                de = (p * sf * aR + ((p + rr) * cf + rr * e) * aS) / h
                drp_rate = (1 - e) * da - a * de
                isS ? (Y[n + 1, m + 1] += 1im * drp_rate * wt) : (Y[n + 1, m + 1] += drp_rate * wt)
            end
        end
    end
    return Y
end

"""Accumulation gain G_m^2 (pre-registered closed form of tesseral_param.py `gain`)."""
function tesseral_gain2(m::Int, nu0, nuend, Kh::Int)
    g = min(1.0 / (2 * sin(π * m * nu0)^2), Kh^2 / 3.0)
    nup = (nuend - nu0) / max(Kh - 1, 1)
    if abs(nup) > 0
        lo, hi = minmax(m * nu0, m * nuend)
        for j in (floor(Int, lo) + 1):(ceil(Int, hi) - 1)
            Kc = abs(m * nu0 - j) / (m * abs(nup)); g += (1 / (m * abs(nup))) * max(0.0, 1 - Kc / Kh)
        end
    end
    return g
end

"""D_T (km): horizon rms of the tesseral periapsis-radius change (Pi_T parameter, tesseral_param.py `cell_param`): kappa_nm = sigma_n |Y_nm| (m), G_m the accumulation gain
of the commensurability sweep nu0 -> nu_end over K_h passes; D_T = sqrt(sum (kappa_nm G_m)^2)."""
function tesseral_DT(t::TesseralField, a1, e1, i1, w1, spin, nu0, nuend, Kh::Int, nodes::Int)
    Y = unit_kicks_rp(t, a1, e1, i1, w1, spin, nodes); tot = 0.0
    for m in 1:t.L
        G2 = tesseral_gain2(m, nu0, nuend, Kh); s = 0.0
        for n in max(2, m):t.L; s += (t.sigma[n + 1] * abs(Y[n + 1, m + 1]))^2; end
        tot += s * G2
    end
    return sqrt(tot) / 1.0e3
end

"""D_T of a v3 prediction `p` (the zonal reference history): first-pass osculating elements, nu from the pass-time differences."""
function dt_of_prediction(t::TesseralField, p::Prediction, spin, nodes::Int)
    Kh = length(p.t_s); Kh >= 1 || return 0.0
    nu = Kh >= 2 ? spin * diff(p.t_s) / (2π) : Float64[spin * sqrt(p.a_m[1]^3 / t.mu)]
    return tesseral_DT(t, p.a_m[1], p.e[1], p.i_rad[1], p.omega_rad[1], spin, nu[1], nu[end], max(Kh, 2), nodes)
end

# ---- stepper -----------------------------------------------------------------------------------------------------------------

"""Add first-order element increments to a mean state (a, k, h, i, Omega); (k, h) = e (cos w, sin w) updated linearly, so zero increments leave the state bit-identical."""
@inline function add_kick(m::SV5, da, de, di, dO, dw)
    k, h = m[2], m[3]; e = hypot(k, h); e > 0 || return SV5(m[1] + da, m[2] + de, m[3], m[4] + di, m[5] + dO)
    c, s = k / e, h / e
    return SV5(m[1] + da, k + de * c - h * dw, h + de * s + k * dw, m[4] + di, m[5] + dO)
end

"""Elements `(a, e, i, Omega, omega)` of a mean state."""
@inline mean_elements(m::SV5) = (m[1], hypot(m[2], m[3]), m[4], m[5], atan(m[3], m[2]))

"""
`predict_history3` plus the tesseral kick. Returns `(prediction, info)`; `info = (ns, D, DT_km, carried)`. `carried == false`: the prediction IS the v3 one (same object).
`tf`: `TesseralField` of the truth gravity model, `tt`: `TesseralOpts`. Self-contained keywords as v3.
"""
function predict_history4(ctx::EstContext, s::EstSettings, o::V2Opts, tab::SPTable, zf::ZonalField, zo::ZonalOpts, tf::TesseralField, tt::TesseralOpts,
                          seed_alt_m, inc_rad, raan_rad, omega_rad, et0;
                          max_passes::Int=typemax(Int), rho_fn=fit_rho(ctx), bodies_rb=(et -> body_positions(ctx, et)), init_override=nothing, DT_override=nothing, theta_times=nothing)
    base, info3 = predict_history3(ctx, s, o, tab, zf, zo, seed_alt_m, inc_rad, raan_rad, omega_rad, et0; max_passes=max_passes, rho_fn=rho_fn, bodies_rb=bodies_rb, init_override=init_override)
    tt.force === false && return base, (info3[1], info3[2], NaN, false)
    spin = norm(ctx.omega_e)
    DT = DT_override !== nothing ? DT_override : (tt.force === true ? NaN : dt_of_prediction(tf, base, spin, tt.nodes_dt))
    carry = tt.force === true ? true : DT > tt.budget_km
    carry || return base, (info3[1], info3[2], DT, false)
    out = Prediction()
    Lt = tt.Lt > 0 ? min(tt.Lt, tf.L) : tf.L; Mt = tt.Mt > 0 ? min(tt.Mt, Lt) : Lt
    rotm = RotationModel(ctx, et0; kind=tt.rotation)
    init = init_override === nothing ? initial_osculating(ctx, seed_alt_m, inc_rad, raan_rad, omega_rad) : init_override
    ns, D = info3[1], info3[2]; Js = [zonal_J(zf, n) for n in ns]
    nz = isempty(ns) ? 0 : max(zo.nodes_min, zo.nodes_per_degree * ns[end])
    ivec = SV5(init.a, init.e * cos(init.omega), init.e * sin(init.omega), init.i, init.raan)
    rbs0 = s.tides ? bodies_rb(et0) : SV3[]
    d0 = zero(SV5); m = ivec; mm = init
    for _ in 1:(o.sp_tide === :none ? 1 : 3)
        el = vec_el(ivec - d0)
        mm = osculating_to_mean((a=el.a, e=el.e, i=el.i, raan=el.raan, omega=el.omega, nu=init.nu); mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements
        m = mean_state(mm, ctx)
        d0 = tide_sp(m, ctx, o, tab, rbs0, :apo)
    end
    M0 = true_to_mean_anomaly(mm.e, mm.nu)
    f_start0 = mm.nu - 2π; f_start0 > 0 && (f_start0 -= 2π)
    t_apo = 0.0
    for kpass in 1:max_passes
        mdot = mean_motion_rate(m, ctx)
        T = 2π / mdot
        dt_in = kpass == 1 ? (2π - M0) / mdot : 0.5 * T
        dt_out = 0.5 * T
        t_k0 = t_apo + dt_in
        t_k0 > s.horizon_s && break
        t_rot = (theta_times !== nothing && kpass <= length(theta_times)) ? theta_times[kpass] : t_k0     # UNFROZEN diagnostic: rotation angle at the TRUTH pass time
        kk = try tesseral_pass_kick(mean_elements(m), t_rot, tf, rotm; nodes=tt.nodes, Lt=Lt, Mt=Mt, f_start=kpass == 1 ? f_start0 : -float(π), e_min=tt.e_min) catch; nothing end
        if kk === nothing || !isfinite(kk.drp)
            out.decayed = true; out.t_decay_s = t_k0; out.ended = :invalid; break
        end
        dtk_in = tt.dt_terms ? kk.dt_in : 0.0; dtk_arc = tt.dt_terms ? kk.dt_arc : 0.0
        et_body = et0 + t_k0 + (s.body_epoch === :apoapsis ? dt_out : 0.0)
        rbs = bodies_rb(et_body)
        i1, dO1, dw1 = mean_rates_r3(m, ctx, s, rbs, ns, Js, nz)
        i2, dO2, dw2 = mean_rates_r3(advance(m, i1, dO1, dw1, dt_in), ctx, s, rbs, ns, Js, nz)
        m_pre = advance(m, i2, dO2, dw2, dt_in)
        m_pre = add_kick(m_pre, kk.dah, kk.deh, kk.dih, kk.dOh, kk.dwh)               # inbound half of the tesseral kick (physical periapsis radius)
        t_k = t_k0 + dtk_in
        d3 = tide_sp(m_pre, ctx, o, tab, s.tides ? rbs : SV3[], :peri)
        ov = try osc_vec2(m_pre, ctx, o, tab, d3) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        osc = vec_el(ov)
        din = zero(SV5); dout = zero(SV5)
        if s.drag; din, dout = drag_deltas(osc, ctx, s, rho_fn); end
        fin = rp_form(ov, din); fout = rp_form(ov, dout)
        m_peri = try shift_mean2(m_pre, fin, ctx, o, tab, d3) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        osc = try vec_el(osc_vec2(m_peri, ctx, o, tab, d3)) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        rp = osc.a * (1 - osc.e)
        if rp - ctx.Rp_e <= s.impact_altitude_m; out.decayed = true; out.t_decay_s = t_k; out.ended = :impact; break; end
        rv = elements_to_cartesian((a=osc.a, e=osc.e, i=osc.i, raan=osc.raan, omega=osc.omega, nu=0.0); mu=ctx.mu)
        r = SV3(rv[1]); v = SV3(rv[2])
        alt, lat = geodetic_altitude_latitude(r, ctx.planet)
        vrel = norm(v - cross(ctx.omega_e, r)); rho = corridor_density(ctx.fit, alt)
        hq, hl = band_edges(ctx, vrel)
        for (arr, val) in zip((out.t_s, out.alt_m, out.radius_m, out.lat_rad, out.vrel_mps, out.rho, out.q_w_cm2, out.a_m, out.e, out.i_rad, out.raan_rad, out.omega_rad, out.h_qhigh_m, out.h_qlow_m),
                              (t_k, alt, rp, lat, vrel, rho, 0.5 * rho * vrel^3 / ctx.w2m, osc.a, osc.e, osc.i, osc.raan, osc.omega, hq, hl))
            push!(arr, val)
        end
        m = try
            m_post = shift_mean2(m_peri, fout, ctx, o, tab, d3)
            # timing (addendum 1 of the plan): the outbound leg runs on the POST-drag ellipse, so its duration is half the post-drag period; v3 uses the pre-pass period
            # for both legs, which accumulates 0.5 (T_k - T_k+1) per pass (about 450 s per pass at Mars under strong drag, 6,000 s after 20 d)
            # the kick's own timing terms (dt_in, dt_arc) already contain the effect of the transient osculating change of a along the arc, so the base period must not
            # contain the inbound half kick (addendum 2: with it the leg is 1.5 da_half/a * T/2 = 70 s too long at Mars)
            tt.post_drag_period && (dt_out = 0.5 * 2π / mean_motion_rate(tt.base_period_without_kick ? add_kick(m_post, -kk.dah, -kk.deh, -kk.dih, -kk.dOh, -kk.dwh) : m_post, ctx))
            advance(m_post, i2, dO2, dw2, dt_out)
        catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        m = add_kick(m, kk.da - kk.dah, kk.de - kk.deh, kk.di - kk.dih, kk.dOmega - kk.dOh, kk.domega - kk.dwh)         # outbound half
        t_apo = t_k + dt_out + (dtk_arc - dtk_in)
    end
    return out, (ns, D, DT, true)
end

"""Self-calibrated seed with the v4 predictor: the v3 seed decides the carry rule; when carried, the seed is iterated with the kick (inbound half) included in the first pass."""
function calibrate_seed4(ctx::EstContext, s::EstSettings, o::V2Opts, tab::SPTable, zf::ZonalField, zo::ZonalOpts, tf::TesseralField, tt::TesseralOpts, inc_rad, raan_rad, omega_rad, et0;
                         seed0=140.0e3, tol_m=1.0, max_iter=8)
    s3 = calibrate_seed3(ctx, s, o, tab, zf, zo, inc_rad, raan_rad, omega_rad, et0; seed0=seed0, tol_m=tol_m, max_iter=max_iter)
    isfinite(s3) || return s3
    p3, info = predict_history4(ctx, s, o, tab, zf, zo, tf, tt, s3, inc_rad, raan_rad, omega_rad, et0)
    info[4] || return s3
    DT = info[3]; q_target = sqrt(ctx.q_low * ctx.q_high); seed = s3
    for _ in 1:max_iter
        p, _ = predict_history4(ctx, s, o, tab, zf, zo, tf, with_force(tt, true), seed, inc_rad, raan_rad, omega_rad, et0; max_passes=1)
        isempty(p.t_s) && return NaN
        h_target = altitude_of_density_fast(ctx.fit, 2 * q_target * ctx.w2m / p.vrel_mps[1]^3)
        abs(h_target - p.alt_m[1]) < tol_m && return seed
        seed += h_target - p.alt_m[1]
    end
    return seed
end
