# An in-situ density-scale-factor filter, in the style of Tracy, Falcone and
# Manchester, "Robust Entry Guidance with Atmospheric Adaptation".
#
# That work estimates a single scalar
#
#     k_rho = rho_observed / rho_nominal                              (their Eq 22)
#
# with a square-root EKF whose state is `[r, v, sigma, k_rho]`, driven by
# navigation position/velocity measurements, and feeds `rho = k_rho * rho_approx`
# into CPEG. Estimating the ratio rather than the density itself is what keeps
# the filter numerically sane across five orders of magnitude of density.
#
# What is reproduced here, and what is not
# ----------------------------------------
# This study's EDL "corridor" is a kinematic path in (lat, lon, alt) — it does
# not respond to density, so a filter driven by *trajectory deviation* has
# nothing to observe. The comparison therefore drives the filter with drag
# acceleration directly:
#
#     a_drag = 0.5 * rho * V^2 / beta,      beta = m / (C_D * A)
#
# which is what an accelerometer measures and what the vehicle's response
# integrates. This is a **simplification that favours the filter**: a direct
# acceleration measurement is strictly more informative about instantaneous
# density than the position/velocity measurements the paper uses, which see
# density only through accumulated dynamics. Read the filter numbers below as an
# upper bound on what the paper's architecture would achieve here.
#
# Retained from the paper: the scalar `k_rho` state, the log-domain formulation
# for numerical range, the nominal-times-ratio structure, and the fact that
# observability collapses where dynamic pressure is low.
#
# Not reproduced: the full `[r, v, sigma, k_rho]` state, the square-root
# covariance factorisation (a scalar state does not need it), bank-angle
# dynamics, and CPEG itself.

"""
    EntryVehicle(; mass_kg, ref_area_m2, cd)

Blunt-body entry vehicle. Defaults are the MSL-class numbers the paper uses: a
70-degree sphere-cone of `2.25 m` base radius and `2800 kg`, giving a ballistic
coefficient near `120 kg/m^2`.
"""
struct EntryVehicle
    mass_kg::Float64
    ref_area_m2::Float64
    cd::Float64
end

function EntryVehicle(; mass_kg::Float64=2800.0, base_radius_m::Float64=2.25, cd::Float64=1.45)
    (mass_kg > 0.0 && base_radius_m > 0.0 && cd > 0.0) ||
        throw(ArgumentError("Entry vehicle parameters must be positive."))
    return EntryVehicle(mass_kg, pi * base_radius_m^2, cd)
end

"""
    ballistic_coefficient(vehicle) -> Float64

`m / (C_D A)`, in kg/m^2.
"""
@inline ballistic_coefficient(v::EntryVehicle)::Float64 = v.mass_kg / (v.cd * v.ref_area_m2)

"""
    DensityFilterConfig(; ...)

`accel_noise_mps2` is the accelerometer noise floor. It is what makes density
unobservable high up, where `a_drag` falls below it — the same effect the paper
reports as large `k_rho` uncertainty at high altitude.

`process_std_per_sqrt_s` is the random-walk rate on `log k_rho`: how fast the
filter is willing to believe the density ratio changes.
"""
struct DensityFilterConfig
    accel_noise_mps2::Float64
    accel_rel_noise::Float64
    process_std_per_sqrt_s::Float64
    init_log_k_std::Float64
end

function DensityFilterConfig(;
    accel_noise_mps2::Float64=1.0e-3,
    accel_rel_noise::Float64=0.02,
    process_std_per_sqrt_s::Float64=0.01,
    init_log_k_std::Float64=0.15,
)
    all(>(0.0), (accel_noise_mps2, accel_rel_noise, process_std_per_sqrt_s, init_log_k_std)) ||
        throw(ArgumentError("Filter noise parameters must be positive."))
    return DensityFilterConfig(accel_noise_mps2, accel_rel_noise, process_std_per_sqrt_s, init_log_k_std)
end

"""
    corridor_speeds(points) -> Vector{Float64}

Path speed implied by the kinematic corridor, from great-circle plus vertical
displacement between consecutive points. The default EDL corridor works out near
`4 km/s`, which is a physically sensible entry speed.
"""
function corridor_speeds(points::Vector{TrajectoryPoint})::Vector{Float64}
    n = length(points)
    n >= 2 || return fill(0.0, n)
    speeds = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        a = points[max(i - 1, 1)]
        b = points[min(i + 1, n)]
        dt = (b.elapsed_s - a.elapsed_s)
        dt <= 0.0 && (dt = 1.0)
        dh_km = chordal_distance_km(a.lat_deg, a.lon_deg, b.lat_deg, b.lon_deg)
        dv_km = (b.alt_m - a.alt_m) * 1e-3
        speeds[i] = 1.0e3 * sqrt(dh_km^2 + dv_km^2) / dt
    end
    return speeds
end

"""
    run_density_filter(cfg, vehicle, points, speeds, rho_truth, rho_nominal; rng)
        -> (k_hat, k_std, observable)

Scalar EKF on `x = log k_rho`, marched along the corridor. Returns the estimated
ratio, its one-sigma spread, and whether each step carried usable information
(drag acceleration above the noise floor).

`rho_nominal` is the onboard model the ratio is taken against — the GRAM prior
in the baseline comparison, or a GP-corrected profile when the two are combined.
"""
function run_density_filter(
    cfg::DensityFilterConfig,
    vehicle::EntryVehicle,
    points::Vector{TrajectoryPoint},
    speeds::Vector{Float64},
    rho_truth::Vector{Float64},
    rho_nominal::Vector{Float64};
    rng::AbstractRNG=MersenneTwister(0),
)
    n = length(points)
    (n == length(speeds) == length(rho_truth) == length(rho_nominal)) ||
        throw(ArgumentError("Filter inputs must share length."))
    beta = ballistic_coefficient(vehicle)

    x = 0.0                      # log k_rho, initialised at "model is right"
    P = cfg.init_log_k_std^2
    k_hat = Vector{Float64}(undef, n)
    k_std = Vector{Float64}(undef, n)
    observable = falses(n)

    @inbounds for i in 1:n
        if i > 1
            dt = max(points[i].elapsed_s - points[i - 1].elapsed_s, 0.0)
            P += cfg.process_std_per_sqrt_s^2 * dt          # random-walk predict
        end

        a_true = 0.5 * rho_truth[i] * speeds[i]^2 / beta
        sigma_meas = sqrt((cfg.accel_rel_noise * a_true)^2 + cfg.accel_noise_mps2^2)
        z = a_true + sigma_meas * randn(rng)

        # h(x) = 0.5 * exp(x) * rho_nominal * V^2 / beta,  dh/dx = h(x)
        h = 0.5 * exp(x) * rho_nominal[i] * speeds[i]^2 / beta
        if h > 0.0 && a_true > cfg.accel_noise_mps2
            observable[i] = true
            S = h * P * h + sigma_meas^2
            K = P * h / S
            x += K * (z - h)
            P = max((1.0 - K * h) * P, 1.0e-12)
        end
        k_hat[i] = exp(x)
        k_std[i] = exp(x) * sqrt(P)
    end
    return k_hat, k_std, observable
end

"""
    forward_prediction_error(rho_pred_at, rho_truth, weights; lookahead)
        -> Float64

Weighted RMS log-density error of the profile a guidance call at index `j` would
predict for the points still ahead of it.

This is the quantity a predictor-corrector actually consumes: CPEG re-simulates
to parachute deploy at every call, so it needs density along the *whole
remaining* trajectory, not just at the vehicle. A reactive filter can only offer
its current estimate held forward, which is why this metric separates it from a
spatial model. `rho_pred_at(j, i)` gives the density a call at `j` predicts for
point `i`.
"""
function forward_prediction_error(
    rho_pred_at,
    rho_truth::Vector{Float64},
    weights::Vector{Float64};
    lookahead::Int=0,
)::Float64
    n = length(rho_truth)
    num = 0.0
    den = 0.0
    for j in 1:(n - 1)
        stop = lookahead > 0 ? min(j + lookahead, n) : n
        for i in (j + 1):stop
            p = rho_pred_at(j, i)
            (isfinite(p) && p > 0.0) || continue
            num += weights[i] * (log(p / rho_truth[i]))^2
            den += weights[i]
        end
    end
    return den > 0.0 ? sqrt(num / den) : NaN
end
