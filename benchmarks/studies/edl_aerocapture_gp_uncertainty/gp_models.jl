struct Parameterization
    name::String
end

const PARAMETERIZATIONS = (
    Parameterization("log_density"),
    Parameterization("density_scale_factor"),
    Parameterization("log_density_scale_factor"),
)

# `gram_exponential` is GRAM's own correlation model rather than a generic
# kernel: a separable exponential on physical distance,
# `exp(-(dx/Lh + dz/Lz))`, with `Lh` and `Lz` read from GRAM's tabulated
# altitude-dependent scales (see gram_correlation.jl). It ignores the
# `--lat-scale-deg` / `--lon-scale-deg` / `--alt-scale-km` knobs, which the
# other three kernels use, and takes its scales from `gram_ref_alt_km`.
const KERNEL_NAMES = ("squared_exponential", "matern32", "matern52", "gram_exponential")

# Explicit mean-function bases for the residual GP.
#
# A zero-mean stationary GP decays to "no correction" more than a length scale
# from its data, so it structurally cannot carry information from the sensed
# band into the predicted band. The most actionable thing an aerocapture pass
# can report -- a bulk density offset relative to the prior -- is exactly such a
# constant. `constant` and `linear_alt` fit that component as a parametric mean
# by generalized least squares (universal kriging), so it extrapolates, and
# inflate the predictive variance where the basis itself is extrapolated.
const MEAN_BASIS_NAMES = ("none", "constant", "linear_alt")

struct GPConfig
    kernel::String
    lat_scale_deg::Float64
    lon_scale_deg::Float64
    alt_scale_km::Float64
    jitter::Float64
    mean_basis::String
    gram_ref_alt_km::Float64
    time_scale_s::Float64   # Inf ignores the time axis entirely
    amplitude_mode::String
end

# Backwards-compatible constructor: no mean function reproduces the original
# zero-mean residual GP exactly.
GPConfig(kernel, lat_scale_deg, lon_scale_deg, alt_scale_km, jitter) =
    GPConfig(kernel, lat_scale_deg, lon_scale_deg, alt_scale_km, jitter, "none")

GPConfig(kernel, lat_scale_deg, lon_scale_deg, alt_scale_km, jitter, mean_basis) =
    GPConfig(kernel, lat_scale_deg, lon_scale_deg, alt_scale_km, jitter, mean_basis, DEFAULT_GRAM_REF_ALT_KM)

GPConfig(kernel, lat_scale_deg, lon_scale_deg, alt_scale_km, jitter, mean_basis, gram_ref_alt_km) =
    GPConfig(
        kernel, lat_scale_deg, lon_scale_deg, alt_scale_km, jitter, mean_basis,
        gram_ref_alt_km, gram_time_scale_s(gram_ref_alt_km),
    )

GPConfig(kernel, lat_scale_deg, lon_scale_deg, alt_scale_km, jitter, mean_basis, gram_ref_alt_km, time_scale_s) =
    GPConfig(
        kernel, lat_scale_deg, lon_scale_deg, alt_scale_km, jitter, mean_basis,
        gram_ref_alt_km, time_scale_s, DEFAULT_AMPLITUDE_MODE,
    )

# How the kernel gets its scale.
#
# `stationary` is the original: one global amplitude, `var(detrended residuals)`,
# and the GRAM prior variance re-introduced afterwards by the fusion heuristic in
# `predict_density_profile`.
#
# `prior_scaled` folds the prior variance into the kernel instead,
#
#     k(x, x") = lambda^2 * s(x) * s(x") * rho(x, x")
#
# with `s(x)` the GRAM ensemble relative sigma at `x` and `rho` the correlation
# function. `D K D` with `D = diag(s)` is positive definite whenever `K` is. The
# GP prior variance at `x` is then exactly `lambda^2 s(x)^2`, so the posterior
# returns to the GRAM prior far from data and tightens only near it, and no
# fusion step is needed at all: the predictive variance *is* the GP posterior
# variance. `lambda` is a single inflation factor for GRAM ensemble
# mis-dispersion, fitted by marginal likelihood.
#
# The two differ in the mean as well as the variance. In the noise-free limit
# `prior_scaled` gives `mu(x*) = s(x*) * interp(y_i / s(x_i))`: the GP transfers
# residuals in units of the local prior sigma rather than in absolute
# log-density, so "the atmosphere is +1.5 sigma" propagates instead of
# "+13% at 90 km" propagating unchanged to 35 km where the real spread is 1.2%.
const AMPLITUDE_MODE_NAMES = ("stationary", "prior_scaled")

# The low-level default preserves the original behaviour so existing callers and
# the published stationary results are reproducible; the study driver defaults to
# `prior_scaled` through `--amplitude-mode`.
const DEFAULT_AMPLITUDE_MODE = "stationary"

# Bounds on the fitted inflation factor.
#
# The floor is not cosmetic. `lambda` is meant to correct GRAM ensemble
# mis-dispersion, but the marginal likelihood cannot separate that from the
# fraction of the training pass that carries no signal, and on this study the
# two are confounded: measured on real MERRA-2 truth, `lambda` fits to 0.55 when
# 19% of the pass lies in the band where truth is defined, 1.77 at 43%, and 3.02
# at 85%. A `lambda` below 1 asserts the atmosphere varies *less* than GRAM's own
# ensemble says — a claim no pass can support through 5% measurement noise, and
# the direct cause of the +1218 weighted NLPD seen at the 0.55 fit. Clamping at 1
# makes the estimator fall back on the prior's own spread instead.
# That argument holds only for a zero-mean GP. With a mean basis, `lambda`
# describes the residual *about the fitted mean*, which legitimately is smaller
# than the raw ensemble spread because the parametric term has already explained
# part of it — and the universal-kriging variance adds the mean coefficients'
# own uncertainty back, so the total stays honest. The floor is therefore applied
# only when there is no mean function.
const DEFAULT_LAMBDA_BOUNDS = (1.0, 10.0)           # zero-mean GP
const DEFAULT_LAMBDA_BOUNDS_WITH_MEAN = (0.1, 10.0) # a mean basis is fitted

# GRAM's scales vary with altitude, but a kernel with position-dependent scales
# is not positive definite without a Gibbs-style normalization that has no clean
# form for the exponential. They are therefore frozen at one reference altitude,
# which is a mild approximation over the band that matters: GRAM tabulates
# Lh = 2240-2318 km and Lz = 17.4-19.3 km across 30-60 km.
const DEFAULT_GRAM_REF_ALT_KM = 35.0

struct FittedGP
    X::Matrix{Float64}
    alpha::Vector{Float64}
    chol::LinearAlgebra.Cholesky{Float64, Matrix{Float64}}
    amplitude2::Float64         # stationary: kernel amplitude; prior_scaled: lambda^2
    cfg::GPConfig
    H::Matrix{Float64}          # n x p training basis
    Kinv_H::Matrix{Float64}     # n x p, K^-1 H
    A_inv::Matrix{Float64}      # p x p, (H' K^-1 H)^-1
    beta::Vector{Float64}       # p, GLS mean-function coefficients
    center_alt_km::Float64
    prior_sigma::Vector{Float64}  # training-point GRAM sigma; empty when stationary
    log_marginal::Float64
    signal_fraction::Float64      # share of residual variance above the noise floor
    lambda_at_bound::Bool
end

@inline function _clip_density(rho::Float64)::Float64
    return max(rho, 1.0e-18)
end

@inline function _prior_target(par::Parameterization, gram::GramSample)::Tuple{Float64, Float64}
    mu = _clip_density(gram.mean_density)
    rel_sigma = gram.std_density / mu
    prior_var = max(rel_sigma^2, 1.0e-10)
    if par.name == "log_density"
        return log(mu), prior_var
    elseif par.name == "density_scale_factor"
        return 1.0, prior_var
    elseif par.name == "log_density_scale_factor"
        return 0.0, prior_var
    end
    throw(ArgumentError("Unsupported parameterization $(par.name)."))
end

@inline function _measurement_target(par::Parameterization, rho_meas::Float64, gram::GramSample)::Tuple{Float64, Float64}
    mu = _clip_density(gram.mean_density)
    rho = _clip_density(rho_meas)
    if par.name == "log_density"
        return log(rho), 0.05^2
    elseif par.name == "density_scale_factor"
        sigma = 0.05 * rho / mu
        return rho / mu, max(sigma^2, 1.0e-10)
    elseif par.name == "log_density_scale_factor"
        return log(rho / mu), 0.05^2
    end
    throw(ArgumentError("Unsupported parameterization $(par.name)."))
end

@inline function _physical_moments(par::Parameterization, mu_target::Float64, var_target::Float64, gram::GramSample)::Tuple{Float64, Float64}
    sigma2 = max(var_target, 1.0e-12)
    if par.name == "log_density"
        mean_rho = exp(mu_target + 0.5 * sigma2)
        var_rho = (exp(sigma2) - 1.0) * exp(2.0 * mu_target + sigma2)
        return mean_rho, sqrt(max(var_rho, 1.0e-18))
    elseif par.name == "density_scale_factor"
        mu = max(gram.mean_density * mu_target, 1.0e-18)
        sigma = abs(gram.mean_density) * sqrt(sigma2)
        return mu, max(sigma, 1.0e-18)
    elseif par.name == "log_density_scale_factor"
        scale_mean = exp(mu_target + 0.5 * sigma2)
        scale_var = (exp(sigma2) - 1.0) * exp(2.0 * mu_target + sigma2)
        mean_rho = gram.mean_density * scale_mean
        sigma_rho = abs(gram.mean_density) * sqrt(max(scale_var, 1.0e-18))
        return max(mean_rho, 1.0e-18), max(sigma_rho, 1.0e-18)
    end
    throw(ArgumentError("Unsupported parameterization $(par.name)."))
end

# Absolute seconds, so aerocapture and EDL points share one time origin.
# `TrajectoryPoint.elapsed_s` restarts at zero on each corridor and cannot be
# used here.
@inline function _feature_row(point::TrajectoryPoint)::NTuple{4, Float64}
    return (point.lat_deg, point.lon_deg, point.alt_m * 1e-3, Dates.value(point.dt) * 1e-3)
end

@inline function _scaled_distance2(cfg::GPConfig, x::NTuple{4, Float64}, y::NTuple{4, Float64})::Float64
    dlat = (x[1] - y[1]) / cfg.lat_scale_deg
    dlon = (x[2] - y[2]) / cfg.lon_scale_deg
    dalt = (x[3] - y[3]) / cfg.alt_scale_km
    r2 = dlat * dlat + dlon * dlon + dalt * dalt
    if isfinite(cfg.time_scale_s)
        dt = (x[4] - y[4]) / cfg.time_scale_s
        r2 += dt * dt
    end
    return r2
end

@inline function _correlation(cfg::GPConfig, x::NTuple{4, Float64}, y::NTuple{4, Float64})::Float64
    if cfg.kernel == "gram_exponential"
        lh_km, lz_km = gram_correlation_scales(cfg.gram_ref_alt_km)
        dx_km = chordal_distance_km(x[1], x[2], y[1], y[2])
        dz_km = abs(x[3] - y[3])
        # GRAM's full space-time form: exp(-(dx/Lh + dz/Lz + dt/Lt)). Each factor
        # is a 1-D exponential kernel, so the product stays positive definite.
        u = dx_km / lh_km + dz_km / lz_km
        if isfinite(cfg.time_scale_s)
            u += abs(x[4] - y[4]) / cfg.time_scale_s
        end
        return exp(-u)
    end
    r2 = _scaled_distance2(cfg, x, y)
    if cfg.kernel == "squared_exponential"
        return exp(-0.5 * r2)
    end
    r = sqrt(r2)
    if cfg.kernel == "matern32"
        a = sqrt(3.0) * r
        return (1.0 + a) * exp(-a)
    elseif cfg.kernel == "matern52"
        a = sqrt(5.0) * r
        return (1.0 + a + 5.0 * r2 / 3.0) * exp(-a)
    end
    throw(ArgumentError("Unsupported kernel $(cfg.kernel)."))
end

@inline function _kernel(cfg::GPConfig, x::NTuple{4, Float64}, y::NTuple{4, Float64}, amplitude2::Float64)::Float64
    return amplitude2 * _correlation(cfg, x, y)
end

@inline function _basis_dim(mean_basis::String)::Int
    if mean_basis == "none"
        return 0
    elseif mean_basis == "constant"
        return 1
    elseif mean_basis == "linear_alt"
        return 2
    end
    throw(ArgumentError("Unsupported mean basis $mean_basis. Use one of: $(join(MEAN_BASIS_NAMES, ", "))."))
end

@inline function _basis_row(cfg::GPConfig, alt_km::Float64, center_alt_km::Float64)::Vector{Float64}
    p = _basis_dim(cfg.mean_basis)
    p == 0 && return Float64[]
    p == 1 && return [1.0]
    return [1.0, (alt_km - center_alt_km) / cfg.alt_scale_km]
end

@inline function _safe_var(v::Vector{Float64})::Float64
    length(v) >= 2 && return var(v)
    return isempty(v) ? 0.0 : v[1]^2
end

@inline function _is_prior_scaled(cfg::GPConfig)::Bool
    cfg.amplitude_mode in AMPLITUDE_MODE_NAMES ||
        throw(ArgumentError("Unsupported amplitude mode $(cfg.amplitude_mode). Use one of: $(join(AMPLITUDE_MODE_NAMES, ", "))."))
    return cfg.amplitude_mode == "prior_scaled"
end

# One solve at a fixed amplitude, returning everything the predictor needs plus
# the (restricted) log marginal likelihood used to choose lambda.
function _solve_gp(
    rows::Vector{NTuple{4, Float64}},
    cfg::GPConfig,
    residuals::Vector{Float64},
    noise_var::Vector{Float64},
    prior_sigma::Vector{Float64},
    H::Matrix{Float64},
    amplitude2::Float64,
    prior_scaled::Bool,
)
    n = length(rows)
    p = size(H, 2)
    K = Matrix{Float64}(undef, n, n)
    @inbounds for i in 1:n
        si = prior_scaled ? prior_sigma[i] : 1.0
        K[i, i] = amplitude2 * si * si + noise_var[i] + cfg.jitter
        for j in (i + 1):n
            sj = prior_scaled ? prior_sigma[j] : 1.0
            v = amplitude2 * si * sj * _correlation(cfg, rows[i], rows[j])
            K[i, j] = v
            K[j, i] = v
        end
    end

    chol = cholesky(Symmetric(K))
    log_det_k = 2.0 * sum(log, diag(chol.U))
    if p == 0
        alpha = chol \ residuals
        loglik = -0.5 * (dot(residuals, alpha) + log_det_k)
        return (chol, Matrix{Float64}(undef, n, 0), Matrix{Float64}(undef, 0, 0), Float64[], alpha, loglik)
    end

    kinv_h = chol \ H
    a_mat = Matrix(Symmetric(H' * kinv_h))
    for j in 1:p
        a_mat[j, j] += cfg.jitter
    end
    a_chol = cholesky(Symmetric(a_mat))
    a_inv = inv(a_chol)
    beta = a_inv * (kinv_h' * residuals)
    alpha = chol \ (residuals - H * beta)
    # Restricted likelihood: the -0.5 log|H' K^-1 H| term accounts for the
    # degrees of freedom spent on the GLS mean, so lambda is not biased small.
    loglik = -0.5 * (dot(residuals - H * beta, alpha) + log_det_k + 2.0 * sum(log, diag(a_chol.U)))
    return (chol, kinv_h, a_inv, beta, alpha, loglik)
end

"""
    fit_gp_residual(points, residuals, noise_var, cfg; prior_sigma) -> FittedGP

`prior_sigma` is the GRAM ensemble relative sigma at each training point and is
required when `cfg.amplitude_mode == "prior_scaled"`, where the inflation factor
`lambda` is chosen by maximizing the restricted log marginal likelihood over a
log-spaced grid. In `stationary` mode it is ignored and the amplitude is the
variance of the residuals about the mean function.
"""
function fit_gp_residual(
    points::Vector{TrajectoryPoint},
    residuals::Vector{Float64},
    noise_var::Vector{Float64},
    cfg::GPConfig;
    prior_sigma::Vector{Float64}=Float64[],
    lambda_bounds::Union{Nothing, Tuple{Float64, Float64}}=nothing,
)::FittedGP
    n = length(points)
    (n == length(residuals) && n == length(noise_var)) || throw(ArgumentError("GP inputs must share length."))
    prior_scaled = _is_prior_scaled(cfg)
    if prior_scaled
        length(prior_sigma) == n ||
            throw(ArgumentError("amplitude_mode=prior_scaled needs one prior sigma per training point."))
        all(>(0.0), prior_sigma) || throw(ArgumentError("Prior sigmas must be positive."))
    end

    X = Matrix{Float64}(undef, n, 4)
    rows = NTuple{4, Float64}[]
    for i in 1:n
        feat = _feature_row(points[i])
        push!(rows, feat)
        X[i, 1] = feat[1]
        X[i, 2] = feat[2]
        X[i, 3] = feat[3]
        X[i, 4] = feat[4]
    end

    p = _basis_dim(cfg.mean_basis)
    center_alt_km = n == 0 ? 0.0 : mean(@view X[:, 3])
    H = Matrix{Float64}(undef, n, p)
    for i in 1:n
        h = _basis_row(cfg, X[i, 3], center_alt_km)
        for j in 1:p
            H[i, j] = h[j]
        end
    end

    amplitude2 = 0.0
    solved = nothing
    if prior_scaled
        # lambda^2 by marginal likelihood: a coarse log-spaced sweep, then a
        # finer one bracketing the best point, both inside `lambda_bounds`.
        lo, hi = lambda_bounds === nothing ?
            (p == 0 ? DEFAULT_LAMBDA_BOUNDS : DEFAULT_LAMBDA_BOUNDS_WITH_MEAN) : lambda_bounds
        (0.0 < lo <= hi) || throw(ArgumentError("lambda bounds must satisfy 0 < lo <= hi."))
        best_ll = -Inf
        best_lambda = clamp(1.0, lo, hi)
        for pass in 1:2
            grid = pass == 1 ? exp.(range(log(lo), log(hi); length=33)) :
                exp.(range(log(clamp(best_lambda / 1.6, lo, hi)), log(clamp(best_lambda * 1.6, lo, hi)); length=25))
            for lambda in grid
                candidate = _solve_gp(rows, cfg, residuals, noise_var, prior_sigma, H, lambda^2, true)
                if candidate[6] > best_ll
                    best_ll = candidate[6]
                    best_lambda = lambda
                    solved = candidate
                end
            end
        end
        amplitude2 = best_lambda^2
    else
        # The kernel amplitude must describe the residual *about* the mean
        # function, otherwise the parametric trend is double-counted as GP signal.
        detrended = p > 0 ? residuals - H * (H \ residuals) : residuals
        amplitude2 = max(_safe_var(detrended), 1.0e-6)
        solved = _solve_gp(rows, cfg, residuals, noise_var, prior_sigma, H, amplitude2, false)
    end

    chol, kinv_h, a_inv, beta, alpha, loglik = solved

    # How much of the residual scatter is above the measurement noise floor. Near
    # zero means the pass carried almost no signal, so `lambda` is unidentifiable
    # and the predictive variance should not be believed.
    detrended_final = p > 0 ? residuals - H * beta : residuals
    resid_var = _safe_var(detrended_final)
    noise_mean = isempty(noise_var) ? 0.0 : mean(noise_var)
    signal_fraction = resid_var > 0.0 ? clamp(1.0 - noise_mean / resid_var, 0.0, 1.0) : 0.0

    eff_bounds = lambda_bounds === nothing ?
        (p == 0 ? DEFAULT_LAMBDA_BOUNDS : DEFAULT_LAMBDA_BOUNDS_WITH_MEAN) : lambda_bounds
    at_bound = prior_scaled && (
        isapprox(sqrt(amplitude2), eff_bounds[1]; rtol=1.0e-6) ||
        isapprox(sqrt(amplitude2), eff_bounds[2]; rtol=1.0e-6)
    )
    return FittedGP(
        X, alpha, chol, amplitude2, cfg, H, kinv_h, a_inv, beta, center_alt_km,
        prior_scaled ? copy(prior_sigma) : Float64[], loglik, signal_fraction, at_bound,
    )
end

"""
    predict_gp_residual(model, point[, prior_sigma]) -> (mean, variance)

`prior_sigma` is the GRAM ensemble relative sigma at `point` and is required in
`prior_scaled` mode, where the kernel is `lambda^2 s(x) s(x") rho(x, x")`. The
returned variance is the full predictive variance of the residual: in
`prior_scaled` mode it tends to `lambda^2 * prior_sigma^2` far from data, so no
further combination with the prior is needed or correct.
"""
function predict_gp_residual(
    model::FittedGP, point::TrajectoryPoint, prior_sigma::Float64=NaN
)::Tuple{Float64, Float64}
    prior_scaled = model.cfg.amplitude_mode == "prior_scaled"
    if prior_scaled
        isfinite(prior_sigma) && prior_sigma > 0.0 ||
            throw(ArgumentError("amplitude_mode=prior_scaled needs a positive prior sigma at the prediction point."))
    end
    s_x = prior_scaled ? prior_sigma : 1.0

    x = _feature_row(point)
    n = size(model.X, 1)
    k = Vector{Float64}(undef, n)
    for i in 1:n
        xi = (model.X[i, 1], model.X[i, 2], model.X[i, 3], model.X[i, 4])
        s_i = prior_scaled ? model.prior_sigma[i] : 1.0
        k[i] = model.amplitude2 * s_x * s_i * _correlation(model.cfg, xi, x)
    end
    v = model.chol.L \ k
    kxx = model.amplitude2 * s_x * s_x + model.cfg.jitter
    mu = dot(k, model.alpha)
    sigma2 = kxx - dot(v, v)

    if !isempty(model.beta)
        h = _basis_row(model.cfg, x[3], model.center_alt_km)
        mu += dot(h, model.beta)
        # Universal-kriging variance: uncertainty in the GLS mean coefficients,
        # which grows where the basis is evaluated away from the training span.
        r = h - model.Kinv_H' * k
        sigma2 += dot(r, model.A_inv * r)
    end
    return mu, max(sigma2, 1.0e-12)
end

"""
    mean_function_report(model) -> NamedTuple

The fitted parametric mean, for reporting what the aerocapture pass concluded
about the bulk state of the atmosphere independently of the GP's local
corrections.
"""
function mean_function_report(model::FittedGP)
    return (
        mean_basis=model.cfg.mean_basis,
        beta_constant=isempty(model.beta) ? 0.0 : model.beta[1],
        beta_alt_slope=length(model.beta) >= 2 ? model.beta[2] : 0.0,
        center_alt_km=model.center_alt_km,
        amplitude2=model.amplitude2,
        amplitude_mode=model.cfg.amplitude_mode,
        lambda=model.cfg.amplitude_mode == "prior_scaled" ? sqrt(model.amplitude2) : NaN,
        log_marginal=model.log_marginal,
        signal_fraction=model.signal_fraction,
        lambda_at_bound=model.lambda_at_bound,
    )
end

function predict_density_profile(
    par::Parameterization,
    gp::FittedGP,
    pred_points::Vector{TrajectoryPoint},
    pred_gram::Vector{GramSample}
)
    prior_scaled = gp.cfg.amplitude_mode == "prior_scaled"
    pred_mean = Vector{Float64}(undef, length(pred_points))
    pred_std = Vector{Float64}(undef, length(pred_points))
    for i in eachindex(pred_points)
        prior_mu, prior_var = _prior_target(par, pred_gram[i])
        if prior_scaled
            # The GRAM prior variance is already inside the kernel, so the GP
            # posterior variance is the answer. Combining it with `prior_var`
            # again would double-count.
            delta_mu, delta_var = predict_gp_residual(gp, pred_points[i], sqrt(prior_var))
            fused_mu = prior_mu + delta_mu
            fused_var = delta_var
        else
            delta_mu, delta_var = predict_gp_residual(gp, pred_points[i])
            fused_mu = prior_mu + delta_mu
            # Legacy heuristic: an additive mean combined with an
            # inverse-variance fusion. The two come from different models, and
            # the result is capped at the kernel amplitude everywhere. Retained
            # only to reproduce the published stationary results.
            fused_var = 1.0 / (1.0 / max(prior_var, 1.0e-12) + 1.0 / max(delta_var, 1.0e-12))
        end
        pred_mean[i], pred_std[i] = _physical_moments(par, fused_mu, fused_var, pred_gram[i])
    end
    return pred_mean, pred_std
end
