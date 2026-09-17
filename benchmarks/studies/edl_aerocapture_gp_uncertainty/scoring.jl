struct MetricSummary
    rmse::Float64
    weighted_rmse::Float64
    rmse_log::Float64
    weighted_rmse_log::Float64
    nlpd::Float64
    weighted_nlpd::Float64
    coverage_1sigma::Float64
    coverage_2sigma::Float64
    weighted_coverage_1sigma::Float64
    weighted_coverage_2sigma::Float64
end

@inline function _gaussian_nlpd(mu::Float64, sigma::Float64, y::Float64)::Float64
    sigma2 = max(sigma * sigma, 1.0e-18)
    return 0.5 * log(2.0 * π * sigma2) + 0.5 * (y - mu)^2 / sigma2
end

@inline function _weighted_mean(values::Vector{Float64}, weights::Vector{Float64})::Float64
    return dot(values, weights) / max(sum(weights), 1.0e-18)
end

function score_predictions(
    truth::Vector{Float64},
    pred_mean::Vector{Float64},
    pred_std::Vector{Float64},
    weights::Vector{Float64}
)::MetricSummary
    n = length(truth)
    (n == length(pred_mean) && n == length(pred_std) && n == length(weights)) || throw(ArgumentError("Prediction vectors must share length."))

    sq_err = Vector{Float64}(undef, n)
    sq_err_log = Vector{Float64}(undef, n)
    nlpd = Vector{Float64}(undef, n)
    c1 = Vector{Float64}(undef, n)
    c2 = Vector{Float64}(undef, n)
    for i in 1:n
        sq_err[i] = (pred_mean[i] - truth[i])^2
        # Absolute-density error is dominated by the densest altitudes: on the
        # default corridor the 10-20 km band carries 87% of the weighted squared
        # error while the q-max band near 35 km carries 0.8%. The log-density
        # error weights every altitude comparably and is the metric that can
        # actually resolve an improvement in the sensed band.
        sq_err_log[i] = (log(max(pred_mean[i], 1.0e-18)) - log(max(truth[i], 1.0e-18)))^2
        nlpd[i] = _gaussian_nlpd(pred_mean[i], pred_std[i], truth[i])
        c1[i] = abs(truth[i] - pred_mean[i]) <= pred_std[i] ? 1.0 : 0.0
        c2[i] = abs(truth[i] - pred_mean[i]) <= 2.0 * pred_std[i] ? 1.0 : 0.0
    end

    return MetricSummary(
        sqrt(mean(sq_err)),
        sqrt(_weighted_mean(sq_err, weights)),
        sqrt(mean(sq_err_log)),
        sqrt(_weighted_mean(sq_err_log, weights)),
        mean(nlpd),
        _weighted_mean(nlpd, weights),
        mean(c1),
        mean(c2),
        _weighted_mean(c1, weights),
        _weighted_mean(c2, weights)
    )
end
