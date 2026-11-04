function _parameter_space(cfg::TunerConfig)::Vector{ParameterSpec}
    return ParameterSpec[
        ParameterSpec(name="epoch_shift_s", kind=:integer, lower=-180.0, upper=180.0),
        ParameterSpec(name="ra_scale", kind=:continuous, lower=0.995, upper=1.005),
        ParameterSpec(name="rp_altitude_offset_m", kind=:continuous, lower=-3000.0, upper=3000.0),
        ParameterSpec(name="i_offset_deg", kind=:continuous, lower=-0.03, upper=0.03),
        ParameterSpec(name="aop_offset_deg", kind=:continuous, lower=-0.2, upper=0.2),
        ParameterSpec(name="raan_offset_deg", kind=:continuous, lower=-0.2, upper=0.2),
        ParameterSpec(name="ta_offset_deg", kind=:continuous, lower=-0.5, upper=0.5),
        ParameterSpec(name="bus_mass_scale", kind=:continuous, lower=0.99, upper=1.01),
        ParameterSpec(name="prop_mass_scale", kind=:continuous, lower=0.95, upper=1.05),
        ParameterSpec(name="panel_mass_scale", kind=:continuous, lower=0.95, upper=1.05),
        ParameterSpec(name="bus_dims_scale", kind=:continuous, lower=0.995, upper=1.005),
        ParameterSpec(name="panel_dims_scale", kind=:continuous, lower=0.995, upper=1.005),
        ParameterSpec(name="panel_offset_scale", kind=:continuous, lower=0.99, upper=1.01),
        ParameterSpec(name="srp_cr_scale", kind=:continuous, lower=0.95, upper=1.05),
        ParameterSpec(name="srp_area_scale", kind=:continuous, lower=0.95, upper=1.05),
        ParameterSpec(name="gravity_variant", kind=:categorical, choices=copy(GRAVITY_VARIANTS)),
        ParameterSpec(name="dv_global_scale", kind=:continuous, lower=0.8, upper=1.2),
        ParameterSpec(name="dv_early_scale", kind=:continuous, lower=0.8, upper=1.2),
        ParameterSpec(name="dv_orbit7_bias_mps", kind=:continuous, lower=-0.3, upper=0.3),
        ParameterSpec(name="solver_mode", kind=:categorical, choices=copy(cfg.solver_modes)),
        ParameterSpec(name="dt_max_orbit_s", kind=:categorical, choices=string.(cfg.dt_max_orbit_values)),
        ParameterSpec(name="dt_max_atm_s", kind=:categorical, choices=string.(cfg.dt_max_atm_values))
    ]
end

@inline function _sample_value(rng::AbstractRNG, p::ParameterSpec)
    if p.kind == :continuous
        return rand(rng) * (p.upper - p.lower) + p.lower
    elseif p.kind == :integer
        lo = ceil(Int, p.lower)
        hi = floor(Int, p.upper)
        lo <= hi || throw(ArgumentError("integer parameter $(p.name) has empty range."))
        return rand(rng, lo:hi)
    elseif p.kind == :categorical
        return rand(rng, p.choices)
    end
    throw(ArgumentError("Unsupported parameter kind $(p.kind)."))
end

function _candidate_signature(c::TuneCandidate)::String
    chunks = String[]
    for key in PARAMETER_ORDER
        push!(chunks, string(key, "=", c.values[key]))
    end
    return join(chunks, "|")
end

function _candidate_row_nt(c::TuneCandidate)
    return (
        candidate_id=c.id,
        epoch_shift_s=Int(c.values["epoch_shift_s"]),
        ra_scale=Float64(c.values["ra_scale"]),
        rp_altitude_offset_m=Float64(c.values["rp_altitude_offset_m"]),
        i_offset_deg=Float64(c.values["i_offset_deg"]),
        aop_offset_deg=Float64(c.values["aop_offset_deg"]),
        raan_offset_deg=Float64(c.values["raan_offset_deg"]),
        ta_offset_deg=Float64(c.values["ta_offset_deg"]),
        bus_mass_scale=Float64(c.values["bus_mass_scale"]),
        prop_mass_scale=Float64(c.values["prop_mass_scale"]),
        panel_mass_scale=Float64(c.values["panel_mass_scale"]),
        bus_dims_scale=Float64(c.values["bus_dims_scale"]),
        panel_dims_scale=Float64(c.values["panel_dims_scale"]),
        panel_offset_scale=Float64(c.values["panel_offset_scale"]),
        srp_cr_scale=Float64(c.values["srp_cr_scale"]),
        srp_area_scale=Float64(c.values["srp_area_scale"]),
        gravity_variant=String(c.values["gravity_variant"]),
        dv_global_scale=Float64(c.values["dv_global_scale"]),
        dv_early_scale=Float64(c.values["dv_early_scale"]),
        dv_orbit7_bias_mps=Float64(c.values["dv_orbit7_bias_mps"]),
        solver_mode_requested=String(c.values["solver_mode"]),
        dt_max_orbit_requested_s=Float64(c.values["dt_max_orbit_s"]),
        dt_max_atm_requested_s=Float64(c.values["dt_max_atm_s"])
    )
end

@inline function _candidate_vector(params::Vector{ParameterSpec}, c::TuneCandidate)::Vector{Float64}
    vec = Vector{Float64}(undef, length(params))
    for i in eachindex(params)
        p = params[i]
        v = c.values[p.name]
        if p.kind == :continuous || p.kind == :integer
            span = p.upper - p.lower
            vec[i] = span > 0 ? clamp((Float64(v) - p.lower) / span, 0.0, 1.0) : 0.5
        else
            idx = findfirst(==(String(v)), p.choices)
            if idx === nothing || length(p.choices) <= 1
                vec[i] = 0.0
            else
                vec[i] = (idx - 1) / (length(p.choices) - 1)
            end
        end
    end
    return vec
end

function _candidate_from_vector(
    params::Vector{ParameterSpec},
    u::Vector{Float64},
    id::Int,
    stage::String
)::TuneCandidate
    values = Dict{String, Any}()
    for i in eachindex(params)
        p = params[i]
        ui = clamp(u[i], 0.0, 1.0)
        if p.kind == :continuous
            values[p.name] = p.lower + ui * (p.upper - p.lower)
        elseif p.kind == :integer
            lo = ceil(Int, p.lower)
            hi = floor(Int, p.upper)
            x = lo + floor(Int, ui * (hi - lo + 1))
            values[p.name] = clamp(x, lo, hi)
        else
            idx = clamp(1 + floor(Int, ui * length(p.choices)), 1, length(p.choices))
            values[p.name] = p.choices[idx]
        end
    end

    values["dt_max_orbit_s"] = values["dt_max_orbit_s"] isa AbstractString ? parse(Float64, values["dt_max_orbit_s"]) : Float64(values["dt_max_orbit_s"])
    values["dt_max_atm_s"] = values["dt_max_atm_s"] isa AbstractString ? parse(Float64, values["dt_max_atm_s"]) : Float64(values["dt_max_atm_s"])
    values["epoch_shift_s"] = Int(values["epoch_shift_s"])
    return TuneCandidate(id=id, values=values, stage=stage)
end

function _lhs_unit_matrix(rng::AbstractRNG, n::Int, d::Int)::Matrix{Float64}
    m = Matrix{Float64}(undef, n, d)
    for j in 1:d
        perm = randperm(rng, n)
        for i in 1:n
            m[i, j] = (perm[i] - rand(rng)) / n
        end
    end
    return m
end

function _initial_design(
    cfg::TunerConfig,
    params::Vector{ParameterSpec},
    n::Int;
    start_id::Int=1,
    stage::String="quick_global_init"
)::Vector{TuneCandidate}
    rng = MersenneTwister(hash((cfg.seed, stage, "initial_design")))
    out = TuneCandidate[]
    seen = Set{String}()

    if cfg.init_design == :lhs
        u = _lhs_unit_matrix(rng, n, length(params))
        next_id = start_id
        for i in 1:n
            cand = _candidate_from_vector(params, vec(u[i, :]), next_id, stage)
            sig = _candidate_signature(cand)
            if sig in seen
                continue
            end
            push!(seen, sig)
            push!(out, cand)
            next_id += 1
        end
        while length(out) < n
            vals = Dict{String, Any}()
            for p in params
                vals[p.name] = _sample_value(rng, p)
            end
            vals["dt_max_orbit_s"] = Float64(vals["dt_max_orbit_s"])
            vals["dt_max_atm_s"] = Float64(vals["dt_max_atm_s"])
            vals["epoch_shift_s"] = Int(vals["epoch_shift_s"])
            cand = TuneCandidate(id=start_id + length(out), values=vals, stage=stage)
            sig = _candidate_signature(cand)
            if sig in seen
                continue
            end
            push!(seen, sig)
            push!(out, cand)
        end
        return out
    end

    next_id = start_id
    while length(out) < n
        vals = Dict{String, Any}()
        for p in params
            vals[p.name] = _sample_value(rng, p)
        end
        vals["dt_max_orbit_s"] = Float64(vals["dt_max_orbit_s"])
        vals["dt_max_atm_s"] = Float64(vals["dt_max_atm_s"])
        vals["epoch_shift_s"] = Int(vals["epoch_shift_s"])
        cand = TuneCandidate(id=next_id, values=vals, stage=stage)
        sig = _candidate_signature(cand)
        if sig in seen
            next_id += 1
            continue
        end
        push!(seen, sig)
        push!(out, cand)
        next_id += 1
    end
    return out
end

@inline function _rbf_kernel(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, l::Float64)::Float64
    d2 = 0.0
    @inbounds for i in eachindex(x)
        d = x[i] - y[i]
        d2 += d * d
    end
    return exp(-d2 / (2.0 * l * l))
end

@inline function _phi(z::Float64)::Float64
    return exp(-0.5 * z * z) / sqrt(2.0 * π)
end

@inline function _Phi(z::Float64)::Float64
    return 0.5 * (1.0 + tanh(sqrt(2.0 / π) * (z + 0.044715 * z^3)))
end

function _fit_surrogate(
    cfg::TunerConfig,
    params::Vector{ParameterSpec},
    candidates::Vector{TuneCandidate},
    y::Vector{Float64}
)::Union{Nothing, SurrogateModel}
    n = length(candidates)
    n == 0 && return nothing

    X = Matrix{Float64}(undef, n, length(params))
    for i in 1:n
        X[i, :] = _candidate_vector(params, candidates[i])
    end

    l = cfg.surrogate_length_scale
    noise = max(cfg.surrogate_noise, 1.0e-9)

    K = Matrix{Float64}(undef, n, n)
    for i in 1:n
        xi = @view X[i, :]
        for j in i:n
            xj = @view X[j, :]
            kval = _rbf_kernel(xi, xj, l)
            K[i, j] = kval
            K[j, i] = kval
        end
    end

    A = K + noise * I
    Ainv = Matrix(inv(Matrix(A)))
    alpha = Ainv * y

    diagA = diag(Ainv)
    loo_pred = similar(y)
    for i in eachindex(y)
        if abs(diagA[i]) < 1.0e-12
            loo_pred[i] = y[i]
        else
            loo_pred[i] = y[i] - alpha[i] / diagA[i]
        end
    end
    residual_sigma = length(y) > 2 ? std(y .- loo_pred) : (length(y) == 1 ? 0.0 : std(y))

    return SurrogateModel(
        X=X,
        y=y,
        Ainv=Ainv,
        alpha=alpha,
        length_scale=l,
        residual_sigma=max(residual_sigma, 1.0e-9)
    )
end

function _predict(model::SurrogateModel, x::Vector{Float64})::Tuple{Float64, Float64}
    n = size(model.X, 1)
    k = Vector{Float64}(undef, n)
    dmin2 = Inf
    for i in 1:n
        xi = @view model.X[i, :]
        k[i] = _rbf_kernel(x, xi, model.length_scale)
        d2 = sum((x .- xi) .^ 2)
        d2 < dmin2 && (dmin2 = d2)
    end

    mu = dot(k, model.alpha)
    var = 1.0 - dot(k, model.Ainv * k)
    base_sigma = sqrt(max(var, 1.0e-12))
    density = exp(-dmin2 / (2.0 * model.length_scale * model.length_scale))
    sigma = base_sigma + model.residual_sigma * (1.0 - density)
    return mu, max(sigma, 1.0e-12)
end

@inline function _ei(mu::Float64, sigma::Float64, best::Float64, xi::Float64)::Float64
    sigma <= 1.0e-12 && return 0.0
    improve = best - mu - xi
    z = improve / sigma
    return improve * _Phi(z) + sigma * _phi(z)
end

@inline function _kappa_schedule(cfg::TunerConfig, iter::Int)::Float64
    cfg.n_bo_iters <= 1 && return cfg.kappa_end
    t = (iter - 1) / max(cfg.n_bo_iters - 1, 1)
    return cfg.kappa_start + (cfg.kappa_end - cfg.kappa_start) * t
end

function _propose_batch(
    cfg::TunerConfig,
    params::Vector{ParameterSpec},
    model::SurrogateModel,
    observed_scores::Vector{Float64},
    seen::Set{String},
    next_id::Int,
    iter::Int,
    q::Int
)::Vector{TuneCandidate}
    proposed = TuneCandidate[]
    local_seen = Set{String}()
    best_obs = minimum(observed_scores)
    kappa = _kappa_schedule(cfg, iter)

    for j in 1:q
        rng = MersenneTwister(hash((cfg.seed, "bo_pool", iter, j, next_id)))
        pool = TuneCandidate[]
        attempts = 0
        max_attempts = max(cfg.bo_pool_size * 40, 1024)
        while length(pool) < cfg.bo_pool_size && attempts < max_attempts
            vals = Dict{String, Any}()
            for p in params
                vals[p.name] = _sample_value(rng, p)
            end
            vals["dt_max_orbit_s"] = Float64(vals["dt_max_orbit_s"])
            vals["dt_max_atm_s"] = Float64(vals["dt_max_atm_s"])
            vals["epoch_shift_s"] = Int(vals["epoch_shift_s"])
            cand = TuneCandidate(id=next_id + j - 1, values=vals, stage="quick_global_bo")
            sig = _candidate_signature(cand)
            if sig in seen || sig in local_seen
                attempts += 1
                continue
            end
            push!(pool, cand)
            push!(local_seen, sig)
            attempts += 1
        end
        isempty(pool) && break

        best_idx = 1
        best_val = cfg.acquisition == :ei ? -Inf : Inf
        best_sig = _candidate_signature(pool[1])

        for i in eachindex(pool)
            cand = pool[i]
            x = _candidate_vector(params, cand)
            mu, sigma = _predict(model, x)
            acq_val = if cfg.acquisition == :ei
                _ei(mu, sigma, best_obs, cfg.ei_xi)
            else
                mu - kappa * sigma
            end

            sig = _candidate_signature(cand)
            if cfg.acquisition == :ei
                if acq_val > best_val || (acq_val == best_val && sig < best_sig)
                    best_idx = i
                    best_val = acq_val
                    best_sig = sig
                end
            else
                if acq_val < best_val || (acq_val == best_val && sig < best_sig)
                    best_idx = i
                    best_val = acq_val
                    best_sig = sig
                end
            end
        end

        chosen = pool[best_idx]
        push!(proposed, chosen)
        push!(seen, _candidate_signature(chosen))
    end

    return proposed
end
