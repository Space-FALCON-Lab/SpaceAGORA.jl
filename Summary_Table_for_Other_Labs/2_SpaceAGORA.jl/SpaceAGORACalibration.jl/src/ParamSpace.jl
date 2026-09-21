module ParamSpace

using Random
using ..Spec: CalibrationSpec, ParameterSpec, continuous, integer, categorical

export Candidate
export sample_random_candidate, sample_lhs_population, sample_initial_population, perturb_candidate
export candidate_signature

Base.@kwdef struct Candidate
    id::Int
    values::Dict{String, Any}
    stage::String = "global_search_quick"
end

@inline function _sample_value(rng::AbstractRNG, p::ParameterSpec)
    if p.kind == continuous
        return rand(rng) * (p.upper - p.lower) + p.lower
    elseif p.kind == integer
        lo = ceil(Int, p.lower)
        hi = floor(Int, p.upper)
        lo <= hi || throw(ArgumentError("integer parameter $(p.name) has empty integer range."))
        return rand(rng, lo:hi)
    elseif p.kind == categorical
        return rand(rng, p.choices)
    end
    throw(ArgumentError("Unsupported parameter kind for $(p.name)."))
end

function sample_random_candidate(
    rng::AbstractRNG,
    spec::CalibrationSpec,
    id::Int;
    stage::String="global_search_quick"
)::Candidate
    values = Dict{String, Any}()
    for p in spec.parameters
        values[p.name] = _sample_value(rng, p)
    end
    return Candidate(id=id, values=values, stage=stage)
end

function _lhs_unit_matrix(rng::AbstractRNG, n::Int, dims::Int)
    m = Matrix{Float64}(undef, n, dims)
    for j in 1:dims
        perm = randperm(rng, n)
        for i in 1:n
            m[i, j] = (perm[i] - rand(rng)) / n
        end
    end
    return m
end

@inline function _value_from_unit(p::ParameterSpec, u::Float64)
    if p.kind == continuous
        return p.lower + clamp(u, 0.0, 1.0) * (p.upper - p.lower)
    elseif p.kind == integer
        lo = ceil(Int, p.lower)
        hi = floor(Int, p.upper)
        lo <= hi || throw(ArgumentError("integer parameter $(p.name) has empty integer range."))
        x = lo + floor(Int, clamp(u, 0.0, 1.0) * (hi - lo + 1))
        return clamp(x, lo, hi)
    elseif p.kind == categorical
        n = length(p.choices)
        n > 0 || throw(ArgumentError("categorical parameter $(p.name) has empty choices."))
        idx = clamp(1 + floor(Int, clamp(u, 0.0, 1.0 - eps(Float64)) * n), 1, n)
        return p.choices[idx]
    end
    throw(ArgumentError("Unsupported parameter kind for $(p.name)."))
end

function sample_lhs_population(
    rng::AbstractRNG,
    spec::CalibrationSpec,
    n::Int;
    start_id::Int=1,
    stage::String="global_search_quick"
)::Vector{Candidate}
    n > 0 || return Candidate[]
    dims = length(spec.parameters)
    dims > 0 || throw(ArgumentError("spec.parameters cannot be empty."))

    units = _lhs_unit_matrix(rng, n, dims)
    out = Vector{Candidate}(undef, n)
    next_id = start_id

    for i in 1:n
        values = Dict{String, Any}()
        for j in 1:dims
            p = spec.parameters[j]
            values[p.name] = _value_from_unit(p, units[i, j])
        end
        out[i] = Candidate(id=next_id, values=values, stage=stage)
        next_id += 1
    end

    return out
end

function sample_initial_population(
    rng::AbstractRNG,
    spec::CalibrationSpec,
    n::Int;
    start_id::Int=1,
    stage::String="global_search_quick",
    design::String="random"
)::Vector{Candidate}
    n > 0 || return Candidate[]

    if lowercase(design) == "lhs"
        return sample_lhs_population(rng, spec, n; start_id=start_id, stage=stage)
    end

    out = Candidate[]
    seen = Set{String}()
    next_id = start_id
    max_attempts = max(n * 20, 64)
    attempts = 0

    while length(out) < n && attempts < max_attempts
        cand = sample_random_candidate(rng, spec, next_id; stage=stage)
        sig = candidate_signature(cand)
        if sig ∉ seen
            push!(seen, sig)
            push!(out, cand)
            next_id += 1
        else
            next_id += 1
        end
        attempts += 1
    end

    while length(out) < n
        cand = sample_random_candidate(rng, spec, next_id; stage=stage)
        push!(out, cand)
        next_id += 1
    end

    return out
end

function perturb_candidate(
    rng::AbstractRNG,
    spec::CalibrationSpec,
    base::Candidate,
    id::Int;
    stage::String="local_refine_full",
    scale::Float64=0.12,
    perturb_discrete::Bool=true
)::Candidate
    values = Dict{String, Any}()

    for p in spec.parameters
        v = base.values[p.name]
        if p.kind == continuous
            span = p.upper - p.lower
            delta = (2.0 * rand(rng) - 1.0) * scale * span
            values[p.name] = clamp(Float64(v) + delta, p.lower, p.upper)
        elseif p.kind == integer && perturb_discrete
            span = max(1, floor(Int, p.upper) - ceil(Int, p.lower) + 1)
            max_delta = max(1, round(Int, scale * span))
            delta = rand(rng, -max_delta:max_delta)
            lo = ceil(Int, p.lower)
            hi = floor(Int, p.upper)
            values[p.name] = clamp(Int(v) + delta, lo, hi)
        elseif p.kind == categorical && perturb_discrete
            if rand(rng) < 0.25
                values[p.name] = rand(rng, p.choices)
            else
                values[p.name] = v
            end
        else
            values[p.name] = v
        end
    end

    return Candidate(id=id, values=values, stage=stage)
end

function candidate_signature(c::Candidate)::String
    keys_sorted = sort(collect(keys(c.values)))
    chunks = String[]
    for k in keys_sorted
        push!(chunks, string(k, "=", c.values[k]))
    end
    return join(chunks, "|")
end

end # module ParamSpace
