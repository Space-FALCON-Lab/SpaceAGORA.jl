function specific_energy(df::DataFrame, mu::Float64)
    r = sqrt.(df.sc1_pos_1.^2 .+ df.sc1_pos_2.^2 .+ df.sc1_pos_3.^2)
    v2 = df.sc1_vel_1.^2 .+ df.sc1_vel_2.^2 .+ df.sc1_vel_3.^2
    return 0.5 .* v2 .- mu ./ r
end

function specific_angular_momentum_magnitude(df::DataFrame)
    h1 = df.sc1_pos_2 .* df.sc1_vel_3 .- df.sc1_pos_3 .* df.sc1_vel_2
    h2 = df.sc1_pos_3 .* df.sc1_vel_1 .- df.sc1_pos_1 .* df.sc1_vel_3
    h3 = df.sc1_pos_1 .* df.sc1_vel_2 .- df.sc1_pos_2 .* df.sc1_vel_1
    return sqrt.(h1.^2 .+ h2.^2 .+ h3.^2)
end

function unwrap_angle_series(values::AbstractVector{<:Real})
    n = length(values)
    n == 0 && return Float64[]
    out = Vector{Float64}(undef, n)
    out[1] = Float64(values[1])
    @inbounds for i in 2:n
        prev = Float64(values[i - 1])
        cur = Float64(values[i])
        Δ = atan(sin(cur - prev), cos(cur - prev))
        out[i] = out[i - 1] + Δ
    end
    return out
end

function linear_regression_slope(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
    length(xs) == length(ys) || throw(ArgumentError("linear_regression_slope requires vectors of equal length."))
    length(xs) >= 2 || throw(ArgumentError("linear_regression_slope requires at least two samples."))
    x̄ = sum(xs) / length(xs)
    ȳ = sum(ys) / length(ys)
    num = 0.0
    den = 0.0
    @inbounds for i in eachindex(xs, ys)
        dx = Float64(xs[i]) - x̄
        dy = Float64(ys[i]) - ȳ
        num += dx * dy
        den += dx * dx
    end
    den > 0.0 || throw(ArgumentError("linear_regression_slope requires non-degenerate x samples."))
    return num / den
end

function harmonics_full_normalization_factor(l::Int, m::Int)::Float64
    l >= 0 || throw(ArgumentError("Degree must be >= 0, got $l."))
    0 <= m <= l || throw(ArgumentError("Order must satisfy 0 <= m <= l, got l=$l, m=$m."))
    δ0m = (m == 0) ? 1.0 : 0.0
    ratio = exp(0.5 * (loggamma(l - m + 1) - loggamma(l + m + 1)))
    return sqrt((2.0 - δ0m) * (2.0 * l + 1.0)) * ratio
end

function write_gravity_harmonics_csv(
    coeff_fn::Function,
    path::AbstractString,
    max_degree::Int,
    max_order::Int
)
    open(path, "w") do io
        println(io, "degree,order,C,S")
        for l in 0:max_degree
            for m in 0:min(max_order, l)
                C, S = coeff_fn(l, m)
                println(io, string(l, ",", m, ",", C, ",", S))
            end
        end
    end
    return path
end

function write_icgem_from_harmonics_model(
    path::AbstractString,
    model::SimulationModel.DynamicEffectors.GravitationalHarmonicsModel;
    model_name::AbstractString="SpaceAGORA_Test"
)
    open(path, "w") do io
        println(io, "product_type gravity_field")
        println(io, "modelname ", model_name)
        println(io, "earth_gravity_constant ", model.planet.μ)
        println(io, "radius ", model.planet.Rp_e)
        println(io, "max_degree ", model.L)
        println(io, "errors no")
        println(io, "norm fully_normalized")
        println(io, "end_of_head")
        for l in 0:model.L
            for m in 0:min(model.M, l)
                i = l + 1
                j = m + 1
                println(io, "gfc ", l, " ", m, " ", model.C[i, j], " ", model.S[i, j])
            end
        end
    end
    return path
end
