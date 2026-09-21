const KS_I4 = Matrix{Float64}(I, 4, 4)

const KS_DLP = let
    m1 = Float64[
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 -1
    ]
    m2 = Float64[
        0 -1 0 0
        1 0 0 0
        0 0 0 1
        0 0 1 0
    ]
    m3 = Float64[
        0 0 -1 0
        0 0 0 -1
        1 0 0 0
        0 -1 0 0
    ]
    m4 = Float64[
        0 0 0 1
        0 0 -1 0
        0 1 0 0
        1 0 0 0
    ]
    (m1, m2, m3, m4)
end

const KS_DSP = Tuple(transpose(m) for m in KS_DLP)

function ks_lambda_matrix(p::AbstractVector{<:Real})
    return Float64[
        p[1] -p[2] -p[3] p[4]
        p[2] p[1] -p[4] -p[3]
        p[3] p[4] p[1] p[2]
    ]
end

function ks_l_matrix(p::AbstractVector{<:Real})
    return Float64[
        p[1] -p[2] -p[3] p[4]
        p[2] p[1] -p[4] -p[3]
        p[3] p[4] p[1] p[2]
        p[4] -p[3] p[2] -p[1]
    ]
end

ks_s_matrix(p::AbstractVector{<:Real}) = transpose(ks_l_matrix(p))
ks_radius_scalar(p::AbstractVector{<:Real}) = dot(p, p)

function ks_position(p::AbstractVector{<:Real})
    lambda_p = ks_lambda_matrix(p)
    return lambda_p * Float64.(p)
end

function ks_velocity(p::AbstractVector{<:Real}, q::AbstractVector{<:Real})
    rho = ks_radius_scalar(p)
    rho > 0.0 || throw(ArgumentError("KS radius must be positive."))
    return (2.0 / rho) .* (ks_lambda_matrix(p) * Float64.(q))
end

function ks_augmented_acceleration(accel::AbstractVector{<:Real})
    accel4 = zeros(Float64, 4)
    accel4[1:3] .= Float64.(accel)
    return accel4
end

function ks_satellite_count(u::AbstractVector{<:Real})
    length(u) % 10 == 0 || throw(ArgumentError("Packed KS state length must be divisible by 10."))
    return length(u) ÷ 10
end

function ks_state_slice(sat_id::Int)
    base = 10 * (sat_id - 1)
    return (base + 1):(base + 10)
end

function ks_state_components(u::AbstractVector{<:Real}, sat_id::Int)
    base = 10 * (sat_id - 1)
    p = Float64.(u[(base + 1):(base + 4)])
    q = Float64.(u[(base + 5):(base + 8)])
    h = float(u[base + 9])
    t = float(u[base + 10])
    return p, q, h, t
end

function ks_cartesian_components(u::AbstractVector{<:Real}, sat_id::Int)
    p, q, h, t = ks_state_components(u, sat_id)
    r = ks_position(p)
    v = ks_velocity(p, q)
    return r, v, h, t
end

function ks_satellite_times(u::AbstractVector{<:Real})
    nsat = ks_satellite_count(u)
    times = zeros(Float64, nsat)
    for sat in 1:nsat
        _, _, _, t_sat = ks_state_components(u, sat)
        times[sat] = t_sat
    end
    return times
end

ks_mean_physical_time(u::AbstractVector{<:Real}) = sum(ks_satellite_times(u)) / ks_satellite_count(u)
ks_min_physical_time(u::AbstractVector{<:Real}) = minimum(ks_satellite_times(u))

function cartesian_to_ks_state(
    r::AbstractVector{<:Real},
    v::AbstractVector{<:Real};
    mu::Real = MU_EARTH_KM3_S2,
    t0::Real = 0.0,
)
    x = float(r[1])
    y = float(r[2])
    z = float(r[3])
    rmag = norm(Float64.(r))
    rmag > 0.0 || throw(ArgumentError("Cartesian position magnitude must be positive."))

    p = zeros(Float64, 4)
    if x > -0.999999999999 * rmag
        p[1] = sqrt(max((rmag + x) / 2.0, 0.0))
        denom = 2.0 * max(p[1], eps())
        p[2] = y / denom
        p[3] = z / denom
        p[4] = 0.0
    else
        p[2] = sqrt(max((rmag - x) / 2.0, 0.0))
        denom = 2.0 * max(p[2], eps())
        p[1] = y / denom
        p[3] = 0.0
        p[4] = z / denom
    end

    lambda_p = ks_lambda_matrix(p)
    q = 0.5 .* (transpose(lambda_p) * Float64.(v))

    energy = 0.5 * dot(Float64.(v), Float64.(v)) - float(mu) / rmag
    h = -2.0 * energy

    x_ks = zeros(Float64, 10)
    x_ks[1:4] .= p
    x_ks[5:8] .= q
    x_ks[9] = h
    x_ks[10] = float(t0)
    return x_ks
end

function ks_pair_indices(nsat::Int)
    nsat >= 1 || throw(ArgumentError("nsat must be positive."))
    pairs = Tuple{Int, Int}[]
    for i in 1:(nsat - 1)
        for j in (i + 1):nsat
            push!(pairs, (i, j))
        end
    end
    return pairs
end
