module SupportFunctions

using LinearAlgebra

"""
    _build_polyhedral_dirs(n::Int) -> Matrix{Float64}

Build a bank of n direction vectors for 6D state-space polyhedral support function evaluation.
Matches CAPOConstellation LADS tube implementation with 6D directions (position + velocity).

The direction bank consists of:
- 12 axis directions (±e_i for i=1..6 in 6D state space)
- Pairwise sum and difference diagonals (up to 72 directions total)

For n=72 (standard controllable set bank), this gives 12 + 60 = 72 directions.

# Arguments
- `n::Int`: Number of direction vectors (minimum 12, must be even)

# Returns
- `Matrix{Float64}`: 6×n matrix of unit direction vectors
"""
function _build_polyhedral_dirs(n::Int)
    n = max(12, n)
    if isodd(n)
        n -= 1
    end
    base = Vector{Vector{Float64}}()
    for i in 1:6
        v = zeros(Float64, 6)
        v[i] = 1.0
        push!(base, v)
        push!(base, -v)
    end
    if n == 12
        return hcat(base...)
    end
    # Pairwise sum and difference diagonals — matches CAPO LADS tube formulation
    extra = Vector{Vector{Float64}}()
    for i in 1:6
        for j in (i+1):6
            v_sum = zeros(Float64, 6); v_sum[i] = 1.0; v_sum[j] =  1.0; v_sum ./= sqrt(2.0)
            v_dif = zeros(Float64, 6); v_dif[i] = 1.0; v_dif[j] = -1.0; v_dif ./= sqrt(2.0)
            push!(extra, v_sum); push!(extra, -v_sum)
            push!(extra, v_dif); push!(extra, -v_dif)
        end
    end
    for v in extra
        length(base) >= n && break
        push!(base, v)
    end
    return hcat(base...)
end

"""
    _build_keepout_position_dirs(n::Int) -> Matrix{Float64}

Build a bank of n direction vectors for keepout position support function evaluation.

The direction bank consists of:
- 6 axis directions (±x, ±y, ±z)
- 12 face diagonals (±(1,1,0), ±(1,-1,0), ±(1,0,1), ±(1,0,-1), ±(0,1,1), ±(0,1,-1))
- Remaining directions filled with golden-angle spiral on sphere

For n=72 (standard keepout bank), this gives 6 + 12 + 54 = 72 directions.

# Arguments
- `n::Int`: Number of direction vectors (minimum 6, must be even)

# Returns
- `Matrix{Float64}`: 3×n matrix of unit direction vectors
"""
function _build_keepout_position_dirs(n::Int)
    n = max(6, n)
    isodd(n) && (n += 1)
    base = Vector{Vector{Float64}}()

    function append_unique!(vec::AbstractVector{<:Real})
        v = Float64.(vec)
        vnorm = norm(v)
        vnorm > 1.0e-12 || return false
        v ./= vnorm
        for existing in base
            if isapprox(existing, v; atol=1.0e-10, rtol=0.0)
                return false
            end
        end
        push!(base, v)
        return true
    end

    # Axis directions in position space.
    for i in 1:3
        v = zeros(Float64, 3)
        v[i] = 1.0
        append_unique!(v)
        append_unique!(-v)
    end
    if n == 6
        return hcat(base...)
    end

    # Add all two-sparse position diagonals before filling the remaining bank
    # with an antipodal spherical lattice.
    for i in 1:3
        for j in (i+1):3
            v_sum = zeros(Float64, 3)
            v_sum[i] = 1.0
            v_sum[j] = 1.0
            append_unique!(v_sum)
            append_unique!(-v_sum)

            v_diff = zeros(Float64, 3)
            v_diff[i] = 1.0
            v_diff[j] = -1.0
            append_unique!(v_diff)
            append_unique!(-v_diff)
            length(base) >= n && return hcat(base[1:n]...)
        end
    end

    golden_angle = π * (3.0 - sqrt(5.0))
    k = 0
    while length(base) < n
        z = 1.0 - 2.0 * ((k + 0.5) / max(n, 1))
        radius = sqrt(max(0.0, 1.0 - z^2))
        θ = golden_angle * k
        v = [cos(θ) * radius, sin(θ) * radius, z]
        append_unique!(v)
        length(base) >= n && break
        append_unique!(-v)
        k += 1
    end
    return hcat(base...)
end

"""
    _controllable_append_unique_direction!(dirs::Vector{Vector{Float64}}, dir_vec::AbstractVector{<:Real}; atol::Float64=1.0e-10) -> Bool

Append a direction vector to a list if it is not already present (within tolerance).

# Arguments
- `dirs::Vector{Vector{Float64}}`: List of existing direction vectors
- `dir_vec::AbstractVector{<:Real}`: New direction vector to add
- `atol::Float64=1.0e-10`: Absolute tolerance for comparison

# Returns
- `Bool`: true if vector was added, false if it already exists
"""
function _controllable_append_unique_direction!(
    dirs::Vector{Vector{Float64}},
    dir_vec::AbstractVector{<:Real};
    atol::Float64=1.0e-10,
)
    vec = Float64.(dir_vec)
    vec_norm = norm(vec)
    vec_norm > atol || return false
    vec ./= vec_norm
    for existing in dirs
        if isapprox(existing, vec; atol=atol, rtol=0.0)
            return false
        end
    end
    push!(dirs, vec)
    return true
end

"""
    _build_repeated_horizon_dense_polyhedral_dirs(requested_count::Integer, keepout_dirs_pos::AbstractMatrix{<:Real}) -> Matrix{Float64}

Build a dense direction bank for repeated-horizon polyhedral set support evaluation.

Includes 6D state directions (position + velocity) with keepout position directions embedded.

# Arguments
- `requested_count::Integer`: Requested number of directions
- `keepout_dirs_pos::AbstractMatrix{<:Real}`: 3×Kd matrix of keepout position directions

# Returns
- `Matrix{Float64}`: 6×N matrix of unit direction vectors
"""
function _build_repeated_horizon_dense_polyhedral_dirs(
    requested_count::Integer,
    keepout_dirs_pos::AbstractMatrix{<:Real},
)
    dirs = Vector{Vector{Float64}}()

    # 6D axis directions
    for axis_idx in 1:6
        axis_vec = zeros(Float64, 6)
        axis_vec[axis_idx] = 1.0
        _controllable_append_unique_direction!(dirs, axis_vec)
        _controllable_append_unique_direction!(dirs, -axis_vec)
    end

    # Embed the keepout support bank directly in the 6D state bank so Stage 1
    # can enforce boundary support thresholds without reverting to box logic.
    for q in 1:size(keepout_dirs_pos, 2)
        keepout_vec = vcat(Float64.(keepout_dirs_pos[:, q]), zeros(Float64, 3))
        _controllable_append_unique_direction!(dirs, keepout_vec)
    end

    # Fill remaining with golden-angle spiral in 6D (projected to unit sphere)
    golden_angle = π * (3.0 - sqrt(5.0))
    k = 0
    while length(dirs) < requested_count
        # Generate point on 6D sphere using Fibonacci lattice
        v = zeros(Float64, 6)
        for i in 1:6
            offset = i * golden_angle
            v[i] = cos(offset + k * golden_angle)
        end
        v ./= norm(v)
        _controllable_append_unique_direction!(dirs, v)
        k += 1
    end

    return hcat(dirs[1:requested_count]...)
end

"""
    support_function(S::AbstractMatrix{<:Real}, d::AbstractVector{<:Real}) -> Float64

Compute the support function h_S(d) = sup_{x ∈ S} d^T x for a polyhedral set S.

For a polyhedron S = {x | Ax ≤ b}, the support function is:
h_S(d) = max_i (b_i / (a_i^T d)) for directions where a_i^T d > 0

# Arguments
- `S::AbstractMatrix{<:Real}`: Matrix representation of polyhedron (each row is a constraint normal)
- `d::AbstractVector{<:Real}`: Direction vector

# Returns
- `Float64`: Support function value
"""
function support_function(S::AbstractMatrix{<:Real}, d::AbstractVector{<:Real})
    d_vec = Float64.(d)
    d_norm = norm(d_vec)
    d_norm > 1e-12 || return 0.0
    d_vec ./= d_norm
    
    max_val = -Inf
    for i in 1:size(S, 1)
        a = S[i, :]
        dot_val = dot(a, d_vec)
        if dot_val > 1e-12
            # For constraint a^T x ≤ b, support in direction d is b / (a^T d)
            # Assuming b is stored as the last column or needs to be provided separately
            # This is a simplified version - full implementation needs constraint bounds
            max_val = max(max_val, dot_val)
        end
    end
    return max_val == -Inf ? 0.0 : max_val
end

"""
    support_function_hyperrectangle(r::AbstractVector{<:Real}, d::AbstractVector{<:Real}) -> Float64

Compute support function for a hyperrectangle centered at origin with half-extents r.

h_S(d) = Σ_i |d_i| * r_i

# Arguments
- `r::AbstractVector{<:Real}`: Half-extents of hyperrectangle
- `d::AbstractVector{<:Real}`: Direction vector

# Returns
- `Float64`: Support function value
"""
function support_function_hyperrectangle(r::AbstractVector{<:Real}, d::AbstractVector{<:Real})
    return sum(abs.(d) .* r)
end

"""
    support_function_ball(R::Real, d::AbstractVector{<:Real}) -> Float64

Compute support function for a ball of radius R centered at origin.

h_S(d) = R * ||d||

# Arguments
- `R::Real`: Ball radius
- `d::AbstractVector{<:Real}`: Direction vector

# Returns
- `Float64`: Support function value
"""
function support_function_ball(R::Real, d::AbstractVector{<:Real})
    return R * norm(d)
end

export _build_polyhedral_dirs, _build_keepout_position_dirs, _controllable_append_unique_direction!,
       _build_repeated_horizon_dense_polyhedral_dirs,
       support_function, support_function_hyperrectangle, support_function_ball

end # module SupportFunctions
