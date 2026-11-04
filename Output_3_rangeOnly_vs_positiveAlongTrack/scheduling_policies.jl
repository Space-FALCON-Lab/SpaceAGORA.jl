module ComparisonPolicies

using LinearAlgebra, StaticArrays

export LaserPolicy, update_schedule!, pair_force, along_track_projection

const POLICY_NAMES = (:naive_next_entering, :positive_along_track)

mutable struct LaserPolicy
    schedule::Symbol
    target_idx::Int
    helper_indices::Vector{Int}
    range_m::Float64
    force_N::Float64
    active_helper_idx::Int
    previous_in_range::Vector{Bool}
    in_range::Vector{Bool}
    distances::Vector{Float64}
    activations::Int
end

function LaserPolicy(schedule::Symbol; target_idx::Int=1,
    helper_indices::AbstractVector{Int}, range_m::Real=200_000.0,
    force_N::Real=100 * 10_000 / 299_792_458)
    schedule in POLICY_NAMES || throw(ArgumentError("Unsupported policy: $schedule"))
    target_idx > 0 || throw(ArgumentError("Target index must be positive"))
    !isempty(helper_indices) && all(>(0), helper_indices) &&
        all(!=(target_idx), helper_indices) &&
        length(unique(helper_indices)) == length(helper_indices) ||
        throw(ArgumentError("Helpers must have unique, positive, non-target indices"))
    isfinite(range_m) && range_m >= 0 ||
        throw(ArgumentError("Range must be finite and nonnegative"))
    isfinite(force_N) && force_N >= 0 ||
        throw(ArgumentError("Force must be finite and nonnegative"))
    n = length(helper_indices)
    LaserPolicy(schedule, target_idx, collect(helper_indices), Float64(range_m),
        Float64(force_N), 0, fill(false, n), fill(false, n), zeros(n), 0)
end

function along_track_projection(target_pos, target_vel, helper_pos)
    displacement = target_pos - helper_pos
    distance = norm(displacement)
    distance == 0 && return 0.0
    radius = norm(target_pos)
    angular_momentum = cross(target_pos, target_vel)
    radius > 0 && norm(angular_momentum) > 0 ||
        throw(ArgumentError("Along-track direction requires a nondegenerate target orbit"))
    transverse = cross(normalize(angular_momentum), target_pos / radius)
    dot(displacement / distance, transverse)
end

function closest_helper(policy::LaserPolicy; entering_only=false)
    best_idx, best_distance = 0, Inf
    for slot in eachindex(policy.helper_indices)
        policy.in_range[slot] || continue
        entering_only && policy.previous_in_range[slot] && continue
        if policy.distances[slot] < best_distance
            best_idx = policy.helper_indices[slot]
            best_distance = policy.distances[slot]
        end
    end
    best_idx
end

"""
Reproduce the two single-target policies in 2_SpaceAGORA.jl.
Call at initialization and after each accepted solver step, never inside the RHS.
Range-only handovers consider newly entering helpers, not every available helper.
Both policies retain their current helper while its policy-specific conditions hold.
"""
function update_schedule!(policy::LaserPolicy, positions, velocities)
    length(positions) == length(velocities) ||
        throw(ArgumentError("Position and velocity counts differ"))
    max(policy.target_idx, maximum(policy.helper_indices)) <= length(positions) ||
        throw(ArgumentError("Spacecraft index is outside the state"))
    target_pos, target_vel = positions[policy.target_idx], velocities[policy.target_idx]
    old_helper = policy.active_helper_idx
    for slot in eachindex(policy.helper_indices)
        distance = norm(positions[policy.helper_indices[slot]] - target_pos)
        isfinite(distance) || throw(ArgumentError("Nonfinite helper distance"))
        policy.distances[slot] = distance
        policy.in_range[slot] = distance <= policy.range_m
    end

    if policy.schedule === :naive_next_entering
        if policy.active_helper_idx > 0
            slot = findfirst(==(policy.active_helper_idx), policy.helper_indices)
            if slot === nothing || !policy.in_range[slot]
                policy.active_helper_idx = closest_helper(policy; entering_only=true)
            end
        else
            policy.active_helper_idx = closest_helper(policy;
                entering_only=any(policy.previous_in_range))
        end
    elseif policy.schedule === :positive_along_track
        if policy.active_helper_idx > 0
            slot = findfirst(==(policy.active_helper_idx), policy.helper_indices)
            valid = slot !== nothing && policy.in_range[slot] &&
                along_track_projection(target_pos, target_vel,
                    positions[policy.active_helper_idx]) > 0
            valid || (policy.active_helper_idx = 0)
        end
        if policy.active_helper_idx == 0
            best_idx, best_projection = 0, 0.0
            for slot in eachindex(policy.helper_indices)
                policy.in_range[slot] || continue
                helper = policy.helper_indices[slot]
                projection = along_track_projection(target_pos, target_vel, positions[helper])
                if projection > best_projection
                    best_idx, best_projection = helper, projection
                end
            end
            policy.active_helper_idx = best_idx
        end
    else
        throw(ArgumentError("Unsupported policy: $(policy.schedule)"))
    end

    if policy.active_helper_idx > 0 && policy.active_helper_idx != old_helper
        policy.activations += 1
    end
    policy.previous_in_range .= policy.in_range
    policy.active_helper_idx
end

function pair_force(policy::LaserPolicy, target_pos, helper_pos)
    displacement = target_pos - helper_pos
    distance = norm(displacement)
    isfinite(distance) || throw(ArgumentError("Nonfinite laser separation"))
    0 < distance <= policy.range_m || return zero(SVector{3,Float64})
    policy.force_N * displacement / distance
end

end
