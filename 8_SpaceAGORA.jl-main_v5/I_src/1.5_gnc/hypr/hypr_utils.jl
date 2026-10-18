"""Shared HYPR helper module for path sampling, PSO bookkeeping, and RRT tree operations reused by RPO and robot-arm planners."""
module HYPRUtils

using LinearAlgebra

export hypr_path_length, hypr_bezier_point, hypr_bezier_point!, hypr_sample_count_path
export hypr_iteration_weights, hypr_material_improvement, hypr_protected_particle_mask
export hypr_rrt_nearest_index, hypr_rrt_near_indices, hypr_rrt_steer
export hypr_rrt_tree_path, hypr_rrt_join_paths, hypr_rrt_refresh_subtree_costs!

"""Return the parent-index vector from either RPO-style or robot-arm-style RRT tree storage."""
_hypr_rrt_parents(tree) = hasproperty(tree, :parents) ? getproperty(tree, :parents) : getproperty(tree, :parent)
"""Return the accumulated-cost vector from either RPO-style or robot-arm-style RRT tree storage."""
_hypr_rrt_costs(tree) = hasproperty(tree, :costs) ? getproperty(tree, :costs) : getproperty(tree, :cost)

"""Return the Euclidean arc length of a waypoint matrix."""
function hypr_path_length(points)
    pts = Matrix{Float64}(points)
    size(pts, 2) < 2 && return 0.0
    total = 0.0
    @inbounds for j in 1:(size(pts, 2) - 1)
        total += norm(pts[:, j + 1] - pts[:, j])
    end
    return total
end

"""Evaluate a Bezier curve at a normalized parameter without mutating caller-owned work buffers."""
function hypr_bezier_point(points, t::Real)
    pts = Matrix{Float64}(points)
    out = zeros(size(pts, 1))
    work = similar(pts)
    return hypr_bezier_point!(out, work, pts, Float64(t))
end

"""Evaluate a Bezier curve in-place using caller-provided output and work buffers."""
function hypr_bezier_point!(out, work, points, t::Float64)
    n = size(points, 2)
    work[:, 1:n] .= points
    @inbounds for r in 1:(n - 1)
        for j in 1:(n - r)
            work[:, j] .= (1 - t) .* work[:, j] .+ t .* work[:, j + 1]
        end
    end
    out .= work[:, 1]
    return out
end

"""Sample either a Bezier or polyline path at a fixed number of points."""
function hypr_sample_count_path(points, n_samples::Int; curve_type::Symbol=:bezier)
    n_samples >= 2 || throw(ArgumentError("n_samples must be at least 2."))
    pts = Matrix{Float64}(points)
    n_dims = size(pts, 1)
    n_points = size(pts, 2)
    samples = zeros(n_dims, n_samples)
    if curve_type == :bezier
        degree = n_points - 1
        @inbounds for k in 1:n_samples
            s = (k - 1) / max(n_samples - 1, 1)
            for i in 0:degree
                coeff = binomial(degree, i) * (1.0 - s)^(degree - i) * s^i
                samples[:, k] .+= coeff .* pts[:, i + 1]
            end
        end
    elseif curve_type == :polyline
        @inbounds for k in 1:n_samples
            u = (n_points - 1) * (k - 1) / max(n_samples - 1, 1)
            seg = clamp(floor(Int, u) + 1, 1, n_points - 1)
            alpha = clamp(u - (seg - 1), 0.0, 1.0)
            samples[:, k] .= (1.0 - alpha) .* pts[:, seg] .+ alpha .* pts[:, seg + 1]
        end
    else
        throw(ArgumentError("curve_type must be :bezier or :polyline."))
    end
    samples[:, 1] .= pts[:, 1]
    samples[:, end] .= pts[:, end]
    return samples
end

"""Interpolate PSO inertia and acceleration weights for the requested iteration."""
function hypr_iteration_weights(
    schedule_enable::Bool,
    n_iters::Int,
    iter::Int,
    w_inertia::Real,
    c1::Real,
    c2::Real,
    transition_fraction::Real,
    w_min::Real,
    w_end_fraction::Real,
    c1_end_fraction::Real,
    c2_end_fraction::Real,
    c_min::Real,
    c_max::Real,
)
    if !schedule_enable || n_iters <= 1
        return (w_inertia=Float64(w_inertia), c1=Float64(c1), c2=Float64(c2))
    end
    transition = clamp(Float64(transition_fraction), 1.0e-6, 1.0)
    progress = clamp((iter - 1) / ((n_iters - 1) * transition), 0.0, 1.0)
    smooth = progress^2 * (3.0 - 2.0 * progress)
    w0 = Float64(w_inertia)
    c10 = Float64(c1)
    c20 = Float64(c2)
    w_end = max(Float64(w_min), Float64(w_end_fraction) * w0)
    c1_end = clamp(Float64(c1_end_fraction) * c10, Float64(c_min), Float64(c_max))
    c2_end = clamp(Float64(c2_end_fraction) * c20, Float64(c_min), Float64(c_max))
    return (
        w_inertia=w0 + smooth * (w_end - w0),
        c1=c10 + smooth * (c1_end - c10),
        c2=c20 + smooth * (c2_end - c20),
    )
end

"""Decide whether a new cost improves enough in absolute or relative terms to reset stagnation logic."""
function hypr_material_improvement(new_cost::Real, reference_cost::Real, min_abs_improvement::Real, min_rel_improvement::Real)::Bool
    new = Float64(new_cost)
    reference = Float64(reference_cost)
    isfinite(new) || return false
    isfinite(reference) || return true
    improvement = reference - new
    threshold = max(
        Float64(min_abs_improvement),
        Float64(min_rel_improvement) * max(abs(reference), 1.0),
    )
    return improvement > threshold
end

"""Mark elite finite-cost particles that should be protected from swarm culling."""
function hypr_protected_particle_mask(costs, elite_fraction)
    n = length(costs)
    mask = falses(n)
    n == 0 && return mask
    elite_count = clamp(Int(ceil(clamp(Float64(elite_fraction), 0.0, 1.0) * n)), 1, n)
    ranked = sortperm(costs)
    for idx in ranked[1:elite_count]
        mask[idx] = true
    end
    return mask
end

"""Return the index of the RRT node nearest to the query state."""
function hypr_rrt_nearest_index(tree, q)
    best_idx = 1
    best_d2 = Inf
    @inbounds for i in eachindex(tree.nodes)
        d2 = sum(abs2, tree.nodes[i] - q)
        if d2 < best_d2
            best_d2 = d2
            best_idx = i
        end
    end
    return best_idx
end

"""Return indices of RRT nodes within a search radius of the query state."""
function hypr_rrt_near_indices(tree, q, radius)
    r2 = Float64(radius)^2
    idxs = Int[]
    @inbounds for i in eachindex(tree.nodes)
        sum(abs2, tree.nodes[i] - q) <= r2 && push!(idxs, i)
    end
    return idxs
end

"""Move from one RRT state toward another by at most the configured step size."""
function hypr_rrt_steer(q_near, q_target, step_size; trap_tol::Real=1.0e-10, step_floor::Real=1.0e-9)
    direction = q_target - q_near
    distance = norm(direction)
    distance <= Float64(trap_tol) && return q_near, :trapped
    step = max(Float64(step_size), Float64(step_floor))
    distance <= step && return q_target, :reached
    return q_near + (step / distance) * direction, :advanced
end

"""Reconstruct a root-to-node path from an RRT tree parent chain."""
function hypr_rrt_tree_path(tree, idx::Integer)
    nodes = typeof(tree.nodes[1])[]
    current = Int(idx)
    parents = _hypr_rrt_parents(tree)
    while current > 0
        push!(nodes, tree.nodes[current])
        current = parents[current]
    end
    reverse!(nodes)
    return reduce(hcat, nodes)
end

"""Join start-side and goal-side RRT paths at a connection point."""
function hypr_rrt_join_paths(start_tree, start_idx::Integer, goal_tree, goal_idx::Integer)
    start_path = hypr_rrt_tree_path(start_tree, start_idx)
    goal_path = hypr_rrt_tree_path(goal_tree, goal_idx)
    return hcat(start_path, reverse(goal_path[:, 1:(end - 1)]; dims=2))
end

"""Recompute accumulated costs for descendants after an RRT parent rewiring."""
function hypr_rrt_refresh_subtree_costs!(tree, parent_idx::Integer)
    queue = [Int(parent_idx)]
    parents = _hypr_rrt_parents(tree)
    costs = _hypr_rrt_costs(tree)
    while !isempty(queue)
        parent = popfirst!(queue)
        @inbounds for idx in eachindex(parents)
            if parents[idx] == parent
                costs[idx] = costs[parent] + norm(tree.nodes[idx] - tree.nodes[parent])
                push!(queue, idx)
            end
        end
    end
    return tree
end

end # module HYPRUtils
