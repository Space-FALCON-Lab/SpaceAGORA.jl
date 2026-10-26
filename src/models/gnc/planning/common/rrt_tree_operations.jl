# Tree operations shared by RPO and robot-arm RRT planning.
# Loaded into HYPRUtils to preserve existing qualified access.
"""Return the parent-index vector from either RPO-style or robot-arm-style RRT tree storage."""
_hypr_rrt_parents(tree) = hasproperty(tree, :parents) ? getproperty(tree, :parents) : getproperty(tree, :parent)
"""Return the accumulated-cost vector from either RPO-style or robot-arm-style RRT tree storage."""
_hypr_rrt_costs(tree) = hasproperty(tree, :costs) ? getproperty(tree, :costs) : getproperty(tree, :cost)

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
