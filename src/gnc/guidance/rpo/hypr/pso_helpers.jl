function rpo_pso_bounds(start_rtn, goal_rtn, cfg::RPOPSOConfig)
    start = SVector{3, Float64}(start_rtn)
    goal = SVector{3, Float64}(goal_rtn)
    lo = min.(start, goal) .- cfg.search_margin_m
    hi = max.(start, goal) .+ cfg.search_margin_m
    span = max(norm(goal - start), cfg.search_margin_m)
    lo = min.(lo, 0.5 .* (start + goal) .- cfg.spread_scale * span)
    hi = max.(hi, 0.5 .* (start + goal) .+ cfg.spread_scale * span)
    return lo, hi
end

function rpo_position_to_path(position, start_rtn, goal_rtn, n_waypoints::Int)
    start = SVector{3, Float64}(start_rtn)
    goal = SVector{3, Float64}(goal_rtn)
    points = zeros(3, n_waypoints + 2)
    points[:, 1] .= start
    points[:, end] .= goal
    for j in 1:n_waypoints
        points[:, j + 1] .= position[(3 * (j - 1) + 1):(3 * j)]
    end
    return points
end

