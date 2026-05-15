function rpo_goal_standoff_point(surface_point_body, surface_normal_body, standoff_m::Real)
    n = SVector{3, Float64}(surface_normal_body)
    n_norm = norm(n)
    n_norm <= eps(Float64) && throw(ArgumentError("surface normal must be nonzero."))
    return SVector{3, Float64}(surface_point_body) + Float64(standoff_m) * (n / n_norm)
end

