# Shared joint-path geometry and clearance. Historical qualified names are retained.

"""Sample a robot-arm HYPR path as a Bezier curve or polyline."""
function robot_arm_sample_hypr_path(points, n_samples::Int; curve_type::Symbol=:bezier)
    return hypr_sample_count_path(points, n_samples; curve_type=curve_type)
end


"""Return joint-space length of sampled robot-arm path points."""
function _robot_arm_path_length(samples)
    return hypr_path_length(samples)
end


"""Compute squared second-difference smoothness for a robot-arm path."""
function _robot_arm_path_smoothness(samples)
    size(samples, 2) < 3 && return 0.0
    total = 0.0
    @inbounds for k in 2:(size(samples, 2) - 1)
        total += sum(abs2, samples[:, k + 1] .- 2.0 .* samples[:, k] .+ samples[:, k - 1])
    end
    return total / max(size(samples, 2) - 2, 1)
end


"""Return the distance from a point to a 3D segment."""
function _robot_arm_segment_distance(p::SVector{3, Float64}, a::SVector{3, Float64}, b::SVector{3, Float64})
    ab = b - a
    denom = dot(ab, ab)
    denom <= eps(Float64) && return norm(p - a)
    t = clamp(dot(p - a, ab) / denom, 0.0, 1.0)
    return norm(p - (a + t * ab))
end

"""Compute obstacle-clearance diagnostics for sampled robot-arm joint paths."""
function robot_arm_clearance_stats_from_samples(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    q_samples,
    obstacles::AbstractVector{RobotArmSphereObstacle},
    safe_distance_m::Real,
)
    isempty(obstacles) && return (
        min_clearance=Inf,
        violation_count=0,
        violation_fraction=0.0,
        clearance_penalty=0.0,
    )
    safe = Float64(safe_distance_m)
    min_clearance = Inf
    violation_count = 0
    clearance_penalty = 0.0
    checks = 0
    @inbounds for k in axes(q_samples, 2)
        pose = cloth_fk(model, base_pose, q_samples[:, k])
        for i in eachindex(model.links)
            a = pose.joint_origins[i]
            b = pose.link_tip_positions[i]
            link_radius = model.links[i].radius_m
            for obs in obstacles
                clearance = _robot_arm_segment_distance(obs.center, a, b) - obs.radius_m - link_radius
                min_clearance = min(min_clearance, clearance)
                checks += 1
                deficit = safe - clearance
                if deficit > 0.0
                    violation_count += 1
                    clearance_penalty += deficit * deficit
                end
            end
        end
    end
    return (
        min_clearance=min_clearance,
        violation_count=violation_count,
        violation_fraction=violation_count / max(checks, 1),
        clearance_penalty=clearance_penalty,
    )
end
