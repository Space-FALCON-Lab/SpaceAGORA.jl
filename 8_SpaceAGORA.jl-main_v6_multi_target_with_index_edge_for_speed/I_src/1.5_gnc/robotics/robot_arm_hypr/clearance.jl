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
