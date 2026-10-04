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
