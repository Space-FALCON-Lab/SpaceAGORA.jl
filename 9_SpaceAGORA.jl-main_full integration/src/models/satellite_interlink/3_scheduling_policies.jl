# POLICY 1. Increase Target Satellite Semimajor Axis Rate
function score_gve_sma!(model::ConstellationInterLinkModel, policy::SchedulingPolicyModel, state, mu::Real;
	candidate_force)
	function semimajor_axis_rate(target, force, mu::Real)
		semimajor_axis = inv(2 / norm(target.pos) - dot(target.vel, target.vel) / mu)
		return 2 * semimajor_axis^2 / mu * dot(target.vel, force) / target.mass
	end
	return score_target_rates!(semimajor_axis_rate, model, policy, state, mu; candidate_force)
end

function score_target_rates!(rate, model::ConstellationInterLinkModel, policy::SchedulingPolicyModel,
	state, mu::Real; candidate_force)
	isfinite(mu) && mu > 0 || throw(ArgumentError("The gravitational parameter must be finite and positive."))
	policy.target_idxs === nothing || all(index -> index <= length(state.sc), policy.target_idxs) ||
		throw(ArgumentError("target_idxs contains an index outside the spacecraft state vector."))
	policy.target_idxs !== nothing || policy.target_idx === nothing || policy.target_idx <= length(state.sc) ||
		throw(ArgumentError("target_idx is outside the spacecraft state vector."))
	work = matching_workspace!(model, policy)
	for index in eachindex(model.graph.edges)
		key = map(v -> model.graph.vertices[v], model.graph.edges[index])
		score = 0.0
		first_satellite, second_satellite = key[1][1], key[2][1]
		if model.graph.available[index]
			for (target_idx, partner_idx) in ((first_satellite, second_satellite), (second_satellite, first_satellite))
				if policy.target_idxs === nothing
					policy.target_idx === nothing || policy.target_idx == target_idx || continue
				else
					target_idx in policy.target_idxs || continue
				end
				target = state.sc[target_idx]
				partner = state.sc[partner_idx]
				force = candidate_force(target, partner, model.graph.edge_properties[index])
				score += rate(target, force, mu)
			end
		end
		isfinite(score) || throw(ArgumentError("Nonfinite $(policy.policy) score for candidate $key."))
		work.score[index] = score
	end
	return nothing
end

"""
Instantaneous inclination rate relative to the inertial xy plane, in radians/s.
Uses atan(norm(h_xy), h_z) to avoid acos cancellation near equatorial orbits.
At exactly equatorial states, returns the forward one-sided rate of inclination:
positive at i=0 and negative at i=pi for a nonzero plane-tilting acceleration.
"""
function inclination_rate(target, force, mu::Real)
	isfinite(target.mass) && target.mass > 0 ||
		throw(ArgumentError("Inclination scoring requires finite, positive target mass."))
	h = cross(target.pos, target.vel)
	h_norm = norm(h)
	isfinite(h_norm) && h_norm > 0 ||
		throw(ArgumentError("Inclination scoring requires finite, nonzero angular momentum."))
	dh = cross(target.pos, force / target.mass)
	all(isfinite, dh) || throw(ArgumentError("Nonfinite angular momentum rate in inclination scoring."))
	rho = hypot(h[1], h[2])
	if rho == 0
		return sign(h[3]) * hypot(dh[1], dh[2]) / h_norm
	end
	drho = (h[1] / rho) * dh[1] + (h[2] / rho) * dh[2]
	return ((h[3] / h_norm) * drho - (rho / h_norm) * dh[3]) / h_norm
end

function score_gve_inc!(model::ConstellationInterLinkModel, policy::SchedulingPolicyModel, state, mu::Real;
	candidate_force)
	return score_target_rates!(inclination_rate, model, policy, state, mu; candidate_force)
end

# POLICY 2. Increase Target Satellite Eccentricity Rate
function score_gve_eccentricity!(model::ConstellationInterLinkModel, policy::SchedulingPolicyModel, state, mu::Real;
	candidate_force)
	function eccentricity_rate(target, force, mu::Real)
		position = target.pos
		velocity = target.vel
		acceleration = force / target.mass
		angular_momentum = cross(position, velocity)
		eccentricity_vector = cross(velocity, angular_momentum) / mu - position / norm(position)
		eccentricity_vector_rate = (
			cross(acceleration, angular_momentum) + cross(velocity, cross(position, acceleration))
		) / mu
		eccentricity = norm(eccentricity_vector)
		return eccentricity > sqrt(eps(Float64)) ?
			dot(eccentricity_vector, eccentricity_vector_rate) / eccentricity : norm(eccentricity_vector_rate)
	end

	isfinite(mu) && mu > 0 || throw(ArgumentError("The gravitational parameter must be finite and positive."))
	policy.target_idx === nothing || policy.target_idx <= length(state.sc) ||
		throw(ArgumentError("target_idx is outside the spacecraft state vector."))
	work = matching_workspace!(model, policy)
	for index in eachindex(model.graph.edges)
		key = map(v -> model.graph.vertices[v], model.graph.edges[index])
		score = 0.0
		first_satellite, second_satellite = key[1][1], key[2][1]
		if model.graph.available[index]
			for (target_idx, partner_idx) in ((first_satellite, second_satellite), (second_satellite, first_satellite))
				policy.target_idx === nothing || policy.target_idx == target_idx || continue
				target = state.sc[target_idx]
				partner = state.sc[partner_idx]
				force = candidate_force(target, partner, model.graph.edge_properties[index])
				score += eccentricity_rate(target, force, mu)
			end
		end
		isfinite(score) || throw(ArgumentError("Nonfinite $(policy.policy) score for candidate $key."))
		work.score[index] = score
	end
	return nothing
end

# POLICY 3: Increase the SMA of the lower-index satellite in each interlink pair.
"""
Score each available link by the SMA rate of its lower spacecraft vector index,
regardless of altitude or registration order. The target is chosen per pair,
so `policy.target_idx` is ignored. Nonpositive scores are excluded by the matcher.
"""
function score_gve_sma_first_index!(model::ConstellationInterLinkModel, policy::SchedulingPolicyModel, state, mu::Real;
	candidate_force)
	isfinite(mu) && mu > 0 || throw(ArgumentError("The gravitational parameter must be finite and positive."))
	work = matching_workspace!(model, policy)
	for index in eachindex(model.graph.edges)
		key = map(v -> model.graph.vertices[v], model.graph.edges[index])
		score = 0.0
		if model.graph.available[index]
			target_idx, partner_idx = minmax(key[1][1], key[2][1])
			target = state.sc[target_idx]
			partner = state.sc[partner_idx]
			force = candidate_force(target, partner, model.graph.edge_properties[index])
			semimajor_axis = inv(2 / norm(target.pos) - dot(target.vel, target.vel) / mu)
			score = 2 * semimajor_axis^2 / mu * dot(target.vel, force) / target.mass
		end
		isfinite(score) || throw(ArgumentError("Nonfinite $(policy.policy) score for candidate $key."))
		work.score[index] = score
	end
	return nothing
end
