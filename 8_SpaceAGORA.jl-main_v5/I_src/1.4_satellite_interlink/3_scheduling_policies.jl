# POLICY 1. Increase Target Satellite Semimajor Axis Rate
function score_gve_sma!(model::ConstellationInterLinkModel, policy::SchedulingPolicyModel, state, mu::Real;
	candidate_force)
	function semimajor_axis_rate(target, force, mu::Real)
		semimajor_axis = inv(2 / norm(target.pos) - dot(target.vel, target.vel) / mu)
		return 2 * semimajor_axis^2 / mu * dot(target.vel, force) / target.mass
	end

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
				score += semimajor_axis_rate(target, force, mu)
			end
		end
		isfinite(score) || throw(ArgumentError("Nonfinite $(policy.policy) score for candidate $key."))
		work.score[index] = score
	end
	return nothing
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
