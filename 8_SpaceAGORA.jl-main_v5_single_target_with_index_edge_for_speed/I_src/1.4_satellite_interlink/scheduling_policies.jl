module SchedulingPolicies

using LinearAlgebra
using DiffEqBase
using ..InterLinkModels
using ..InterLinkModels: InterLinkMatchingWorkspace

export SchedulingPolicyModel, semimajor_axis_rate, score_candidates!
export select_interlinks!, schedule_interlinks!, interlink_scheduler_callback

struct SchedulingPolicyModel
	policy::Symbol
	target_idx::Int

	function SchedulingPolicyModel(policy::Symbol; target_idx::Integer=1)
		policy === :gve_sma || throw(ArgumentError("Supported scheduling policy: :gve_sma."))
		target_idx > 0 || throw(ArgumentError("target_idx must be a positive spacecraft vector index."))
		new(policy, Int(target_idx))
	end
end

SchedulingPolicyModel(; policy::Symbol=:gve_sma, target_idx::Integer=1) = SchedulingPolicyModel(policy; target_idx)

function semimajor_axis_rate(state, force, mu::Real)
	semimajor_axis = inv(2 / norm(state.pos) - dot(state.vel, state.vel) / mu)
	return 2 * semimajor_axis^2 / mu * dot(state.vel, force) / state.mass
end

function score_candidates!(model::InterLinkModel, policy::SchedulingPolicyModel, state, mu::Real)
	model.link_type === :laser || throw(ArgumentError("Radio force physics is not implemented."))
	isfinite(mu) && mu > 0 || throw(ArgumentError("The gravitational parameter must be finite and positive."))
	policy.target_idx <= length(state.sc) || throw(ArgumentError("target_idx is outside the spacecraft state vector."))
	for index in eachindex(model.connections)
		key = model.edge_keys[index]
		connection = model.connections[index]
		score = 0.0
		first_satellite, second_satellite = key[1][1], key[2][1]
		if connection.state.available && policy.target_idx in (first_satellite, second_satellite)
			target = state.sc[policy.target_idx]
			partner = state.sc[policy.target_idx == first_satellite ? second_satellite : first_satellite]
			score = semimajor_axis_rate(target, force_on_endpoint(target, partner, connection.parameters), mu)
		end
		isfinite(score) || throw(ArgumentError("Nonfinite $(policy.policy) score for candidate $key."))
		connection.state.score = score
	end
	return nothing
end

function matching_workspace!(model::InterLinkModel)
	model.matching !== nothing && return model.matching
	terminals = sort!(unique([endpoint for key in model.edge_keys for endpoint in key]))
	indices = Dict(endpoint => index for (index, endpoint) in enumerate(terminals))
	n, m = length(terminals), length(model.edge_keys)
	endpoints = Matrix{Int}(undef, 2, m)
	neighbors = [Int[] for _ in 1:n]
	for edge in sortperm(model.edge_keys)
		for side in 1:2
			vertex = indices[model.edge_keys[edge][side]]
			endpoints[side, edge] = vertex
			push!(neighbors[vertex], edge)
		end
	end
	model.matching = InterLinkMatchingWorkspace(endpoints, neighbors, zeros(Int, n),
		sizehint!(Int[], n), sizehint!(Int[], n), zeros(n),
		zeros(n), zeros(Int, n), fill(false, n), zeros(m), zeros(Int, n))
	return model.matching
end

@inline other_endpoint(work::InterLinkMatchingWorkspace, edge, vertex) =
	work.endpoints[1, edge] == vertex ? work.endpoints[2, edge] : work.endpoints[1, edge]

function select_tree!(model, work, first_index, last_index, selected)
	for index in last_index:-1:first_index
		vertex = work.order[index]
		base = 0.0
		for edge in work.neighbors[vertex]
			work.weights[edge] > 0 || continue
			child = other_endpoint(work, edge, vertex)
			work.parent[child] == vertex || continue
			base += work.free_score[child]
		end
		work.blocked_score[vertex] = base
		gain, chosen = 0.0, 0
		for edge in work.neighbors[vertex]
			work.weights[edge] > 0 || continue
			child = other_endpoint(work, edge, vertex)
			work.parent[child] == vertex || continue
			candidate = work.weights[edge] + work.blocked_score[child] - work.free_score[child]
			if candidate > gain
				gain, chosen = candidate, edge
			end
		end
		work.free_score[vertex] = base + gain
		work.choice[vertex] = chosen
	end
	for index in first_index:last_index
		vertex = work.order[index]
		edge = work.choice[vertex]
		(work.blocked[vertex] || edge == 0) && continue
		work.blocked[other_endpoint(work, edge, vertex)] = true
		model.connections[edge].state.active = true
		push!(selected, model.edge_keys[edge])
	end
	return nothing
end

function select_cyclic_component!(model, work, first_index, last_index, selected)
	vertices = sort!(work.order[first_index:last_index])
	for (index, vertex) in enumerate(vertices)
		work.local_index[vertex] = index
	end
	bits = [big(1) << (index - 1) for index in eachindex(vertices)]
	memo = Dict{BigInt, Tuple{Float64, Int}}()
	function best_matching(mask::BigInt)
		iszero(mask) && return (0.0, 0)
		haskey(memo, mask) && return memo[mask]
		index = Int(trailing_zeros(mask)) + 1
		vertex = vertices[index]
		remaining = mask & ~bits[index]
		best_score, _ = best_matching(remaining)
		chosen = 0
		for edge in work.neighbors[vertex]
			work.weights[edge] > 0 || continue
			partner = work.local_index[other_endpoint(work, edge, vertex)]
			iszero(remaining & bits[partner]) && continue
			score, _ = best_matching(remaining & ~bits[partner])
			score += work.weights[edge]
			if score > best_score
				best_score, chosen = score, edge
			end
		end
		return memo[mask] = (best_score, chosen)
	end
	mask = (big(1) << length(vertices)) - 1
	best_matching(mask)
	while !iszero(mask)
		index = Int(trailing_zeros(mask)) + 1
		_, edge = memo[mask]
		mask &= ~bits[index]
		edge == 0 && continue
		partner = work.local_index[other_endpoint(work, edge, vertices[index])]
		mask &= ~bits[partner]
		model.connections[edge].state.active = true
		push!(selected, model.edge_keys[edge])
	end
	return nothing
end

function select_interlinks!(model::InterLinkModel)
	work = matching_workspace!(model)
	fill!(work.parent, 0)
	fill!(work.blocked, false)
	empty!(work.order)
	empty!(work.stack)
	for (edge, connection) in enumerate(model.connections)
		connection.state.active = false
		work.weights[edge] = connection.state.available &&
			connection.state.score > model.active_link_penalty ?
			connection.state.score - model.active_link_penalty : 0.0
	end
	selected = InterLinkKey[]
	for root in eachindex(work.neighbors)
		work.parent[root] == 0 || continue
		work.parent[root] = -1
		first_index = length(work.order) + 1
		degree_sum = 0
		push!(work.stack, root)
		while !isempty(work.stack)
			vertex = pop!(work.stack)
			push!(work.order, vertex)
			for edge in work.neighbors[vertex]
				work.weights[edge] > 0 || continue
				degree_sum += 1
				partner = other_endpoint(work, edge, vertex)
				work.parent[partner] == 0 || continue
				work.parent[partner] = vertex
				push!(work.stack, partner)
			end
		end
		last_index = length(work.order)
		degree_sum == 0 && continue
		if degree_sum ÷ 2 == last_index - first_index
			select_tree!(model, work, first_index, last_index, selected)
		else
			select_cyclic_component!(model, work, first_index, last_index, selected)
		end
	end
	return sort!(selected)
end

function schedule_interlinks!(model::InterLinkModel, policy::SchedulingPolicyModel,
	spacecraft, state, mu::Real; is_active=nothing)
	update_availability!(model, spacecraft, state; is_active=is_active)
	score_candidates!(model, policy, state, mu)
	return select_interlinks!(model)
end

function interlink_scheduler_callback()
	function affect!(integrator)
		args = integrator.p.args
		model = args.interlink_model
		selected = schedule_interlinks!(model, args.scheduling_policy_model,
			args.dynamics_model.spacecraft, integrator.u, args.environment_model.planet.μ;
			is_active=integrator.p.is_active)
		sample = InterLinkScheduleSample(Float64(integrator.t), selected)
		if !isempty(model.history) && model.history[end].time == sample.time
			model.history[end] = sample
		else
			push!(model.history, sample)
		end
		DiffEqBase.u_modified!(integrator, true)
		return nothing
	end
	initialize = (callback, state, time, integrator) -> affect!(integrator)
	return DiffEqBase.DiscreteCallback((state, time, integrator) -> true, affect!;
		initialize=initialize, save_positions=(false, false))
end

end
