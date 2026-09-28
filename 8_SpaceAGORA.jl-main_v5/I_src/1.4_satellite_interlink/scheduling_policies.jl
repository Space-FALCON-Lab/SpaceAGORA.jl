module SchedulingPolicies

using LinearAlgebra
using DiffEqBase
using ..InterLinkModels

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
	for (key, connection) in model.linkgraph
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

function select_interlinks!(model::InterLinkModel)
	candidates = sort!([key for (key, connection) in model.linkgraph
		if connection.state.available && connection.state.score > model.active_link_penalty])
	endpoints = sort!(unique([endpoint for key in candidates for endpoint in key]))
	bits = Dict(endpoint => big(1) << (index - 1) for (index, endpoint) in enumerate(endpoints))
	neighbors = Dict(endpoint => InterLinkKey[] for endpoint in endpoints)
	for key in candidates, endpoint in key
		push!(neighbors[endpoint], key)
	end
	memo = Dict{BigInt, Tuple{Float64, Vector{InterLinkKey}}}()

	function best_matching(mask::BigInt)
		iszero(mask) && return (0.0, InterLinkKey[])
		haskey(memo, mask) && return memo[mask]
		endpoint = endpoints[findfirst(endpoint -> !iszero(mask & bits[endpoint]), endpoints)]
		remaining = mask & ~bits[endpoint]
		best_score, best_keys = best_matching(remaining)
		for key in neighbors[endpoint]
			partner = key[1] == endpoint ? key[2] : key[1]
			iszero(remaining & bits[partner]) && continue
			score, selected = best_matching(remaining & ~bits[partner])
			score += model.linkgraph[key].state.score - model.active_link_penalty
			if score > best_score
				best_score, best_keys = score, vcat(InterLinkKey[key], selected)
			end
		end
		memo[mask] = (best_score, best_keys)
		return memo[mask]
	end

	_, selected = best_matching((big(1) << length(endpoints)) - 1)
	selected = sort!(copy(selected))
	selected_set = Set(selected)
	for (key, connection) in model.linkgraph
		connection.state.active = key in selected_set
	end
	return selected
end

function schedule_interlinks!(model::InterLinkModel, policy::SchedulingPolicyModel,
	spacecraft, state, mu::Real; is_active=nothing)
	for key in keys(model.linkgraph)
		update_availability!(model, key, spacecraft, state; is_active=is_active)
	end
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
