module SchedulingPolicies

using LinearAlgebra
using DiffEqBase
import PythonCall
using ..InterLinkModels

export SchedulingPolicyModel, score_candidates!
export select_interlinks!, schedule_interlinks!, interlink_scheduler_callback

# 1. Define the scheduling policy model structure
struct SchedulingPolicyModel
	policy::Symbol
	target_idx::Union{Nothing, Int}

	function SchedulingPolicyModel(policy::Symbol; target_idx::Union{Nothing, Integer}=1)
		policy in (:gve_sma, :gve_eccentricity) ||
			throw(ArgumentError("Supported scheduling policies: :gve_sma and :gve_eccentricity."))
		target_idx === nothing || target_idx > 0 ||
			throw(ArgumentError("target_idx must be nothing or a positive spacecraft vector index."))
		new(policy, target_idx === nothing ? nothing : Int(target_idx))
	end
end

include(joinpath(@__DIR__, "scheduling_policies.jl"))

# 2.1. Define the function to score candidate interlinks based on the scheduling policy
function score_candidates!(model::InterLinkModel, policy::SchedulingPolicyModel, state, mu::Real;
	# compute the score of each interlink based on the scheduling policy / scoring policy
	# resulting scores will be stored in interlink_connection.state.score
	candidate_force)
	if policy.policy === :gve_sma
		return score_gve_sma!(model, policy, state, mu; candidate_force)
	elseif policy.policy === :gve_eccentricity
		return score_gve_eccentricity!(model, policy, state, mu; candidate_force)
	end
	throw(ArgumentError("Unsupported scheduling policy: $(policy.policy)."))
end

# 2.2. Define the function to select the best set of interlinks based on the current scores
"""
	select_interlinks!(model::InterLinkModel)

Activate a maximum-total-score matching of available, positive-score links,
with at most one link per terminal. Uses NetworkX's weighted Blossom algorithm
through PythonCall, with O(N^3) arithmetic operations for N candidate endpoints.
Float64 scores are converted to exact rationals and scaled by their largest
denominator (all are powers of two), preserving the objective with integer weights.
Sorted graph insertion makes ties repeatable for a given backend version, but
does not preserve the previous recursive selector's choice among equal optima.
PythonCall provisions NetworkX from CondaPkg.toml on first use; setup needs network access.
"""
function select_interlinks!(model::InterLinkModel)
	# select the best set of interlinks based on the current scores, then set them to be active links
	# In the selected links will have their interlink_connection.state.active set to true
	candidates = sort!([key for (key, connection) in model.linkgraph
		if connection.state.available && connection.state.score > 0.0])
	endpoints = sort!(unique([endpoint for key in candidates for endpoint in key]))
	selected = InterLinkKey[]
	if !isempty(candidates)
		endpoint_indices = Dict(endpoint => index for (index, endpoint) in enumerate(endpoints))
		scores = [model.linkgraph[key].state.score for key in candidates]
		all(isfinite, scores) || throw(ArgumentError("Interlink matching scores must be finite."))
		exact_scores = Rational{BigInt}.(scores)
		common_denominator = maximum(denominator, exact_scores)
		networkx = PythonCall.pyimport("networkx")
		graph = networkx.Graph()
		for (key, score) in zip(candidates, exact_scores)
			weight = numerator(score) * div(common_denominator, denominator(score))
			graph.add_edge(endpoint_indices[key[1]], endpoint_indices[key[2]];
				weight=PythonCall.pybuiltins.int(string(weight)))
		end
		matching = PythonCall.pyconvert(Vector{Tuple{Int, Int}},
			PythonCall.pybuiltins.list(networkx.max_weight_matching(graph; maxcardinality=false)))
		for (first_index, second_index) in matching
			push!(selected, minmax(endpoints[first_index], endpoints[second_index]))
		end
		sort!(selected)
	end
	selected_set = Set(selected)
	for (key, connection) in model.linkgraph
		connection.state.active = key in selected_set
	end
	return selected
end

# 2.3. Define the function to schedule interlinks based on the current state and scoring policy
function schedule_interlinks!(model::InterLinkModel, policy::SchedulingPolicyModel,
	spacecraft, state, mu::Real; candidate_force, is_active=nothing)
	# step 1. update the availability of each interlink based on the current state
	for key in keys(model.linkgraph)
		update_availability!(model, key, spacecraft, state; is_active=is_active) 
	end

	# step 2. compute the score of each interlink based on the scheduling policy / scoring policy
	score_candidates!(model, policy, state, mu; candidate_force) 

	# step 3. select the best set of interlinks based on the current scores
	return select_interlinks!(model)
end

# 3. Define the callback function for the interlink graph update with scheduler
function interlink_scheduler_callback(; candidate_force)
	function affect!(integrator)
		args = integrator.p.args
		model = args.interlink_model
		selected = schedule_interlinks!(model, args.scheduling_policy_model,
			args.dynamics_model.spacecraft, integrator.u, args.environment_model.planet.μ;
			candidate_force, is_active=integrator.p.is_active)
		sample = InterLinkScheduleSample(Float64(integrator.t), selected)
		if !isempty(model.history) && model.history[end].time == sample.time
			model.history[end] = sample
		else
			push!(model.history, sample)
		end
		DiffEqBase.u_modified!(integrator, true)
		return nothing
	end
	initialize = (callback, state, time, integrator) -> affect!(integrator) # update link graph when the simulation is initialized
	return DiffEqBase.DiscreteCallback((state, time, integrator) -> true, affect!; # update link graph at every internal simulation step
		initialize=initialize, save_positions=(false, false))
end

end
