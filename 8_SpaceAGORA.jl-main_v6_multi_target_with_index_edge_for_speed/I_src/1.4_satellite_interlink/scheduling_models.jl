module SchedulingPolicies

using LinearAlgebra
using DiffEqBase
import PythonCall
using ..InterLinkModels
using ..InterLinkModels: InterLinkMatchingWorkspace

export SchedulingPolicyModel, score_candidates!
export select_interlinks!, schedule_interlinks!, interlink_scheduler_callback

# 1. Define the scheduling policy model structure
struct SchedulingPolicyModel
	policy::Symbol
	target_idx::Union{Nothing, Int}

	function SchedulingPolicyModel(policy::Symbol; target_idx::Union{Nothing, Integer}=1)
		policy in (:gve_sma, :gve_eccentricity, :gve_sma_first_index) ||
			throw(ArgumentError("Supported scheduling policies: :gve_sma, :gve_eccentricity, and :gve_sma_first_index."))
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
	elseif policy.policy === :gve_sma_first_index
		return score_gve_sma_first_index!(model, policy, state, mu; candidate_force)
	end
	throw(ArgumentError("Unsupported scheduling policy: $(policy.policy)."))
end

function matching_workspace!(model::InterLinkModel)
	model.matching !== nothing && return model.matching
	terminals = sort!(unique([endpoint for key in model.edge_keys for endpoint in key]))
	indices = Dict(endpoint => index for (index, endpoint) in enumerate(terminals))
	endpoints = Matrix{Int}(undef, 2, length(model.edge_keys))
	edge_order = sortperm(model.edge_keys)
	neighbors = [Int[] for _ in terminals]
	for edge in edge_order, side in 1:2
		vertex = indices[model.edge_keys[edge][side]]
		endpoints[side, edge] = vertex
		push!(neighbors[vertex], edge)
	end
	vertices, edges = length(terminals), length(model.edge_keys)
	model.matching = InterLinkMatchingWorkspace(terminals, endpoints, edge_order, neighbors,
		sizehint!(Int[], edges), sizehint!(Rational{BigInt}[], edges), zeros(BigInt, edges),
		zeros(Int, vertices), sizehint!(Int[], vertices), sizehint!(Int[], vertices),
		zeros(BigInt, vertices), zeros(BigInt, vertices), zeros(Int, vertices),
		fill(false, vertices), fill(false, vertices))
	return model.matching
end

other_endpoint(work, edge, vertex) = work.endpoints[1, edge] == vertex ?
	work.endpoints[2, edge] : work.endpoints[1, edge]

function select_tree!(model, work, first_index, last_index, selected)
	for index in last_index:-1:first_index
		vertex = work.order[index]
		base = BigInt(0)
		for edge in work.neighbors[vertex]
			work.weights[edge] > 0 || continue
			child = other_endpoint(work, edge, vertex)
			work.parent[child] == vertex || continue
			base += work.free_score[child]
		end
		work.blocked_score[vertex] = base
		gain, chosen = BigInt(0), 0
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
		push!(selected, model.edge_keys[edge])
	end
	return nothing
end

# 2.2. Define the function to select the best set of interlinks based on the current scores
"""
	select_interlinks!(model::InterLinkModel)

Activate a maximum-total-score matching of available, positive-score links,
with at most one link per terminal. Cached topology and iterative DFS identify
tree components, solved with linear-time dynamic programming. Cyclic components
use NetworkX's weighted Blossom algorithm through PythonCall, with O(N^3)
arithmetic operations for N cyclic endpoints. Both solvers use exact BigInt
weights scaled from Float64 scores; no floating-point score sums are compared.
Ties are deterministic for a given backend version, but tree ties may choose a
different optimum than Blossom alone. Register topology through register_candidate!;
mutating existing connection states remains supported and is reread on every call.
PythonCall provisions NetworkX from CondaPkg.toml on first use; setup needs network access.
"""
function select_interlinks!(model::InterLinkModel)
	# select the best set of interlinks based on the current scores, then set them to be active links
	# In the selected links will have their interlink_connection.state.active set to true
	work = matching_workspace!(model)
	empty!(work.candidates)
	empty!(work.exact_scores)
	for edge in work.edge_order
		connection = model.connections[edge]
		model.scores[edge] = connection.state.score
		if connection.state.available && connection.state.score > 0.0
			isfinite(connection.state.score) || throw(ArgumentError("Interlink matching scores must be finite."))
			push!(work.candidates, edge)
			push!(work.exact_scores, Rational{BigInt}(connection.state.score))
		end
	end
	selected = InterLinkKey[]
	if !isempty(work.candidates)
		common_denominator = maximum(denominator, work.exact_scores)
		fill!(work.weights, 0)
		for (edge, score) in zip(work.candidates, work.exact_scores)
			work.weights[edge] = numerator(score) * div(common_denominator, denominator(score))
		end
		fill!(work.parent, 0)
		fill!(work.blocked, false)
		fill!(work.cyclic, false)
		empty!(work.order)
		empty!(work.stack)
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
				for index in first_index:last_index
					work.cyclic[work.order[index]] = true
				end
			end
		end
	end
	if !isempty(work.candidates) && any(work.cyclic)
		networkx = PythonCall.pyimport("networkx")
		graph = networkx.Graph()
		for edge in work.candidates
			first_vertex, second_vertex = work.endpoints[1, edge], work.endpoints[2, edge]
			work.cyclic[first_vertex] || continue
			graph.add_edge(first_vertex, second_vertex;
				weight=PythonCall.pybuiltins.int(string(work.weights[edge])))
		end
		matching = PythonCall.pyconvert(Vector{Tuple{Int, Int}},
			PythonCall.pybuiltins.list(networkx.max_weight_matching(graph; maxcardinality=false)))
		for (first_index, second_index) in matching
			push!(selected, minmax(work.terminals[first_index], work.terminals[second_index]))
		end
	end
	for connection in model.connections
		connection.state.active = false
	end
	for key in selected
		model.connections[model.edge_index[key]].state.active = true
	end
	return sort!(selected)
end

# 2.3. Define the function to schedule interlinks based on the current state and scoring policy
function schedule_interlinks!(model::InterLinkModel, policy::SchedulingPolicyModel,
	spacecraft, state, mu::Real; candidate_force, is_active=nothing)
	# step 1. update the availability of each interlink based on the current state
	update_availability!(model, spacecraft, state; is_active)

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
