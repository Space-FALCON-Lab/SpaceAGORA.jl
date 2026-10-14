module SchedulingPolicies

using LinearAlgebra
using DiffEqBase
import PythonCall
using ..InterLinkModels

export SchedulingPolicyModel, score_candidates!
export select_interlinks!, schedule_interlinks!, interlink_scheduler_callback

# 1. Matching workspace for dynamic programming
# 1.1. Scheduler-owned scores and matching scratch (reused between scheduling calls)
struct InterLinkMatchingWorkspace
	graph::InterLinkGraph
	score::Vector{Float64}
	candidates::Vector{Int}
	exact_scores::Vector{Rational{BigInt}}
	weights::Vector{BigInt}
	parent::Vector{Int}
	order::Vector{Int}
	stack::Vector{Int}
	free_score::Vector{BigInt}
	blocked_score::Vector{BigInt}
	choice::Vector{Int}
	blocked::Vector{Bool}
	cyclic::Vector{Bool}
end

# 2. Scheduling policy for interlink selection
# 2.1. Scheduling policy model
"""
	SchedulingPolicyModel(policy; target_idx=1, target_idxs=nothing)

Create a scheduling policy. For `:gve_sma` and `:gve_eccentricity`, `target_idx`
selects one target and `target_idx=nothing` scores both link endpoints.
For `:gve_sma`, `target_idxs` optionally selects a set of targets and takes
precedence over `target_idx`; all selected indices must exist in the simulation state.
"""
mutable struct SchedulingPolicyModel
	const policy::Symbol
	const target_idx::Union{Nothing, Int}
	const target_idxs::Union{Nothing, Vector{Int}}
	matching::Union{Nothing, InterLinkMatchingWorkspace}

	function SchedulingPolicyModel(policy::Symbol; target_idx::Union{Nothing, Integer}=1,
		target_idxs::Union{Nothing, AbstractVector{<:Integer}}=nothing)
		policy in (:gve_sma, :gve_eccentricity, :gve_sma_first_index) ||
			throw(ArgumentError("Supported scheduling policies: :gve_sma, :gve_eccentricity, and :gve_sma_first_index."))
		target_idxs === nothing || policy === :gve_sma ||
			throw(ArgumentError("target_idxs is supported only by the :gve_sma policy."))
		target_idx === nothing || target_idx > 0 ||
			throw(ArgumentError("target_idx must be nothing or a positive spacecraft vector index."))
		target_indices = target_idxs === nothing ? nothing : Int[target_idxs...]
		if target_indices !== nothing
			isempty(target_indices) && throw(ArgumentError("target_idxs must not be empty."))
			all(>(0), target_indices) || throw(ArgumentError("target_idxs must contain only positive indices."))
			length(unique(target_indices)) == length(target_indices) ||
				throw(ArgumentError("target_idxs must not contain duplicates."))
			sort!(target_indices)
		end
		new(policy, target_idx === nothing ? nothing : Int(target_idx), target_indices, nothing)
	end
end

# 1.2. Create and maintain the dynamic-programming workspace.
function matching_workspace!(model::ConstellationInterLinkModel, policy::SchedulingPolicyModel)
	graph = model.graph
	vertices, edges = length(graph.vertices), length(graph.edges)
	work = policy.matching
	if work !== nothing && work.graph === graph &&
		length(work.parent) == vertices && length(work.weights) == edges
		return work
	end
	scores = zeros(Float64, edges)
	if work !== nothing && work.graph === graph
		copyto!(scores, 1, work.score, 1, min(edges, length(work.score)))
	end
	policy.matching = InterLinkMatchingWorkspace(
		graph, scores, sizehint!(Int[], edges), sizehint!(Rational{BigInt}[], edges), zeros(BigInt, edges),
		zeros(Int, vertices), sizehint!(Int[], vertices), sizehint!(Int[], vertices),
		zeros(BigInt, vertices), zeros(BigInt, vertices), zeros(Int, vertices),
		fill(false, vertices), fill(false, vertices))
	return policy.matching
end

include(joinpath(@__DIR__, "3_scheduling_policies.jl"))

# 2.2. Define the function to score candidate interlinks based on the scheduling policy
function score_candidates!(model::ConstellationInterLinkModel, policy::SchedulingPolicyModel, state, mu::Real;
	# compute the score of each interlink based on the scheduling policy / scoring policy
	# resulting scores will be stored in the policy's workspace
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

# 3. Dynamic programming:  Select the best set of interlinks based on the current scores
# Helper Function: given one vertex and a related edge, this helper function of find the other vertex of the edge
other_endpoint(model, edge, vertex) = model.graph.edges[edge][1] == vertex ?
	model.graph.edges[edge][2] : model.graph.edges[edge][1]

# Helper function: select the best set of interlinks in a tree structure
function select_tree!(model, work, first_index, last_index, selected)
	for index in last_index:-1:first_index
		vertex = work.order[index]
		base = BigInt(0)
		for edge in model.graph.vertex_neighbors[vertex]
			work.weights[edge] > 0 || continue
			child = other_endpoint(model, edge, vertex)
			work.parent[child] == vertex || continue
			base += work.free_score[child]
		end
		work.blocked_score[vertex] = base
		gain, chosen = BigInt(0), 0
		for edge in model.graph.vertex_neighbors[vertex]
			work.weights[edge] > 0 || continue
			child = other_endpoint(model, edge, vertex)
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
		work.blocked[other_endpoint(model, edge, vertex)] = true
		push!(selected, map(v -> model.graph.vertices[v], model.graph.edges[edge]))
	end
	return nothing
end


"""
	select_interlinks!(model::ConstellationInterLinkModel, policy::SchedulingPolicyModel)

Activate a maximum-total-score matching of available, positive-score links,
with at most one link per terminal. Topology belongs to `model.graph`; `policy.matching`
contains raw per-edge `score` values and reusable solver scratch arrays. Pass the same
policy used by `score_candidates!`; selection uses its stored scores, not its scoring
objective. A fresh workspace starts with zero scores. Appending links preserves existing
scores and initializes new scores to zero when the workspace is refreshed. Switching
graphs or discarding the workspace resets scores; rescore before selecting again.
Use separate policies for independent scheduling runs. Deep-copy the model and policy
together to preserve their workspace binding. Edge order is sorted locally on each
selection, not stored in either structure. Registered topology and iterative DFS
identify tree components, solved with linear-time dynamic programming. Cyclic components
use NetworkX's weighted Blossom algorithm through PythonCall, with O(N^3)
arithmetic operations for N cyclic endpoints. Both solvers use exact BigInt
weights scaled from Float64 scores; no floating-point score sums are compared.
Ties are deterministic for a given backend version, but tree ties may choose a
different optimum than Blossom alone. Register topology through `register_allowed_link!`.
For manual scores, refresh with `matching_workspace!(model, policy)` before editing
its `score` vector. Graph availability and workspace scores are reread on every call.
PythonCall provisions NetworkX from CondaPkg.toml on first use; setup needs network access.
"""
function select_interlinks!(model::ConstellationInterLinkModel, policy::SchedulingPolicyModel)
	# Step 1. Collect available, positive-score edges in canonical endpoint order.
	work = matching_workspace!(model, policy)
	empty!(work.candidates)
	empty!(work.exact_scores)
	edge_order = sortperm([map(v -> model.graph.vertices[v], endpoints) for endpoints in model.graph.edges])
	for edge in edge_order
		if model.graph.available[edge] && work.score[edge] > 0.0
			isfinite(work.score[edge]) || throw(ArgumentError("Interlink matching scores must be finite."))
			push!(work.candidates, edge)
			push!(work.exact_scores, Rational{BigInt}(work.score[edge]))
		end
	end
	selected = InterLinkKey[]
	if !isempty(work.candidates)
		# Step 2. Scale exact Float64 rationals to integer weights without rounding score sums.
		# Their denominators are powers of two, so the largest is a common denominator.
		common_denominator = maximum(denominator, work.exact_scores)
		fill!(work.weights, 0)
		for (edge, score) in zip(work.candidates, work.exact_scores)
			work.weights[edge] = numerator(score) * div(common_denominator, denominator(score))
		end
		# Step 3. Traverse positive-weight components with depth-first search; solve trees and mark cycles if detected.
		fill!(work.parent, 0)
		fill!(work.blocked, false)
		fill!(work.cyclic, false)
		empty!(work.order)
		empty!(work.stack)
		for root in sortperm(model.graph.vertices)
			work.parent[root] == 0 || continue
			work.parent[root] = -1
			first_index = length(work.order) + 1
			degree_sum = 0
			push!(work.stack, root)
			while !isempty(work.stack)
				vertex = pop!(work.stack)
				push!(work.order, vertex)
				for edge in model.graph.vertex_neighbors[vertex]
					work.weights[edge] > 0 || continue
					degree_sum += 1
					partner = other_endpoint(model, edge, vertex)
					work.parent[partner] == 0 || continue
					work.parent[partner] = vertex
					push!(work.stack, partner)
				end
			end
			last_index = length(work.order)
			degree_sum == 0 && continue
			# Each undirected edge is counted twice; a connected tree has E = V - 1.
			if degree_sum ÷ 2 == last_index - first_index
				select_tree!(model, work, first_index, last_index, selected)
			else
				for index in first_index:last_index
					work.cyclic[work.order[index]] = true
				end
			end
		end
	end
	# Step 4. Solve cyclic components with NetworkX's maximum-weight Blossom matching.
	if !isempty(work.candidates) && any(work.cyclic)
		PythonCall.C.PyGILState_Check() != 0 || throw(ArgumentError(
			"Cyclic interlink matching requires the Python GIL. Run on Python's owning thread or use process workers. " *
			"For threaded calls, release the GIL on its owner and acquire it with PythonCall.GIL.@lock."))
		networkx = PythonCall.pyimport("networkx")
		graph = networkx.Graph()
		for edge in work.candidates
			first_vertex, second_vertex = model.graph.edges[edge][1], model.graph.edges[edge][2]
			work.cyclic[first_vertex] || continue
			graph.add_edge(first_vertex, second_vertex;
				weight=PythonCall.pybuiltins.int(string(work.weights[edge])))
		end
		matching = PythonCall.pyconvert(Vector{Tuple{Int, Int}},
			PythonCall.pybuiltins.list(networkx.max_weight_matching(graph; maxcardinality=false)))
		for (first_index, second_index) in matching
			push!(selected, minmax(model.graph.vertices[first_index], model.graph.vertices[second_index]))
		end
	end
	# Step 5. Replace the active-link flags with the combined tree and cyclic selections.
	fill!(model.graph.active, false)
	for key in selected
		model.graph.active[model.graph.edge_index_lookup[key]] = true
	end
	return sort!(selected)
end

# 4. Overall scheduling function to schedule interlinks based on the current state and scoring policy
function schedule_interlinks!(model::ConstellationInterLinkModel, policy::SchedulingPolicyModel,
	spacecraft, state, mu::Real; candidate_force, is_active=nothing)
	# step 1. update the availability of each interlink based on the current state
	update_availability!(model, spacecraft, state; is_active)

	# step 2. compute the score of each interlink based on the scheduling policy / scoring policy
	score_candidates!(model, policy, state, mu; candidate_force) 

	# step 3. select the best set of interlinks based on the current scores
	return select_interlinks!(model, policy)
end

# 5. Define the callback function for the interlink graph update with scheduler
function interlink_scheduler_callback(; candidate_force)
	function affect!(integrator)
		args = integrator.p.args
		model = args.interlink_model
		selected = schedule_interlinks!(model, args.scheduling_policy_model,
			args.dynamics_model.spacecraft, integrator.u, args.environment_model.planet.μ;
			candidate_force, is_active=integrator.p.is_active)
		sample = InterLinkHistory(Float64(integrator.t), selected)
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
