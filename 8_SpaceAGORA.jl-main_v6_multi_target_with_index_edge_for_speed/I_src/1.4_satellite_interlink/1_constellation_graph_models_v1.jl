module InterLinkModels

using LinearAlgebra
using StaticArrays

export InterLinkModel, InterLinkGraph, InterLinkParameters, InterLinkState, InterLinkConnection
export TerminalEndpoint, InterLinkKey, InterLinkScheduleSample
export register_candidate!, get_connection, update_availability!

const TerminalEndpoint = Tuple{Int, Int}
const InterLinkKey = Tuple{TerminalEndpoint, TerminalEndpoint}

# 1.1. Interlink Property & State Struct
struct InterLinkParameters
	P::Float64
	B::Float64
	range::Float64

	function InterLinkParameters(; P::Real=10_000.0, B::Real=100.0, range::Real=200e3)
		isfinite(P) && P >= 0 || throw(ArgumentError("P must be finite and nonnegative."))
		isfinite(B) && B > 0 || throw(ArgumentError("B must be finite and positive."))
		isfinite(range) && range > 0 || throw(ArgumentError("range must be finite and positive."))
		new(Float64(P), Float64(B), Float64(range))
	end
end

Base.@kwdef mutable struct InterLinkState
	available::Bool = false
	active::Bool = false
	score::Float64 = 0.0
end

struct InterLinkConnection
	parameters::InterLinkParameters
	state::InterLinkState
end

# 1.2. Interlink History Struct
struct InterLinkScheduleSample
	time::Float64
	active::Vector{InterLinkKey}
end

# 1.3. Matching Scratch Workspace
struct InterLinkMatchingWorkspace
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

# 1.4. Candidate Graph
"""
Candidate terminal graph with indexed edge properties and state.
`connections[e].state` is the sole source of availability, active status and score.
Terminal topology is built lazily. Since `register_candidate!` only appends edges,
the cache is current when `size(endpoints, 2) == length(edge_keys)`, independently
of the scheduler workspace. Direct topology edits are unsupported.
"""
mutable struct InterLinkGraph
	connections::Vector{InterLinkConnection}

	edge_keys::Vector{InterLinkKey}
	edge_index::Dict{InterLinkKey, Int}

	satellite_edges::Vector{Vector{Int}}
	terminals::Vector{TerminalEndpoint}
	endpoints::Matrix{Int}
	neighbors::Vector{Vector{Int}}

	InterLinkGraph() = new(InterLinkConnection[], InterLinkKey[], Dict{InterLinkKey, Int}(),
		Vector{Int}[], TerminalEndpoint[], Matrix{Int}(undef, 2, 0), Vector{Int}[])
end

function ensure_topology!(graph::InterLinkGraph)
	size(graph.endpoints, 2) == length(graph.edge_keys) && return graph
	terminals = sort!(unique([endpoint for key in graph.edge_keys for endpoint in key]))
	indices = Dict(endpoint => index for (index, endpoint) in enumerate(terminals))
	endpoints = Matrix{Int}(undef, 2, length(graph.edge_keys))
	neighbors = [Int[] for _ in terminals]
	for edge in sortperm(graph.edge_keys), side in 1:2
		vertex = indices[graph.edge_keys[edge][side]]
		endpoints[side, edge] = vertex
		push!(neighbors[vertex], edge)
	end
	graph.terminals = terminals
	graph.neighbors = neighbors
	graph.endpoints = endpoints
	return graph
end

# 1.5. Interlink Simulation Model
"""
Interlink configuration, candidate graph, scheduling history and solver scratch space.
Use `get_connection(model, key)` for keyed access to the graph's authoritative connection.
Registering candidates invalidates graph topology and discards matching scratch storage.
"""
mutable struct InterLinkModel
	link_type::Symbol
	battery_energy_threshold::Float64
	tempurature_threshold::Float64
	forbidden_pairs::Set{Tuple{Int, Int}}
	graph::InterLinkGraph
	history::Vector{InterLinkScheduleSample}
	matching::Union{Nothing, InterLinkMatchingWorkspace}

	function InterLinkModel(; link_type::Symbol=:laser, battery_energy_threshold::Real=50.0,
		tempurature_threshold::Real=0.0, forbidden_pairs=Tuple{Int, Int}[])
		link_type in (:laser, :radio) || throw(ArgumentError("link_type must be :laser or :radio."))
		0 <= battery_energy_threshold <= 100 || throw(ArgumentError("battery_energy_threshold must be in [0, 100]."))
		0 <= tempurature_threshold <= 100 || throw(ArgumentError("tempurature_threshold must be in [0, 100]."))
		pairs = Set{Tuple{Int, Int}}(minmax(pair...) for pair in forbidden_pairs)
		new(link_type, Float64(battery_energy_threshold), Float64(tempurature_threshold),
			pairs, InterLinkGraph(), InterLinkScheduleSample[], nothing)
	end
end

"""
	get_connection(model::InterLinkModel, key::InterLinkKey)

Return the registered connection for a canonical endpoint key, or throw `KeyError`.
Mutating its state is supported; do not edit graph topology arrays directly.
"""
get_connection(model::InterLinkModel, key::InterLinkKey) =
	model.graph.connections[model.graph.edge_index[key]]

# 2. Find all possible end pairs (without considering constraaints)
function register_candidate!(model::InterLinkModel, spacecraft,
	first_endpoint::TerminalEndpoint, second_endpoint::TerminalEndpoint;
	parameters::InterLinkParameters=InterLinkParameters())
	
	# i) rule 1: The interlink must join different spacecraft
	first_endpoint[1] != second_endpoint[1] || throw(ArgumentError("An interlink must join different spacecraft."))
	
	# ii) rule 2: The terminals must be valid for the given spacecraft
	for (satellite, terminal) in (first_endpoint, second_endpoint)
		1 <= satellite <= length(spacecraft) || throw(ArgumentError("Invalid spacecraft vector index $satellite."))
		1 <= terminal <= spacecraft[satellite].n_terminal || throw(ArgumentError("Invalid terminal $terminal on spacecraft $satellite."))
	end

	# iii) rule 3: The key for the candidate pair is the sorted tuple of endpoints
	key = minmax(first_endpoint, second_endpoint)

	# iv) rule 4: The candidate must not already be registered
	graph = model.graph
	haskey(graph.edge_index, key) && throw(ArgumentError("Candidate $key is already registered."))
	connection = InterLinkConnection(parameters, InterLinkState())
	push!(graph.edge_keys, key)
	push!(graph.connections, connection)
	index = length(graph.edge_keys)
	graph.edge_index[key] = index
	while length(graph.satellite_edges) < length(spacecraft)
		push!(graph.satellite_edges, Int[])
	end
	for (satellite, _) in key
		push!(graph.satellite_edges[satellite], index)
	end
	model.matching = nothing
	return key
end

# 3. Update availability for one candidate pair
function update_availability!(model::InterLinkModel, key::InterLinkKey, spacecraft, state; is_active=nothing)
	return _update_availability!(model, model.graph.edge_index[key], spacecraft, state, is_active)
end

"""
	update_availability!(model::InterLinkModel, spacecraft, state; is_active=nothing)

Refresh all connection availability states. Return a fresh `Vector{Bool}`
in edge-index order; the returned snapshot is not retained by the model.
"""
function update_availability!(model::InterLinkModel, spacecraft, state; is_active=nothing)
	return Bool[_update_availability!(model, index, spacecraft, state, is_active)
		for index in eachindex(model.graph.connections)]
end

function _update_availability!(model::InterLinkModel, index::Int, spacecraft, state, is_active)
	key = model.graph.edge_keys[index]
	connection = model.graph.connections[index]
	first_satellite, second_satellite = key[1][1], key[2][1]

	# i) rule 1: The interlink must be within range
	separation = norm(SVector{3, Float64}(state.sc[first_satellite].pos) -
		SVector{3, Float64}(state.sc[second_satellite].pos))
	available = 0 < separation <= connection.parameters.range

	# ii) rule 2: The interlink must not be forbidden
	available &= !(minmax(first_satellite, second_satellite) in model.forbidden_pairs)

	for (satellite, terminal) in key
		vehicle = spacecraft[satellite]

		# iii) rule 3: The vehicle must have sufficient battery energy and temperature levels
		available &= vehicle.battery_energy_index >= model.battery_energy_threshold
		available &= vehicle.tempurature_index >= model.tempurature_threshold

		# iv) rule 4: The spacecraft must be active when activity flags are supplied (not below 50 km)
		available &= (is_active === nothing || is_active[satellite])
	end

	connection.state.available = available # Update availability status of the interlink connection
	return available
end

end
