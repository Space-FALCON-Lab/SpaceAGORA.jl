module InterLinkModels

	using LinearAlgebra
	using StaticArrays

	export InterLinkModel, InterLinkGraph, InterLinkParameters, InterLinkState, InterLinkConnection
	export vertex, InterLinkKey, InterLinkHistory
	export register_candidate!, get_connection, update_availability!

	const vertex = Tuple{Int, Int}
	const InterLinkKey = Tuple{vertex, vertex}

	# 1.1. Candidate Graph
	"""
	Candidate terminal graph with indexed edge properties and state.
	`connections[e].state` is the sole source of availability, active status and score.
	Terminal topology is built lazily. Since `register_candidate!` only appends edges,
	the cache is current when `size(endpoints, 2) == length(edge_keys)`, independently
	of the scheduler workspace. Direct topology edits are unsupported.
	"""
	mutable struct InterLinkGraph
		# i) Core graph data structures
		vertices::Vector{vertex} # 1x2 vertex_id = (satellite_id, terminal_id)
		edges::Matrix{Int} # 2x1 edge_id = (vertex_id; vertex_id)
		
		# ii) Edge properties
		P::Float64 # Power of the interlink
		B::Float64 # Magnification factor of the interlink (for open-cavity laser interlinks)
		range::Float64 # Maximum range of the interlink

		# iii) Edge state
		available::Bool = false # Availability status of the interlink (within_range, not forbidden pair, satellites above 50km)
		active::Bool = false # Active status of the interlink
		score::Float64 = 0.0 # Score of the interlink for scheduler

		# iv) Auxiliary data structures for efficient graph operations (optional)
		vertex_neighbors::Vector{Vector{Int}} # Adjacent edge list for each vertex
		edge_index_lookup::Dict{InterLinkKey, Int} # A lookup from a vertex pair to its edge number (vertex 1, vertex 2) => edge_index
		satellite_edges::Vector{Vector{Int}} # List of edges associated with each satellite (which include multiple terminals / vertices)

		InterLinkGraph() = new(vertex[], Matrix{Int}(undef, 2, 0), 
			0.0, 0.0, 0.0,
			Vector{Vector{Int}}(), InterLinkKey[], Dict{InterLinkKey, Int}(), Vector{Vector{Int}}())
	end

	# 1.2. InterLink History Struct
	struct InterLinkHistory
		time::Float64
		active::Vector{InterLinkKey}
	end

	# 1.3. InterLink Matching Scratch Workspace
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

	# 1.4. InterLink Simulation Model
	mutable struct ConstellationInterLinkModel
		# i) Core model parameters
		constellation_link_type::Symbol
		forbidden_pairs::Set{Tuple{Int, Int}}
		graph::InterLinkGraph

		# ii) Auxiliary data structures for tracking history and dynamics programming matching (optional)
		history::Vector{InterLinkHistory}
		matching::Union{Nothing, InterLinkMatchingWorkspace}

		function ConstellationInterLinkModel(; constellation_link_type::Symbol=:laser, battery_energy_threshold::Real=50.0,
			tempurature_threshold::Real=0.0, forbidden_pairs=Tuple{Int, Int}[])
			constellation_link_type in (:laser, :radio) || throw(ArgumentError("constellation_link_type must be :laser or :radio."))
			pairs = Set{Tuple{Int, Int}}(minmax(pair...) for pair in forbidden_pairs)
			new(constellation_link_type, pairs, InterLinkGraph(), InterLinkHistory[], nothing)
		end
	end

	# 2. Find all possible end pairs (without considering constraaints)
	# function register_candidate!(spacecraft,
	# 	first_endpoint::vertex, second_endpoint::vertex;
		
	# 	# 2.1. Validate the candidate interlink according to the rules.
	# 	# i) rule 1: The interlink must join different spacecraft
	# 	first_endpoint[1] != second_endpoint[1] || throw(ArgumentError("An interlink must join different spacecraft."))
		
	# 	# ii) rule 2: The terminals must be valid for the given spacecraft
	# 	for (satellite, terminal) in (first_endpoint, second_endpoint)
	# 		1 <= satellite <= length(spacecraft) || throw(ArgumentError("Invalid spacecraft vector index $satellite."))
	# 		1 <= terminal <= spacecraft[satellite].n_terminal || throw(ArgumentError("Invalid terminal $terminal on spacecraft $satellite."))
	# 	end

	# 	# iii) rule 3: The key for the candidate pair is the sorted tuple of endpoints
	# 	key = minmax(first_endpoint, second_endpoint)

	# 	# iv) rule 4: The candidate must not already be registered
	# 	graph = InterLinkGraph()  # Assuming a new graph is created for the registration
	# 	haskey(graph.edge_index, key) && throw(ArgumentError("Candidate $key is already registered."))
		
	# 	# 2.2. Register the connection and its key-to-index lookup.
	# 	push!(graph.edges, [graph.vertices .== first_endpoint; graph.vertices .== second_endpoint])
	# 	index = size(graph.edges, 2)
	# 	graph.edge_index[key] = index

	# 	# 2.3. Add the edge index to both spacecraft's adjacency lists.
	# 	while length(graph.satellite_edges) < length(spacecraft)
	# 		push!(graph.satellite_edges, Int[])
	# 	end
	# 	for (satellite, _) in key
	# 		push!(graph.satellite_edges[satellite], index)
	# 	end

	# 	return key 		# Return the key of the registered interlink.
	
	# size(graph.endpoints, 2) == length(graph.edge_keys) && return graph
	# terminals = sort!(unique([endpoint for key in graph.edge_keys for endpoint in key]))
	# indices = Dict(endpoint => index for (index, endpoint) in enumerate(terminals))
	# endpoints = Matrix{Int}(undef, 2, length(graph.edge_keys))
	# neighbors = [Int[] for _ in terminals]
	# for edge in sortperm(graph.edge_keys), side in 1:2
	# 	vertex = indices[graph.edge_keys[edge][side]]
	# 	endpoints[side, edge] = vertex
	# 	push!(neighbors[vertex], edge)
	# end
	# graph.terminals = terminals
	# graph.neighbors = neighbors
	# graph.endpoints = endpoints
	# return graph

	# end

	# 3. Update link availability
	"""
		update_availability!(model::ConstellationInterLinkModel, spacecraft, state; is_active=nothing)

	Check range, forbidden pairs, battery, temperature and optional activity flags.
	Update connection availability and return a fresh `Vector{Bool}`
	in edge-index order; the returned snapshot is not retained by the model.
	"""
	function update_availability!(model::ConstellationInterLinkModel, spacecraft, state; is_active=nothing)
		results = Vector{Bool}(undef, length(model.graph.connections))
		for index in eachindex(model.graph.connections)
			key = model.graph.edge_keys[index]
			connection = model.graph.connections[index]
			first, second = key[1][1], key[2][1]

			# i) rule 1 & 2: The interlink must be within range & not forbidden
			separation = norm(SVector{3, Float64}(state.sc[first].pos) -
				SVector{3, Float64}(state.sc[second].pos))
			available = (0 < separation <= connection.parameters.range) &
				!(minmax(first, second) in model.forbidden_pairs)

			# ii) rule 3: The spacecraft must be active when activity flags are supplied (not below 50 km)
			for satellite in (first, second)
				available &= (is_active === nothing || is_active[satellite])
			end

			connection.state.available = available
			results[index] = available
		end
		return results
	end

end
