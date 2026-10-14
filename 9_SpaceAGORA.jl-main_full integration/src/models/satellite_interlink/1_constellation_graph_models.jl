module InterLinkModels

	using LinearAlgebra
	using StaticArrays

	export ConstellationInterLinkModel, InterLinkGraph, InterLinkParameters
	export vertex, InterLinkKey, InterLinkHistory
	export register_allowed_link!, update_availability!

	const vertex = Tuple{Int, Int} # (spacecraft vector index, terminal number)
	const InterLinkKey = Tuple{vertex, vertex} # Canonically ordered pair of terminal endpoints

	# 1.1 Candidate graph (Vertices are satellite terminals; Edges are allowed interlinks)
	struct InterLinkParameters
		P::Float64 # Interlink power [W]
		B::Float64 # Magnification factor for open-cavity laser interlinks
		range::Float64 # Maximum interlink range [m]

		function InterLinkParameters(; P::Real=10_000.0, B::Real=100.0, range::Real=200e3)
			isfinite(P) && P >= 0 || throw(ArgumentError("P must be finite and nonnegative."))
			isfinite(B) && B > 0 || throw(ArgumentError("B must be finite and positive."))
			isfinite(range) && range > 0 || throw(ArgumentError("range must be finite and positive."))
			new(Float64(P), Float64(B), Float64(range))
		end
	end

	mutable struct InterLinkGraph
		# i) Core graph data structures
		vertices::Vector{vertex} # vertices[v] = (satellite, terminal)
		edges::Vector{Tuple{Int, Int}} # edges[e] = (first vertex index, second vertex index)

		# ii) Edge properties (one entry per edge)
		edge_properties::Vector{InterLinkParameters}

		# iii) Edge state (availability does not imply scheduler selection)
		available::Vector{Bool} # Within range, not forbidden, and both spacecraft active if flags are supplied
		active::Vector{Bool} # Links selected by the scheduler

		# iv) Auxiliary indexes for efficient graph operations
		vertex_neighbors::Vector{Vector{Int}} # Incident edge indexes for each terminal vertex
		edge_index_lookup::Dict{InterLinkKey, Int} # Canonical endpoint pair => edge index
		satellite_edges::Vector{Vector{Int}} # Incident edge indexes across all terminals of each spacecraft

		InterLinkGraph() = new(vertex[], Tuple{Int, Int}[], InterLinkParameters[],
			Bool[], Bool[], Vector{Int}[], Dict{InterLinkKey, Int}(), Vector{Int}[])
	end

	# 1.2. Interlink history
	struct InterLinkHistory
		time::Float64
		active::Vector{InterLinkKey}
	end

	# 1.3. Constellation interlink model
	mutable struct ConstellationInterLinkModel
		# i) Core model parameters
		constellation_link_type::Symbol
		forbidden_pairs::Set{Tuple{Int, Int}}
		graph::InterLinkGraph

		# ii) Scheduling history (independent of the selection algorithm)
		history::Vector{InterLinkHistory}

		function ConstellationInterLinkModel(; constellation_link_type::Symbol=:laser,
			forbidden_pairs=Tuple{Int, Int}[])
			constellation_link_type in (:laser, :radio) || throw(ArgumentError("constellation_link_type must be :laser or :radio."))
			pairs = Set{Tuple{Int, Int}}(minmax(pair...) for pair in forbidden_pairs)
			new(constellation_link_type, pairs, InterLinkGraph(), InterLinkHistory[])
		end
	end

	# 2. Register a allowed pair as for vertices and edges in graph without checking its current availability
	function register_allowed_link!(model::ConstellationInterLinkModel, spacecraft,
		first_endpoint::vertex, second_endpoint::vertex;
		parameters::InterLinkParameters=InterLinkParameters())

		# 2.1. Validate the terminals before changing the graph.
		# i) Rule 1: The interlink must join different spacecraft.
		first_endpoint[1] != second_endpoint[1] || throw(ArgumentError("An interlink must join different spacecraft."))

		# ii) Rule 2: Both terminal numbers must be valid for their spacecraft.
		for (satellite, terminal) in (first_endpoint, second_endpoint)
			1 <= satellite <= length(spacecraft) || throw(ArgumentError("Invalid spacecraft vector index $satellite."))
			1 <= terminal <= spacecraft[satellite].n_terminal || throw(ArgumentError("Invalid terminal $terminal on spacecraft $satellite."))
		end

		# iii) Rule 3: Sort endpoints so reversing a pair produces the same key.
		key = minmax(first_endpoint, second_endpoint)

		# iv) Rule 4: The candidate must not already be registered.
		graph = model.graph
		haskey(graph.edge_index_lookup, key) && throw(ArgumentError("Candidate $key is already registered."))

		# 2.2. Connect terminal-pairs, and give the new link its own parameters and initial state.
		# i) Locate the terminals, put them into the graph's vertex list
		indices = map(key) do endpoint
			index = findfirst(==(endpoint), graph.vertices)
			if index === nothing
				push!(graph.vertices, endpoint)
				push!(graph.vertex_neighbors, Int[])
				index = length(graph.vertices)
			end
			index
		end

		# ii) Add the new edge (interlink) to the graph's edge list (interlink list) and initialize the properties and state of that edge
		push!(graph.edges, indices)
		push!(graph.edge_properties, parameters)
		push!(graph.available, false)
		push!(graph.active, false)
		edge = length(graph.edges)

		# iii) Record the lookup
		graph.edge_index_lookup[key] = edge

		# 2.3. Keep terminal adjacency sorted by endpoint key for deterministic matching.
		for index in indices
			neighbors = graph.vertex_neighbors[index]
			position = searchsortedfirst(neighbors, edge;
				by=e -> map(v -> graph.vertices[v], graph.edges[e]))
			insert!(neighbors, position, edge)
		end

		# 2.4. Add the edge to both spacecraft's adjacency lists, including empty lists for unlinked spacecraft.
		while length(graph.satellite_edges) < length(spacecraft)
			push!(graph.satellite_edges, Int[])
		end
		for (satellite, _) in key
			push!(graph.satellite_edges[satellite], edge)
		end

		return key
	end

	# 3. Update link availability
	"""
		update_availability!(model::ConstellationInterLinkModel, spacecraft, state; is_active=nothing)

	Check positive separation, per-edge range, forbidden pairs and optional spacecraft
activity flags. Battery and temperature do not gate availability.
	Return a fresh `Vector{Bool}` in registration order.
	"""
	function update_availability!(model::ConstellationInterLinkModel, spacecraft, state; is_active=nothing)
		graph = model.graph
		for edge in eachindex(graph.edges)
			first, second = map(v -> graph.vertices[v], graph.edges[edge])

			# i) Rule 1: Separation must be positive and within this edge's range.
			separation = norm(SVector{3, Float64}(state.sc[first[1]].pos) -
				SVector{3, Float64}(state.sc[second[1]].pos))
			graph.available[edge] = (0 < separation <= graph.edge_properties[edge].range) &&

				# ii) Rule 2: The spacecraft pair must not be forbidden.
				!((first[1], second[1]) in model.forbidden_pairs) &&

				# iii) Rule 3: Both spacecraft must be active when activity flags are supplied.
				(is_active === nothing || (is_active[first[1]] && is_active[second[1]]))
		end
		return copy(graph.available)
	end

end
