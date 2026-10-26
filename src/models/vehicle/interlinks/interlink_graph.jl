module InterLinkModels

	using LinearAlgebra
	using StaticArrays
	using ..Components: LaserTerminal, PowerReceiver, RadioTerminal, _SPEED_OF_LIGHT_MPS
	using ..AbstractTypes: AbstractActuator
	using ..EffectorSampling: StateSample
	using ..HookDeprecation: _hook_deprecation, _NOT_PASSED
	using ..Geodesy: geodetic_to_pp
	using ..QuaternionMath: rot

	export AbstractInterlink, LinkContext, LaserMomentumLink, link_range_m, link_available, link_force
	export interlink_wrench, link_applies_torque, terminal_pose
	export link_directed, link_terminal_levels, link_terminal, terminal_supports, terminal_shed_tier, terminal_count
	export register_links!, all_pairs, same_plane_pairs, ring_pairs, grid_pairs
	export PowerBeamLink, link_received_power_w, beam_spot_radius_m, beam_capture_fraction
	export LaserBeamLink, AbstractLaserBeam, GaussianBeam, TerminalPointing, beam_force_factor
	export DataLink, AbstractDataBudget, FixedRateBudget, ShannonBudget, link_data_rate_bps, budget_rate_bps
	export has_data_links, contact_windows, contact_time_s, available_links
	export InterlinkModel, InterLinkGraph, register_link!, update_availability!, interlink_force_on_spacecraft
	export vertex, InterLinkKey, InterLinkHistory
	export ConstellationInterLinkModel, InterLinkParameters, register_allowed_link!
	export GroundStation, add_ground_station!, station_pairs

	include(joinpath(@__DIR__, "links.jl"))
	include(joinpath(@__DIR__, "power_links.jl"))
	include(joinpath(@__DIR__, "data_links.jl"))
	include(joinpath(@__DIR__, "ground_stations.jl"))
	include(joinpath(@__DIR__, "laser_beam_links.jl"))

	const vertex = Tuple{Int, Int} # (node, terminal number): a spacecraft vector index, or minus a ground station's number
	const InterLinkKey = Tuple{vertex, vertex} # Canonically ordered pair of terminal endpoints

	"""
	    InterLinkParameters

	Deprecated, use `LaserMomentumLink`: the same type under its old name, which
	will be removed in SpaceAGORA 1.0. The name is an alias and cannot warn on its
	own: `InterLinkParameters()`, `InterLinkParameters(terminal)` and the positional
	`InterLinkParameters(P, B, range)` run silently. Its keywords `P`, `B` and
	`range` and the fields of those names (now `power_w`, `magnification` and
	`range_m`) warn.
	"""
	const InterLinkParameters = LaserMomentumLink

	# 1.1 Candidate graph (Vertices are satellite terminals; Edges are allowed interlinks)
	mutable struct InterLinkGraph
		# i) Core graph data structures
		vertices::Vector{vertex} # vertices[v] = (satellite, terminal)
		edges::Vector{Tuple{Int, Int}} # edges[e] = (first vertex index, second vertex index)

		# ii) The link of each edge (its kind and parameters)
		edge_properties::Vector{AbstractInterlink}

		# iii) Edge state (availability does not imply scheduler selection)
		available::Vector{Bool} # Within range, not forbidden, and both spacecraft active if flags are supplied
		active::Vector{Bool} # Links selected by the scheduler

		# iv) Auxiliary indexes for efficient graph operations
		vertex_neighbors::Vector{Vector{Int}} # Incident edge indexes for each terminal vertex
		edge_index_lookup::Dict{InterLinkKey, Int} # Canonical endpoint pair => edge index
		satellite_edges::Vector{Vector{Int}} # Incident edge indexes across all terminals of each spacecraft

		InterLinkGraph() = new(vertex[], Tuple{Int, Int}[], AbstractInterlink[],
			Bool[], Bool[], Vector{Int}[], Dict{InterLinkKey, Int}(), Vector{Int}[])
	end

	# 1.2. Interlink history: at each accepted step, the links the scheduler
	# selected, the score that selected each, and, when the model records them,
	# the links that were available (could close, `record_available`) and the
	# scores of a weighted mix's terms (`record_terms`). A sample built from a time
	# and the selected links records neither scores nor available links.
	struct InterLinkHistory
		time::Float64
		active::Vector{InterLinkKey}
		available::Vector{InterLinkKey}
		scores::Vector{Float64}   # the score of each active link, in `active` order; empty when not recorded
		terms::Matrix{Float64}    # each top-level term's score (rows) for each active link (columns); no rows when not recorded

		function InterLinkHistory(time::Real, active::AbstractVector, available::AbstractVector=InterLinkKey[],
			scores::AbstractVector=Float64[], terms::AbstractMatrix=zeros(0, length(active)))
			isempty(scores) || length(scores) == length(active) || throw(ArgumentError(
				"An InterLinkHistory sample needs one score for each active link, or none; got $(length(scores)) " *
				"scores for $(length(active)) links."))
			size(terms, 2) == length(active) || throw(ArgumentError(
				"An InterLinkHistory sample's terms need one column for each active link; got $(size(terms, 2)) " *
				"columns for $(length(active)) links."))
			new(Float64(time), convert(Vector{InterLinkKey}, active), convert(Vector{InterLinkKey}, available),
				convert(Vector{Float64}, scores), convert(Matrix{Float64}, terms))
		end
	end

	# 1.3. Interlink model
	"""
	    InterlinkModel(; forbidden_pairs=(), record_available=nothing, record_terms=false)

	The interlinks of a run: the candidate links between spacecraft terminals,
	each with its own link kind (`register_link!`), spacecraft pairs that may never
	link, and the accepted-step history of the links the scheduler selected. Attach
	it with `SimulationConfiguration(interlink_model=model, ...)`. Endpoints are
	`(spacecraft vector index, terminal number)` pairs: the position of the
	spacecraft in `dynamics_model.spacecraft`, not its `id`. Terminals
	`1:n_terminal` are virtual; mounted terminal hardware (see `link_terminal`)
	follows. Ground stations (`add_ground_station!`) are link ends too: antenna
	`k` of station `s` is the endpoint `(-s, k)`. `forbidden_pairs` lists pairs of
	spacecraft indices that may never link; `(i, -s)` forbids spacecraft `i` and
	station `s`.

	`record_available` sets whether each history sample also records the links
	that were available, which `contact_windows` reads: `nothing`, the default,
	records them in a run with a `DataLink`; `true` always; `false` never. A large
	constellation over a long run keeps a list of available links per accepted
	step, so the default leaves them out of runs that carry no data.

	Each history sample also records the score that selected each active link.
	`record_terms=true` records as well, for a `WeightedPolicy`, each of its
	top-level terms' scores for each active link; a run with any other policy
	refuses it. `link_schedule` turns the history into a table.
	"""
	mutable struct InterlinkModel
		# i) Core model parameters
		forbidden_pairs::Set{Tuple{Int, Int}}
		graph::InterLinkGraph

		# ii) Scheduling history (independent of the selection algorithm)
		history::Vector{InterLinkHistory}

		# iii) The run's hardware link endpoints, set at run setup: for each spacecraft,
		# `(edge, side, mounted actuator position)` of each terminal hardware endpoint.
		terminal_endpoints::Vector{Vector{NTuple{3, Int}}}

		# iv) Whether the history records the available links: `nothing` for runs with a data link.
		record_available::Union{Nothing, Bool}

		# v) Whether the history records a weighted policy's term scores.
		record_terms::Bool

		# vi) The ground stations, in registration order: station `s` is node `-s`.
		ground_stations::Vector{GroundStation}
		# vii) The run's terminal hardware at each edge's two ends (`nothing` for a
		# virtual terminal), set at run setup with `terminal_endpoints`; empty before.
		edge_terminals::Vector{Tuple{Any, Any}}

		function InterlinkModel(; forbidden_pairs=Tuple{Int, Int}[], record_available::Union{Nothing, Bool}=nothing,
			record_terms::Bool=false, constellation_link_type=_NOT_PASSED)
			if constellation_link_type !== _NOT_PASSED
				_interlink_deprecation((InterlinkModel, :constellation_link_type),
					"The constellation_link_type keyword is deprecated and will be removed in SpaceAGORA 1.0: " *
					"each registered link now has its own kind (register_link!).")
				constellation_link_type === :laser || throw(ArgumentError(
					"constellation_link_type=$(repr(constellation_link_type)) is not supported: only laser momentum " *
					"links (LaserMomentumLink) are implemented."))
			end
			pairs = Set{Tuple{Int, Int}}(minmax(pair...) for pair in forbidden_pairs)
			new(pairs, InterLinkGraph(), InterLinkHistory[], Vector{NTuple{3, Int}}[], record_available, record_terms,
				GroundStation[], Tuple{Any, Any}[])
		end
	end

	"""
	    ConstellationInterLinkModel

	Deprecated, use `InterlinkModel`: the same type under its old name, which will
	be removed in SpaceAGORA 1.0. The name is an alias and cannot warn on its own:
	`ConstellationInterLinkModel(; forbidden_pairs)` runs silently. Its
	`constellation_link_type` keyword and property warn: each registered link has
	its own kind.
	"""
	const ConstellationInterLinkModel = InterlinkModel

	# `model.constellation_link_type`: every link of an `InterlinkModel` is a laser
	# momentum link or a custom kind; the old property reads `:laser`.
	@inline function Base.getproperty(model::InterlinkModel, name::Symbol)
		if name === :constellation_link_type
			_interlink_deprecation((InterlinkModel, :constellation_link_type),
				"InterlinkModel.constellation_link_type is deprecated and will be removed in SpaceAGORA 1.0: " *
				"each registered link has its own kind (graph.edge_properties).")
			return :laser
		end
		return getfield(model, name)
	end

	# 2. Register a allowed pair as for vertices and edges in graph without checking its current availability
	"""
	    register_link!(model::InterlinkModel, spacecraft, first_endpoint, second_endpoint, link::AbstractInterlink) -> InterLinkKey

	Register a candidate link of kind `link` between two terminals and return its
	key. Each endpoint is `(spacecraft vector index, terminal)`, where the terminal
	is a terminal number or a terminal mounted on that `SpacecraftModel`, or
	`(-s, antenna)` for ground station `s` (`add_ground_station!`), with an
	antenna number or one of its antennas. The terminals must be on different
	spacecraft, or on a spacecraft and a station (two stations cannot link),
	numbered within each spacecraft's `terminal_count` or the station's antennas,
	and able to carry the link (`terminal_supports`); a link requiring hardware
	at both ends (a `DataLink` with a `ShannonBudget`, or a `LaserBeamLink` with
	mounted recoil or terminal pointing) refuses a virtual terminal. A link registers once. Its key lists the lower endpoint
	first (a station's before a spacecraft's), or, for a directed link
	(`link_directed`), the endpoints in the order given. The link is available
	only when `link_available` allows it, and active only when the scheduler
	selects it.
	"""
	function register_link!(model::InterlinkModel, spacecraft,
		first_endpoint::Tuple{Integer, Any}, second_endpoint::Tuple{Integer, Any}, link::AbstractInterlink)

		# 2.1. Validate the terminals before changing the graph.
		# i) Rule 1: The interlink must join different spacecraft, or a spacecraft and a ground station.
		first_endpoint[1] < 0 && second_endpoint[1] < 0 &&
			throw(ArgumentError("An interlink cannot join two ground stations."))
		first_endpoint[1] != second_endpoint[1] || throw(ArgumentError("An interlink must join different spacecraft."))

		# ii) Rule 2: Both terminals must be valid for their spacecraft or station and able to carry the link.
		sides = link_directed(link) ? (:transmitter, :receiver) : (:either, :either)
		endpoints = map((first_endpoint, second_endpoint), sides) do (satellite, terminal), side
			satellite < 0 && return _station_endpoint(model, satellite, terminal, link, side)
			1 <= satellite <= length(spacecraft) || throw(ArgumentError("Invalid spacecraft vector index $satellite."))
			number = _terminal_number(spacecraft[satellite], terminal)
			1 <= number <= terminal_count(spacecraft[satellite]) ||
				throw(ArgumentError("Invalid terminal $number on spacecraft $satellite."))
			device = _terminal_device(spacecraft[satellite], number)
			device === nothing || terminal_supports(device, link, side) || throw(ArgumentError(
				"Terminal $number of spacecraft $satellite, a $(nameof(typeof(device))), cannot carry a " *
				"$(nameof(typeof(link)))$(side === :either ? "" : " as its $(side)") (terminal_supports)."))
			device === nothing && _link_needs_terminals(link) && throw(ArgumentError(
				"Terminal $number of spacecraft $satellite is virtual, and this $(nameof(typeof(link))) needs terminal " *
				"hardware at both ends$(_needs_terminals_reason(link))."))
			(Int(satellite), number)
		end

		# iii) Rule 3: An undirected link is keyed lower endpoint first, so reversing a pair gives
		# the same key; a directed one keeps its transmitter first.
		key = link_directed(link) ? endpoints : minmax(endpoints...)

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
		push!(graph.edge_properties, link)
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

		# 2.4. Add the edge to its spacecraft's adjacency lists, including empty lists for unlinked spacecraft.
		while length(graph.satellite_edges) < length(spacecraft)
			push!(graph.satellite_edges, Int[])
		end
		for (satellite, _) in key
			satellite > 0 && push!(graph.satellite_edges[satellite], edge)
		end

		return key
	end

	# A ground station end `(node, terminal)` of a link being registered, node
	# `-s` for station `s`, checked as a spacecraft end is.
	function _station_endpoint(model::InterlinkModel, node::Integer, terminal, link, side::Symbol)
		s = -Int(node)
		s <= length(model.ground_stations) || throw(ArgumentError(
			"Invalid ground station $s (node $node): the model has $(length(model.ground_stations)) ground stations " *
			"(add_ground_station!)."))
		station = model.ground_stations[s]
		number = _station_antenna_number(station, terminal)
		1 <= number <= length(station.antennas) ||
			throw(ArgumentError("Invalid antenna $number on ground station $s ($(station.name))."))
		terminal_supports(station.antennas[number], link, side) || throw(ArgumentError(
			"Antenna $number of ground station $s ($(station.name)), a RadioTerminal, cannot carry a " *
			"$(nameof(typeof(link)))$(side === :either ? "" : " as its $(side)") (terminal_supports)."))
		return (Int(node), number)
	end

	"""
	    add_ground_station!(model::InterlinkModel, station::GroundStation) -> Int

	Add `station` to `model`'s ground stations and return its station number `s`,
	its position in registration order. Its antenna `k` is then the link endpoint
	`(-s, k)` for `register_link!`, and `station_pairs` lists a pair with every
	spacecraft for `register_links!`. Each station's name must be its own.
	"""
	function add_ground_station!(model::InterlinkModel, station::GroundStation)
		any(other -> other.name == station.name, model.ground_stations) && throw(ArgumentError(
			"The interlink model already has a ground station named $(repr(station.name))."))
		push!(model.ground_stations, station)
		return length(model.ground_stations)
	end

	station_pairs(spacecraft, model::InterlinkModel) = station_pairs(spacecraft, model.ground_stations)

	# Whether `model` (or `nothing`) has ground stations.
	has_ground_stations(model) = model !== nothing && !isempty(model.ground_stations)

	# Whether a run with interlink model `model` counts each ground station's bits
	# received and sent: with stations and data links, as the spacecraft's counters.
	has_station_counters(model) = has_ground_stations(model) && has_data_links(model)

	# Whether the link end at node `node` is active: a spacecraft's flag, and always for a station.
	@inline _node_active(is_active, node::Int) = node < 0 || is_active[node]

	# The state sample of the link end at node `node`: spacecraft `node` of the run
	# state `state`, with its model when known, or ground station `-node`, with the
	# planet-fixed rotation `planet_frame` at the state's time (`nothing` for the
	# planet's current `L_PI`).
	@inline function _endpoint_sample(model, node::Int, state, spacecraft, context, planet_frame=nothing)
		node > 0 && return link_endpoint(state.sc[node], _spacecraft_at(spacecraft, node), node)
		return _station_sample(model.ground_stations[-node], -node, context.planet, planet_frame)
	end

	# The terminal hardware at the station end `(node, number)`: its antenna.
	_station_antenna(model, (node, number)) = model.ground_stations[-node].antennas[number]

	"""
	    register_allowed_link!(model, spacecraft, first_endpoint, second_endpoint; parameters=InterLinkParameters())

	Deprecated, use `register_link!`: registers a `LaserMomentumLink` given as
	`parameters`. The old name will be removed in SpaceAGORA 1.0.
	"""
	function register_allowed_link!(model::InterlinkModel, spacecraft,
		first_endpoint::vertex, second_endpoint::vertex; parameters::AbstractInterlink=LaserMomentumLink())
		_interlink_deprecation((InterlinkModel, :register_allowed_link!),
			"register_allowed_link! is deprecated and will be removed in SpaceAGORA 1.0; use " *
			"register_link!(model, spacecraft, first_endpoint, second_endpoint, link) instead.")
		return register_link!(model, spacecraft, first_endpoint, second_endpoint, parameters)
	end

	"""
	    all_pairs(spacecraft) -> Vector{Tuple{Int, Int}}

	Every pair `(i, j)` of distinct spacecraft vector indices with `i < j`, in
	order: the candidate pairs of a fully connected interlink graph, for
	`register_links!`.
	"""
	all_pairs(spacecraft) = [(i, j) for i in eachindex(spacecraft) for j in eachindex(spacecraft) if i < j]

	# The spacecraft with a constellation slot, by plane: plane => [(slot, vector index)], by slot.
	# Two spacecraft in one slot have no single neighbour, so they are refused.
	function _slots_by_plane(spacecraft)
		planes = Dict{Int, Vector{Tuple{Int, Int}}}()
		holder = Dict{Tuple{Int, Int}, Int}()
		for (index, sc) in enumerate(spacecraft)
			slot = hasproperty(sc, :constellation_slot) ? sc.constellation_slot : nothing
			slot === nothing && continue
			other = get!(holder, (slot.plane, slot.slot), index)
			other == index || throw(ArgumentError(
				"Spacecraft $other and $index share plane $(slot.plane), slot $(slot.slot); " *
				"pairs by plane and slot need each slot held once (build_constellation gives each its own)."))
			push!(get!(Vector{Tuple{Int, Int}}, planes, slot.plane), (slot.slot, index))
		end
		foreach(sort!, values(planes))
		return planes
	end

	"""
	    same_plane_pairs(spacecraft) -> Vector{Tuple{Int, Int}}

	Every pair `(i, j)`, `i < j`, of spacecraft in the same orbital plane, by their
	`constellation_slot` (set by `build_constellation`), for `register_links!`.
	Spacecraft without a slot take part in no pair; two spacecraft in one slot are
	refused.
	"""
	function same_plane_pairs(spacecraft)
		pairs = Tuple{Int, Int}[]
		for members in values(_slots_by_plane(spacecraft)), a in eachindex(members), b in a+1:length(members)
			push!(pairs, minmax(members[a][2], members[b][2]))
		end
		return sort!(unique!(pairs))
	end

	"""
	    ring_pairs(spacecraft) -> Vector{Tuple{Int, Int}}

	Each spacecraft with the next slot in its plane, and the last slot with the
	first: the in-plane rings of a constellation, by `constellation_slot`. A plane
	of two gives one pair, a plane of one none. Pairs are `(i, j)` with `i < j`.
	"""
	function ring_pairs(spacecraft)
		pairs = Tuple{Int, Int}[]
		for members in values(_slots_by_plane(spacecraft))
			n = length(members)
			n < 2 && continue
			for k in 1:n
				push!(pairs, minmax(members[k][2], members[mod1(k + 1, n)][2]))
			end
		end
		return sort!(unique!(pairs))
	end

	"""
	    grid_pairs(spacecraft; wrap_planes=true) -> Vector{Tuple{Int, Int}}

	The "+grid" of a constellation, by `constellation_slot`: the in-plane rings
	(`ring_pairs`) and each slot with the same slot of the next plane, planes taken
	in plane-number order. With `wrap_planes`, the last plane links to the first;
	pass `false` for a Walker star, whose first and last planes counter-rotate. A
	slot missing from the next plane links nothing across. Pairs are `(i, j)` with
	`i < j`.
	"""
	function grid_pairs(spacecraft; wrap_planes::Bool=true)
		planes = _slots_by_plane(spacecraft)
		pairs = ring_pairs(spacecraft)
		order = sort!(collect(keys(planes)))
		for (k, plane) in enumerate(order)
			k == length(order) && !wrap_planes && break
			next = order[mod1(k + 1, length(order))]
			next == plane && continue
			across = Dict(planes[next])
			for (slot, index) in planes[plane]
				haskey(across, slot) && push!(pairs, minmax(index, across[slot]))
			end
		end
		return sort!(unique!(pairs))
	end

	"""
	    register_links!(model::InterlinkModel, spacecraft, pairs, link; terminals=nothing) -> Vector{InterLinkKey}

	Register one candidate link per spacecraft pair `(i, j)` in `pairs` and return
	their keys in order. The pairs come from `all_pairs(spacecraft)`, from a
	constellation's planes and slots (`same_plane_pairs`, `ring_pairs`,
	`grid_pairs`), from `station_pairs` for ground stations (node `-s` for station
	`s`), or from any vector of index pairs. `link` is a link kind value, or a function
	`(i, j) -> link` that builds each pair's link.

	Each spacecraft's (and station's) links take the terminals that can carry them
	in turn, starting again from the first when they run out; links that share a
	terminal compete for it when the scheduler selects. Pass
	`terminals=(i, j) -> (terminal_i, terminal_j)` to choose them instead. For a
	directed link, `i` is the transmitter.
	"""
	function register_links!(model::InterlinkModel, spacecraft, pairs, link; terminals=nothing)
		turn = zeros(Int, length(spacecraft))
		station_turn = zeros(Int, length(model.ground_stations))
		keys = InterLinkKey[]
		for (i, j) in pairs
			pair_link = link isa AbstractInterlink ? link : link(i, j)
			pair_link isa AbstractInterlink || throw(ArgumentError(
				"register_links!: the link for pair ($i, $j) is a $(typeof(pair_link)), not a link kind (AbstractInterlink)."))
			sides = link_directed(pair_link) ? (:transmitter, :receiver) : (:either, :either)
			terminal_of(node, side) = node < 0 ? _next_antenna!(station_turn, model, node, pair_link, side) :
				_next_terminal!(turn, spacecraft, node, pair_link, side)
			ends = terminals === nothing ? (terminal_of(i, sides[1]), terminal_of(j, sides[2])) : terminals(i, j)
			push!(keys, register_link!(model, spacecraft, (i, ends[1]), (j, ends[2]), pair_link))
		end
		return keys
	end

	# The next terminal of spacecraft `satellite` that can carry `link`, in turn.
	function _next_terminal!(turn::Vector{Int}, spacecraft, satellite::Integer, link, side::Symbol)
		1 <= satellite <= length(spacecraft) || throw(ArgumentError("Invalid spacecraft vector index $satellite."))
		usable = filter(1:terminal_count(spacecraft[satellite])) do number
			device = _terminal_device(spacecraft[satellite], number)
			device === nothing ? !_link_needs_terminals(link) : terminal_supports(device, link, side)
		end
		isempty(usable) && throw(ArgumentError(
			"Spacecraft $satellite has no terminal that can carry a $(nameof(typeof(link)))."))
		turn[satellite] += 1
		return usable[mod1(turn[satellite], length(usable))]
	end

	# The next antenna of the station at node `node` that can carry `link`, in turn.
	function _next_antenna!(turn::Vector{Int}, model::InterlinkModel, node::Integer, link, side::Symbol)
		s = -Int(node)
		s <= length(model.ground_stations) || throw(ArgumentError(
			"Invalid ground station $s (node $node): the model has $(length(model.ground_stations)) ground stations " *
			"(add_ground_station!)."))
		station = model.ground_stations[s]
		usable = filter(number -> terminal_supports(station.antennas[number], link, side), eachindex(station.antennas))
		isempty(usable) && throw(ArgumentError(
			"Ground station $s ($(station.name)) has no antenna that can carry a $(nameof(typeof(link)))."))
		turn[s] += 1
		return usable[mod1(turn[s], length(usable))]
	end

	# The model of spacecraft `index`, when the run's models are known.
	_spacecraft_at(spacecraft, index::Int) = spacecraft === nothing ? nothing : spacecraft[index]

	# The terminal hardware at the two ends of `edge` (`nothing` for a virtual
	# terminal): the run's, once its setup has recorded them, or else looked up on
	# the spacecraft models when they are known.
	function _edge_terminals(model::InterlinkModel, edge::Int, spacecraft)
		edge <= length(model.edge_terminals) && return model.edge_terminals[edge]
		first, second = map(v -> model.graph.vertices[v], model.graph.edges[edge])
		terminal(endpoint) = endpoint[1] < 0 ? _station_antenna(model, endpoint) :
			spacecraft === nothing ? nothing : _terminal_device(spacecraft[endpoint[1]], endpoint[2])
		return (terminal(first), terminal(second))
	end

	# What the link of `edge` does to its two endpoints, in key order (`interlink_wrench`).
	function _edge_wrench(model::InterlinkModel, edge::Int, state; spacecraft=nothing, context=LinkContext(),
		planet_frame=nothing)
		graph = model.graph
		first, second = map(v -> graph.vertices[v], graph.edges[edge])
		a = _endpoint_sample(model, first[1], state, spacecraft, context, planet_frame)
		b = _endpoint_sample(model, second[1], state, spacecraft, context, planet_frame)
		terminal_a, terminal_b = _edge_terminals(model, edge, spacecraft)
		return interlink_wrench(graph.edge_properties[edge], terminal_a, terminal_b, a, b, context)
	end

	# 3. Update link availability
	"""
		update_availability!(model::InterlinkModel, spacecraft, state; is_active=nothing, context=LinkContext(),
		                     terminal_available=nothing)

	Mark each candidate link available when its spacecraft pair is not forbidden,
	both spacecraft are active (when activity flags are supplied), both terminals
	are available (when `terminal_available` gives a flag per graph vertex: in a
	run, terminal hardware the power budget sheds is not), and its link's
	`link_available` allows it (by default, the endpoints are apart and within
	`link_range_m`). A link with a ground station end also needs the spacecraft at
	least the station's `min_elevation_deg` above its horizon, at the planet's
	current frame (`context.planet.L_PI`); a station is always active. Return a
	fresh `Vector{Bool}` in registration order.
	"""
	function update_availability!(model::InterlinkModel, spacecraft, state; is_active=nothing, context=LinkContext(),
		terminal_available=nothing)
		graph = model.graph
		for edge in eachindex(graph.edges)
			first, second = map(v -> graph.vertices[v], graph.edges[edge])

			# i) Rule 1: The spacecraft pair must not be forbidden.
			allowed = !(minmax(first[1], second[1]) in model.forbidden_pairs) &&

				# ii) Rule 2: Both spacecraft must be active when activity flags are supplied.
				(is_active === nothing || (_node_active(is_active, first[1]) && _node_active(is_active, second[1]))) &&

				# iii) Rule 3: Both terminals must be available when terminal flags are supplied.
				(terminal_available === nothing ||
					(terminal_available[graph.edges[edge][1]] && terminal_available[graph.edges[edge][2]]))

			# iv) Rule 4: A ground station must see the spacecraft above its mask.
			# v) Rule 5: The link itself must be able to close (by default, within its range).
			graph.available[edge] = if first[1] > 0 && second[1] > 0
				allowed && link_available(graph.edge_properties[edge],
					_edge_terminals(model, edge, spacecraft)...,
					link_endpoint(state.sc[first[1]], _spacecraft_at(spacecraft, first[1]), first[1]),
					link_endpoint(state.sc[second[1]], _spacecraft_at(spacecraft, second[1]), second[1]), context)
			else
				allowed && _station_link_available(model, edge, first[1], second[1], state, spacecraft, context)
			end
		end
		return copy(graph.available)
	end

	# Rules 4 and 5 for the link of `edge`, between a spacecraft and a ground station.
	function _station_link_available(model::InterlinkModel, edge::Int, first::Int, second::Int, state, spacecraft, context)
		a = _endpoint_sample(model, first, state, spacecraft, context)
		b = _endpoint_sample(model, second, state, spacecraft, context)
		station, craft = first < 0 ? (a, b) : (b, a)
		elevation = _station_elevation_rad(station.spacecraft, context.planet, nothing, craft.pos_ii)
		elevation >= deg2rad(station.spacecraft.min_elevation_deg) || return false
		return link_available(model.graph.edge_properties[edge], _edge_terminals(model, edge, spacecraft)..., a, b, context)
	end

	"""
		available_links(model::InterlinkModel) -> Vector{InterLinkKey}

	The keys of the candidate links marked available at the last
	`update_availability!`, sorted. The run's scheduler records them in each
	`InterLinkHistory` sample, beside the links it selected, when the model
	records them (`records_available`).
	"""
	function available_links(model::InterlinkModel)
		graph = model.graph
		return sort!([map(v -> graph.vertices[v], graph.edges[edge]) for edge in eachindex(graph.edges) if graph.available[edge]])
	end

	# Whether the run's history records the available links: as `record_available`
	# says, or, when it is `nothing`, in a run with a data link.
	records_available(model::InterlinkModel) =
		model.record_available === nothing ? has_data_links(model) : model.record_available

	# 4. The force of the active links on one spacecraft
	"""
		interlink_force_on_spacecraft(model::InterlinkModel, state, satellite; is_active=nothing, spacecraft=nothing, context=LinkContext())

	The sum of the inertial forces (N) the selected links put on spacecraft
	`satellite` (its vector index), each from its link's `interlink_wrench` (by
	default, its `link_force`). Optional `is_active` flags mask links with an
	inactive endpoint. A ground station reads `planet_frame` when supplied,
	otherwise the planet's current `L_PI`; its force is borne by the planet.
	"""
	interlink_force_on_spacecraft(model::InterlinkModel, state, satellite::Int; is_active=nothing,
		spacecraft=nothing, context=LinkContext(), planet_frame=nothing) =
		first(interlink_wrench_on_spacecraft(model, state, satellite; is_active, spacecraft, context, planet_frame))

	# The sums of the inertial forces (N) and of the body torques (N·m) that the
	# selected links put on spacecraft `satellite`, each from its link's
	# `interlink_wrench`, and whether any torque was nonzero.
	function interlink_wrench_on_spacecraft(model::InterlinkModel, state, satellite::Int; is_active=nothing,
		spacecraft=nothing, context=LinkContext(), planet_frame=nothing)
		force = SVector{3, Float64}(0.0, 0.0, 0.0)
		torque = SVector{3, Float64}(0.0, 0.0, 0.0)
		turned = false
		1 <= satellite <= length(model.graph.satellite_edges) || return (force, torque, turned)
		for index in model.graph.satellite_edges[satellite]
			key = map(v -> model.graph.vertices[v], model.graph.edges[index])
			model.graph.active[index] || continue
			is_active === nothing || (_node_active(is_active, key[1][1]) && _node_active(is_active, key[2][1])) || continue
			first_force, first_torque, second_force, second_torque =
				_edge_wrench(model, index, state; spacecraft, context, planet_frame)
			if key[1][1] == satellite
				force += first_force
				iszero(first_torque) || (torque += first_torque; turned = true)
			elseif key[2][1] == satellite
				force += second_force
				iszero(second_torque) || (torque += second_torque; turned = true)
			end
		end
		return (force, torque, turned)
	end

end
