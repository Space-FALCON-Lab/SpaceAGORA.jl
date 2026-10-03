module InterLinkModels

using LinearAlgebra
using StaticArrays

export InterLinkModel, InterLinkParameters, InterLinkState, InterLinkConnection
export TerminalEndpoint, InterLinkKey, InterLinkScheduleSample
export register_candidate!, update_availability!, force_on_endpoint

const TerminalEndpoint = Tuple{Int, Int}
const InterLinkKey = Tuple{TerminalEndpoint, TerminalEndpoint}

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

struct InterLinkScheduleSample
	time::Float64
	active::Vector{InterLinkKey}
end

struct InterLinkMatchingWorkspace
	endpoints::Matrix{Int}
	neighbors::Vector{Vector{Int}}
	parent::Vector{Int}
	order::Vector{Int}
	stack::Vector{Int}
	free_score::Vector{Float64}
	blocked_score::Vector{Float64}
	choice::Vector{Int}
	blocked::Vector{Bool}
	weights::Vector{Float64}
	local_index::Vector{Int}
end

mutable struct InterLinkModel
	link_type::Symbol
	linkgraph::Dict{InterLinkKey, InterLinkConnection}
	battery_energy_threshold::Float64
	tempurature_threshold::Float64
	forbidden_pairs::Set{Tuple{Int, Int}}
	eligibility::Union{Nothing, Function}
	active_link_penalty::Float64
	history::Vector{InterLinkScheduleSample}
	edge_keys::Vector{InterLinkKey}
	connections::Vector{InterLinkConnection}
	edge_index::Dict{InterLinkKey, Int}
	satellite_edges::Vector{Vector{Int}}
	available::Vector{Bool}
	availability_transition::Vector{Int8}
	matching::Union{Nothing, InterLinkMatchingWorkspace}

	function InterLinkModel(; link_type::Symbol=:laser, battery_energy_threshold::Real=50.0,
		tempurature_threshold::Real=0.0, forbidden_pairs=Tuple{Int, Int}[],
		eligibility::Union{Nothing, Function}=nothing, active_link_penalty::Real=0.0)
		link_type in (:laser, :radio) || throw(ArgumentError("link_type must be :laser or :radio."))
		0 <= battery_energy_threshold <= 100 || throw(ArgumentError("battery_energy_threshold must be in [0, 100]."))
		0 <= tempurature_threshold <= 100 || throw(ArgumentError("tempurature_threshold must be in [0, 100]."))
		isfinite(active_link_penalty) && active_link_penalty >= 0 ||
			throw(ArgumentError("active_link_penalty must be finite and nonnegative."))
		pairs = Set{Tuple{Int, Int}}(minmax(pair...) for pair in forbidden_pairs)
		new(link_type, Dict{InterLinkKey, InterLinkConnection}(), Float64(battery_energy_threshold),
			Float64(tempurature_threshold), pairs, eligibility, Float64(active_link_penalty),
			InterLinkScheduleSample[], InterLinkKey[], InterLinkConnection[],
			Dict{InterLinkKey, Int}(), Vector{Int}[], Bool[], Int8[], nothing)
	end
end

function force_on_endpoint end

function register_candidate!(model::InterLinkModel, spacecraft,
	first_endpoint::TerminalEndpoint, second_endpoint::TerminalEndpoint;
	parameters::InterLinkParameters=InterLinkParameters())
	first_endpoint[1] != second_endpoint[1] || throw(ArgumentError("An interlink must join different spacecraft."))
	for (satellite, terminal) in (first_endpoint, second_endpoint)
		1 <= satellite <= length(spacecraft) || throw(ArgumentError("Invalid spacecraft vector index $satellite."))
		1 <= terminal <= spacecraft[satellite].n_terminal || throw(ArgumentError("Invalid terminal $terminal on spacecraft $satellite."))
	end
	key = minmax(first_endpoint, second_endpoint)
	haskey(model.linkgraph, key) && throw(ArgumentError("Candidate $key is already registered."))
	connection = InterLinkConnection(parameters, InterLinkState())
	model.linkgraph[key] = connection
	push!(model.edge_keys, key)
	push!(model.connections, connection)
	index = length(model.edge_keys)
	model.edge_index[key] = index
	while length(model.satellite_edges) < length(spacecraft)
		push!(model.satellite_edges, Int[])
	end
	for (satellite, _) in key
		push!(model.satellite_edges[satellite], index)
	end
	push!(model.available, false)
	push!(model.availability_transition, 0)
	model.matching = nothing
	return key
end

function _update_availability!(model::InterLinkModel, index::Int, spacecraft, state, is_active)
	key = model.edge_keys[index]
	connection = model.connections[index]
	first_satellite, second_satellite = key[1][1], key[2][1]
	separation = norm(SVector{3, Float64}(state.sc[first_satellite].pos) -
		SVector{3, Float64}(state.sc[second_satellite].pos))
	available = 0 < separation <= connection.parameters.range &&
		!(minmax(first_satellite, second_satellite) in model.forbidden_pairs)
	for (satellite, terminal) in key
		vehicle = spacecraft[satellite]
		available &= 1 <= terminal <= vehicle.n_terminal &&
			vehicle.battery_energy_index >= model.battery_energy_threshold &&
			vehicle.tempurature_index >= model.tempurature_threshold &&
			(is_active === nothing || is_active[satellite])
	end
	if available && model.eligibility !== nothing
		available = model.eligibility(key, spacecraft, state)
	end
	model.availability_transition[index] = Int8(available) - Int8(connection.state.available)
	model.available[index] = available
	connection.state.available = available
	return available
end

function update_availability!(model::InterLinkModel, key::InterLinkKey, spacecraft, state; is_active=nothing)
	return _update_availability!(model, model.edge_index[key], spacecraft, state, is_active)
end

function update_availability!(model::InterLinkModel, spacecraft, state; is_active=nothing)
	for index in eachindex(model.connections)
		_update_availability!(model, index, spacecraft, state, is_active)
	end
	return model.available
end

end
