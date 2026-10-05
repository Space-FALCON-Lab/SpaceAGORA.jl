module InterLinkModels

using LinearAlgebra

export InterLinkModel, InterLinkParameters, InterLinkState, InterLinkConnection
export TerminalEndpoint, InterLinkKey, InterLinkScheduleSample
export register_candidate!, update_availability!

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

# 1.3. Interlink Model Struct
mutable struct InterLinkModel
	link_type::Symbol
	linkgraph::Dict{InterLinkKey, InterLinkConnection}
	battery_energy_threshold::Float64
	tempurature_threshold::Float64
	forbidden_pairs::Set{Tuple{Int, Int}}
	history::Vector{InterLinkScheduleSample}

	function InterLinkModel(; link_type::Symbol=:laser, battery_energy_threshold::Real=50.0,
		tempurature_threshold::Real=0.0, forbidden_pairs=Tuple{Int, Int}[])
		link_type in (:laser, :radio) || throw(ArgumentError("link_type must be :laser or :radio."))
		0 <= battery_energy_threshold <= 100 || throw(ArgumentError("battery_energy_threshold must be in [0, 100]."))
		0 <= tempurature_threshold <= 100 || throw(ArgumentError("tempurature_threshold must be in [0, 100]."))
		pairs = Set{Tuple{Int, Int}}(minmax(pair...) for pair in forbidden_pairs)
		new(link_type, Dict{InterLinkKey, InterLinkConnection}(), Float64(battery_energy_threshold),
			Float64(tempurature_threshold), pairs, InterLinkScheduleSample[])
	end
end

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
	haskey(model.linkgraph, key) && throw(ArgumentError("Candidate $key is already registered."))
	model.linkgraph[key] = InterLinkConnection(parameters, InterLinkState())
	return key
end

# 3. Update availability for one candidate pair
function update_availability!(model::InterLinkModel, key::InterLinkKey, spacecraft, state; is_active=nothing)
	connection = model.linkgraph[key]
	first_satellite, second_satellite = key[1][1], key[2][1]

	# i) rule 1: The interlink must be within range
	separation = norm(state.sc[first_satellite].pos - state.sc[second_satellite].pos)
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
