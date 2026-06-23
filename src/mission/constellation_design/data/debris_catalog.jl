module DebrisCatalog

using CSV
using DataFrames
using Random
using ..PhysicsAdapter
using ..ConstellationConfig

"""
    load_debris_catalog(config_dict::AbstractDict) -> DataFrame

Load debris catalog from file or generate synthetic data.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `DataFrame`: Debris catalog with columns [a_km, e, inc_deg, raan_deg, arg_p_deg, ta_deg, norad_id]
"""
function load_debris_catalog(config_dict::AbstractDict)
    debris_params = config_dict["debris_params"]
    data_source = debris_params["data_source"]
    
    if data_source == "file"
        file_path = debris_params["file_path"]
        @assert isfile(file_path) "Debris file not found: $file_path"
        df = CSV.read(file_path, DataFrame)
        return normalize_debris_catalog(df)
    elseif data_source == "sample"
        return load_sample_debris_catalog()
    elseif data_source == "synthetic"
        return generate_synthetic_debris_catalog(config_dict)
    else
        error("Unknown debris data source: $data_source")
    end
end

"""
    normalize_debris_catalog(df::DataFrame) -> DataFrame

Normalize debris catalog to standard column names and units.

# Arguments
- `df::DataFrame`: Input debris catalog

# Returns
- `DataFrame`: Normalized catalog with standard columns
"""
function normalize_debris_catalog(df::DataFrame)
    # Ensure required columns exist
    required_cols = [:a_km, :e, :inc_deg, :raan_deg, :arg_p_deg, :ta_deg]
    
    # Map common column name variations
    col_map = Dict(
        :semimajor_axis => :a_km,
        :a => :a_km,
        :eccentricity => :e,
        :ecc => :e,
        :inclination => :inc_deg,
        :inc => :inc_deg,
        :raan => :raan_deg,
        :omega => :arg_p_deg,
        :argp => :arg_p_deg,
        :true_anomaly => :ta_deg,
        :ta => :ta_deg,
        :norad => :norad_id,
        :id => :norad_id,
    )
    
    for (old_name, new_name) in col_map
        if old_name in names(df) && !(new_name in names(df))
            rename!(df, old_name => new_name)
        end
    end
    
    # Add norad_id if missing
    if !(:norad_id in names(df))
        df[!, :norad_id] = 1:nrow(df)
    end
    
    # Ensure all required columns exist
    for col in required_cols
        if !(col in names(df))
            error("Missing required column: $col")
        end
    end
    
    return select(df, required_cols..., :norad_id)
end

"""
    load_sample_debris_catalog() -> DataFrame

Load sample debris catalog from CAPOConstellation data.

# Returns
- `DataFrame`: Sample debris catalog
"""
function load_sample_debris_catalog()
    # Check if CAPOConstellation sample data is available
    sample_path = joinpath(@__DIR__, "../../../../../../CAPOConstellation/orbital_object_data")
    
    if isdir(sample_path)
        # Try to load from CAPOConstellation
        csv_files = filter(f -> endswith(f, ".csv"), readdir(sample_path))
        if !isempty(csv_files)
            df = CSV.read(joinpath(sample_path, csv_files[1]), DataFrame)
            return normalize_debris_catalog(df)
        end
    end
    
    # Fallback: generate synthetic data
    @warn "Sample debris data not found, generating synthetic data"
    return generate_synthetic_debris_catalog(Dict(
        "debris_params" => Dict(
            "n_clients" => 10,
            "synthetic_seed" => 42,
            "altitude_range_km" => [400, 800],
            "inclination_range_deg" => [0, 90],
        ),
    ))
end

"""
    generate_synthetic_debris_catalog(config_dict::AbstractDict) -> DataFrame

Generate synthetic debris catalog for testing.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `DataFrame`: Synthetic debris catalog
"""
function generate_synthetic_debris_catalog(config_dict::AbstractDict)
    debris_params = config_dict["debris_params"]
    n_clients = debris_params["n_clients"]
    seed = debris_params["synthetic_seed"]
    alt_range = debris_params["altitude_range_km"]
    inc_range = debris_params["inclination_range_deg"]
    
    rng = MersenneTwister(seed)
    
    # Earth radius in km
    Re = 6371.0
    
    # Generate orbital elements
    a_km = Re .+ rand(rng, Uniform(alt_range[1], alt_range[2]), n_clients)
    e = rand(rng, Uniform(0.0, 0.1), n_clients)
    inc_deg = rand(rng, Uniform(inc_range[1], inc_range[2]), n_clients)
    raan_deg = rand(rng, Uniform(0.0, 360.0), n_clients)
    arg_p_deg = rand(rng, Uniform(0.0, 360.0), n_clients)
    ta_deg = rand(rng, Uniform(0.0, 360.0), n_clients)
    norad_id = 1:n_clients
    
    df = DataFrame(
        a_km = a_km,
        e = e,
        inc_deg = inc_deg,
        raan_deg = raan_deg,
        arg_p_deg = arg_p_deg,
        ta_deg = ta_deg,
        norad_id = norad_id,
    )
    
    return df
end

"""
    debris_to_initial_conditions(df::DataFrame, config_dict::AbstractDict) -> Dict{Symbol, Vector{Float64}}

Convert debris catalog to initial conditions format for constellation optimization.

# Arguments
- `df::DataFrame`: Debris catalog
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `Dict{Symbol, Vector{Float64}}`: Initial conditions with keys :a, :e, :inc, :raan, :arg_p, :ta
"""
function debris_to_initial_conditions(df::DataFrame, config_dict::AbstractDict)
    n = nrow(df)
    
    # Convert to SI units (meters, radians)
    a = df.a_km * 1000.0
    e = df.e
    inc = deg2rad.(df.inc_deg)
    raan = deg2rad.(df.raan_deg)
    arg_p = deg2rad.(df.arg_p_deg)
    ta = deg2rad.(df.ta_deg)
    
    return Dict{Symbol, Vector{Float64}}(
        :a => a,
        :e => e,
        :inc => inc,
        :raan => raan,
        :arg_p => arg_p,
        :ta => ta,
    )
end

"""
    propagate_debris_trajectories(df::DataFrame, config_dict::AbstractDict) -> Array{Float64,3}

Propagate debris trajectories for all clients.

# Arguments
- `df::DataFrame`: Debris catalog
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `Array{Float64,3}`: Client trajectories [3, N, P]
"""
function propagate_debris_trajectories(df::DataFrame, config_dict::AbstractDict)
    sim_params = config_dict["sim_params"]
    dt = sim_params["dt"]
    N = sim_params["N"]
    
    debris_ics = debris_to_initial_conditions(df, config_dict)
    P = nrow(df)
    
    planet = get_planet_from_config(config_dict)
    
    # Pre-allocate trajectory array
    trajectories = zeros(Float64, 3, N + 1, P)
    
    for p in 1:P
        a = debris_ics[:a][p]
        e = debris_ics[:e][p]
        inc = debris_ics[:inc][p]
        raan = debris_ics[:raan][p]
        arg_p = debris_ics[:arg_p][p]
        ta = debris_ics[:ta][p]
        
        for k in 0:N
            state = kepler_client_state_adapter(k, dt, a, e, inc, raan, arg_p, ta, planet)
            trajectories[1:3, k + 1, p] = state[1:3]
        end
    end
    
    return trajectories
end

"""
    save_debris_catalog(df::DataFrame, path::AbstractString)

Save debris catalog to CSV file.

# Arguments
- `df::DataFrame`: Debris catalog
- `path::AbstractString`: Output file path
"""
function save_debris_catalog(df::DataFrame, path::AbstractString)
    CSV.write(path, df)
end

export load_debris_catalog, normalize_debris_catalog, load_sample_debris_catalog,
       generate_synthetic_debris_catalog, debris_to_initial_conditions,
       propagate_debris_trajectories, save_debris_catalog

end # module DebrisCatalog
