module Stage0CacheLoader

using JLD2
using YAML
using Dates
using ..ConstellationConfig

"""
    load_cached_stage0_summary(path::AbstractString) -> Dict{String,Any}

Load cached Stage 0 summary from YAML file.

# Arguments
- `path::AbstractString`: Path to YAML summary file

# Returns
- `Dict{String,Any}`: Stage 0 summary dictionary
"""
function load_cached_stage0_summary(path::AbstractString)
    return YAML.load_file(path)
end

"""
    load_cached_stage0_jld2(path::AbstractString) -> Dict{String,Any}

Load cached Stage 0 data from JLD2 file.

# Arguments
- `path::AbstractString`: Path to JLD2 cache file

# Returns
- `Dict{String,Any}`: Stage 0 cache data including h_fwd, h_neg, support_lift_coeffs
"""
function load_cached_stage0_jld2(path::AbstractString)
    return load(path)
end

"""
    validate_cache_signature(config_dict::AbstractDict, cache_data::AbstractDict) -> Bool

Validate that cached Stage 0 data matches the current configuration.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `cache_data::AbstractDict`: Cached Stage 0 data

# Returns
- `Bool`: True if cache signature matches
"""
function validate_cache_signature(config_dict::AbstractDict, cache_data::AbstractDict)
    # Check if cache has signature metadata
    if !haskey(cache_data, "cache_signature")
        @warn "Cache data does not contain signature metadata"
        return false
    end
    
    cache_sig = cache_data["cache_signature"]
    
    # Compare key parameters
    n_horizons_match = cache_sig["n_horizons"] == config_dict["optimizer_params"]["n_horizons"]
    n_clients_match = cache_sig["n_clients"] == config_dict["client_bounds"]["n_clients"]
    Kd_match = cache_sig["controllable_set_dirs"] == config_dict["optimizer_params"]["controllable_set_dirs"]
    
    return n_horizons_match && n_clients_match && Kd_match
end

"""
    extract_stage0_coefficients(cache_data::AbstractDict) -> Dict{String,Any}

Extract Stage 0 support coefficients from cached data.

# Arguments
- `cache_data::AbstractDict`: Cached Stage 0 data

# Returns
- `Dict{String,Any}`: Dictionary containing h_fwd_exact_coeffs, h_Wcorr_coeffs, support_lift_coeffs
"""
function extract_stage0_coefficients(cache_data::AbstractDict)
    coeffs = Dict{String,Any}()
    
    if haskey(cache_data, "h_fwd_exact_coeffs")
        coeffs["h_fwd_exact_coeffs"] = cache_data["h_fwd_exact_coeffs"]
    end
    
    if haskey(cache_data, "h_Wcorr_coeffs")
        coeffs["h_Wcorr_coeffs"] = cache_data["h_Wcorr_coeffs"]
    end
    
    if haskey(cache_data, "support_lift_coeffs")
        coeffs["support_lift_coeffs"] = cache_data["support_lift_coeffs"]
    end
    
    if haskey(cache_data, "backward_lift_coeffs")
        coeffs["backward_lift_coeffs"] = cache_data["backward_lift_coeffs"]
    end
    
    return coeffs
end

"""
    generate_cache_signature(config_dict::AbstractDict) -> Dict{String,Any}

Generate a cache signature for the current configuration.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- `Dict{String,Any}`: Cache signature dictionary
"""
function generate_cache_signature(config_dict::AbstractDict)
    return Dict{String,Any}(
        "n_horizons" => config_dict["optimizer_params"]["n_horizons"],
        "n_clients" => config_dict["client_bounds"]["n_clients"],
        "controllable_set_dirs" => config_dict["optimizer_params"]["controllable_set_dirs"],
        "max_sats" => config_dict["constellation_params"]["max_sats"],
        "timestamp" => string(now()),
    )
end

export load_cached_stage0_summary, load_cached_stage0_jld2,
       validate_cache_signature, extract_stage0_coefficients, generate_cache_signature

end # module Stage0CacheLoader
