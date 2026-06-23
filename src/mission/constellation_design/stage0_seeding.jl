module Stage0Seeding

using ..ConstellationUtils
using ..ConstellationSlots
using ..CapoIntegration
using ..Stage0FHSG
using Random
using DataFrames, CSV
using Arrow
using Dates
using LinearAlgebra
using Statistics

# Constants
const EARTH_RADIUS_M = 6371.0e3
const EARTH_MU = 3.986004418e14
const ALT_MARGIN_M = 200_000.0
const INC_MARGIN_DEG = 10.0

"""
    normalize_raan_deg(raan_deg::Real) -> Float64

Normalize RAAN to [0, 360) degrees.
"""
function normalize_raan_deg(raan_deg::Real)
    wrapped = mod(float(raan_deg), 360.0)
    return wrapped < 0 ? wrapped + 360.0 : wrapped
end

"""
    normalize_inclination_deg(inc_deg::Real) -> Float64

Normalize inclination to [0, 360) degrees.
"""
function normalize_inclination_deg(inc_deg::Real)
    wrapped = mod(float(inc_deg), 360.0)
    return wrapped < 0 ? wrapped + 360.0 : wrapped
end

"""
    compute_shell_bounds(cluster_df::DataFrame) -> NamedTuple

Derive orbital element bounds for the shell from the cluster's debris data.
Returns a NamedTuple with a_min, a_max, i_min, i_max, raan_min, raan_max.
"""
function compute_shell_bounds(cluster_df::DataFrame)
    a_min = (minimum(cluster_df.a_km) * 1_000.0) - EARTH_RADIUS_M - ALT_MARGIN_M
    a_max = (maximum(cluster_df.a_km) * 1_000.0) - EARTH_RADIUS_M + ALT_MARGIN_M
    inc_deg = normalize_inclination_deg.(Float64.(cluster_df.inc_deg))
    raan_deg = normalize_raan_deg.(Float64.(cluster_df.raan_deg))
    i_min = clamp(minimum(inc_deg) - INC_MARGIN_DEG, 0.0, 360.0)
    i_max = clamp(maximum(inc_deg) + INC_MARGIN_DEG, 0.0, 360.0)
    raan_min = clamp(minimum(raan_deg) - INC_MARGIN_DEG, 0.0, 360.0)
    raan_max = clamp(maximum(raan_deg) + INC_MARGIN_DEG, 0.0, 360.0)
    return (a_min=a_min, a_max=a_max, i_min=i_min, i_max=i_max, raan_min=raan_min, raan_max=raan_max)
end

"""
    compute_cluster_timing(cluster_df::DataFrame; mu) -> NamedTuple

Compute the slowest debris orbital period in the cluster and convert the fixed
48-hour precompute horizon into an integer number of whole orbits.
"""
function compute_cluster_timing(cluster_df::DataFrame; mu::Real=EARTH_MU)
    a_m = cluster_df.a_km .* 1_000.0
    periods_s = 2π .* sqrt.(Float64.(a_m).^3 ./ Float64(mu))
    largest_period_s = maximum(periods_s)
    horizon_hours = 48.0
    horizon_s = horizon_hours * 3600.0
    num_orbits = max(1, floor(Int, horizon_s / largest_period_s))
    return (
        horizon_hours = horizon_hours,
        horizon_s = horizon_s,
        largest_period_s = largest_period_s,
        largest_period_hours = largest_period_s / 3600.0,
        num_orbits = num_orbits,
    )
end

"""
    subset_tag(ids::Vector{Int}) -> String

Canonical string tag for a subset, e.g. [1,3,5] -> "1_3_5".
"""
subset_tag(ids::Vector{Int}) = join(sort(ids), "_")

"""
    generate_stochastic_seeds(config_dict::AbstractDict) -> Dict{String,Any}

Generate stochastic seed configurations for constellation design.
This is a simplified implementation that will be expanded for full parity.
"""
function generate_stochastic_seeds(config_dict::AbstractDict)
    constellation_log("stage0", "Starting stochastic seed generation")
    
    # Extract parameters from config
    const_params = get(config_dict, "constellation_params", Dict{String,Any}())
    sim_params = get(config_dict, "sim_params", Dict{String,Any}())
    
    # Get shell bounds
    a_min = get(const_params, "a_min", 7.0e6)
    a_max = get(const_params, "a_max", 8.0e6)
    i_min = get(const_params, "i_min", 0.0)
    i_max = get(const_params, "i_max", 90.0)
    raan_min = get(const_params, "raan_min", 0.0)
    raan_max = get(const_params, "raan_max", 360.0)
    
    # Generate random seed points within bounds
    rng_seed = Int(get(sim_params, "rng_seed", 67))
    rng = MersenneTwister(rng_seed)
    
    num_seeds = 100  # Default seed count
    seeds = []
    
    for i in 1:num_seeds
        a = rand(rng, Uniform(a_min, a_max))
        inc = rand(rng, Uniform(deg2rad(i_min), deg2rad(i_max)))
        raan = rand(rng, Uniform(deg2rad(raan_min), deg2rad(raan_max)))
        arg_p = rand(rng, Uniform(0.0, 2π))
        ta = rand(rng, Uniform(0.0, 2π))
        e = 0.0  # Circular orbits for debris
        
        push!(seeds, Dict(
            "a" => a,
            "e" => e,
            "inc" => inc,
            "raan" => raan,
            "arg_p" => arg_p,
            "ta" => ta,
        ))
    end
    
    constellation_log("stage0", "Generated $(length(seeds)) stochastic seeds")
    
    return Dict{String,Any}(
        "seeds" => seeds,
        "num_seeds" => length(seeds),
    )
end

"""
    run_stage0_seeding(config_dict::AbstractString) -> Dict{String,Any}

Stage 0 stochastic seeding entry point. Loads a YAML config file and runs
stochastic seed generation for constellation design.

# Arguments
- `config_path::AbstractString`: Path to YAML configuration file

# Returns
- Dictionary containing seed configuration and metadata
"""
function run_stage0_seeding(config_path::AbstractString)
    config_dict = ingest_yaml(config_path)
    return run_stage0_seeding(config_dict)
end

"""
    run_stage0_seeding(config_dict::AbstractDict) -> Dict{String,Any}

Stage 0 seeding entry point. Runs seeding for constellation design using the
provided configuration dictionary. Supports stochastic, RL-based, and FHSG methods.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `method::String`: Seeding method ("stochastic", "rl", "fhsg", or "stochastic_greedy", default: "stochastic")

# Returns
- Dictionary containing seed configuration and metadata
"""
function run_stage0_seeding(config_dict::AbstractDict; method::String="stochastic")
    constellation_log_init!(config_dict; context="stage0_seeding")
    
    try
        if method == "rl"
            constellation_log("stage0", "Starting Stage 0 RL-based seeding")
            seed_result = run_rl_stage0_seeding(config_dict)
        elseif method == "fhsg" || method == "stochastic_greedy"
            constellation_log("stage0", "Starting Stage 0 FHSG (finite horizon stochastic greedy) seeding")
            seed_result = run_fhsg_stage0(config_dict)
        else
            constellation_log("stage0", "Starting Stage 0 stochastic seeding")
            seed_result = generate_stochastic_seeds(config_dict)
        end
        
        # Store results in config_dict for downstream stages
        if haskey(seed_result, "constellation_orbitals")
            config_dict["stage0_constellation"] = seed_result["constellation_orbitals"]
            config_dict["stage0_n_sats"] = seed_result["n_sats"]
        elseif haskey(seed_result, "seeds")
            config_dict["stage0_seeds"] = seed_result["seeds"]
            config_dict["stage0_num_seeds"] = seed_result["num_seeds"]
        elseif haskey(seed_result, "seed_indices")
            config_dict["stage0_seed_indices"] = seed_result["seed_indices"]
            config_dict["stage0_candidate_bank"] = seed_result["candidate_bank"]
            config_dict["stage0_h_fwd_exact_coeffs"] = seed_result["h_fwd_exact_coeffs"]
            config_dict["stage0_h_Wcorr_coeffs"] = seed_result["h_Wcorr_coeffs"]
            config_dict["stage0_support_lift_coeffs"] = seed_result["support_lift_coeffs"]
            config_dict["stage0_backward_lift_coeffs"] = seed_result["backward_lift_coeffs"]
        end
        
        constellation_log("stage0", "Stage 0 completed successfully")
        
        return seed_result
    catch err
        constellation_log_exception("stage0", err)
        rethrow(err)
    finally
        constellation_log_close!()
    end
end

export run_stage0_seeding, compute_shell_bounds, compute_cluster_timing, subset_tag

end # module Stage0Seeding
