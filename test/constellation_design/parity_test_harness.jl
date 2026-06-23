module ParityTestHarness

using JLD2
using YAML
using LinearAlgebra
using ..ConstellationConfig
using ..Stage0CacheLoader
using ..Stage1Controllable
using ..Stage2OCPVerification

"""
    ParityTestResult

Result of a parity test between CAPOConstellation and SpaceAGORA.
"""
Base.@kwdef struct ParityTestResult
    test_name::String
    passed::Bool
    beta_match::Bool
    z_tube_match::Bool
    stage2_match::Bool
    beta_error::Float64
    z_tube_error::Float64
    stage2_error::Float64
    details::Dict{String,Any}
end

"""
    load_capo_baseline(cluster_id::Int, n_horizons::Int=5; capo_root::AbstractString="/home/robbie/capo/CAPOConstellation",
                       cluster_dir::AbstractString="fresh_stage0_lp_warm_3horizon_clusters_1_2_20260602_175956") -> Dict{String,Any}

Load CAPOConstellation baseline results for a given cluster.

# Arguments
- `cluster_id::Int`: Cluster ID
- `n_horizons::Int=5`: Number of horizons
- `capo_root::AbstractString="/home/robbie/capo/CAPOConstellation"`: CAPOConstellation repository root
- `cluster_dir::AbstractString`: Cluster output directory name

# Returns
- `Dict{String,Any}`: Baseline results dictionary

# Note
This function loads Stage 0 coefficients. For full Stage 1/Stage 2 parity testing,
you need to run the full CAPOConstellation pipeline and save the results to a
non-gitignored location, then load them with `load_capo_stage1_baseline` and
`load_capo_stage2_baseline`.
"""
function load_capo_baseline(cluster_id::Int, n_horizons::Int=5; capo_root::AbstractString="/home/robbie/capo/CAPOConstellation",
                           cluster_dir::AbstractString="fresh_stage0_lp_warm_3horizon_clusters_1_2_20260602_175956")
    # Path to CAPOConstellation cluster output
    stage0_path = joinpath(capo_root, "analysis_tools/cluster_output", cluster_dir, 
                          "cluster_$(cluster_id)/results/stage0_summary.jld2")
    
    if !isfile(stage0_path)
        error("Stage 0 baseline not found at $stage0_path")
    end
    
    stage0_data = load(stage0_path)
    
    # Extract coefficients for Stage 1
    coeffs = extract_stage0_coefficients(stage0_data)
    
    # Load cluster configuration
    config_path = joinpath(capo_root, "analysis_tools/cluster_output", cluster_dir,
                          "cluster_$(cluster_id)/inputs/stage0_config.yaml")
    
    config_dict = nothing
    if isfile(config_path)
        config_dict = YAML.load_file(config_path)
    end
    
    return Dict{String,Any}(
        "stage0" => stage0_data,
        "coefficients" => coeffs,
        "config" => config_dict,
        "cluster_id" => cluster_id,
        "n_horizons" => n_horizons,
        "stage0_path" => stage0_path,
    )
end

"""
    load_capo_stage1_baseline(cluster_id::Int, baseline_path::AbstractString) -> Dict{String,Any}

Load CAPOConstellation Stage 1 baseline results from a JLD2 file.

# Arguments
- `cluster_id::Int`: Cluster ID
- `baseline_path::AbstractString`: Path to Stage 1 baseline JLD2 file

# Returns
- `Dict{String,Any}`: Stage 1 baseline results

# Note
This requires running the full CAPOConstellation Stage 1 pipeline and saving
the results (beta, z_tube, etc.) to a JLD2 file outside the gitignored
cluster_output directory.
"""
function load_capo_stage1_baseline(cluster_id::Int, baseline_path::AbstractString)
    if !isfile(baseline_path)
        error("Stage 1 baseline not found at $baseline_path")
    end
    
    baseline_data = load(baseline_path)
    
    return Dict{String,Any}(
        "beta" => get(baseline_data, "beta", nothing),
        "z_tube" => get(baseline_data, "z_tube", nothing),
        "num_active" => get(baseline_data, "num_active", nothing),
        "objective_value" => get(baseline_data, "objective_value", nothing),
        "cluster_id" => cluster_id,
        "baseline_path" => baseline_path,
    )
end

"""
    load_capo_stage2_baseline(cluster_id::Int, baseline_path::AbstractString) -> Dict{String,Any}

Load CAPOConstellation Stage 2 baseline results from a JLD2 file.

# Arguments
- `cluster_id::Int`: Cluster ID
- `baseline_path::AbstractString`: Path to Stage 2 baseline JLD2 file

# Returns
- `Dict{String,Any}`: Stage 2 baseline results

# Note
This requires running the full CAPOConstellation Stage 2 pipeline and saving
the results (verified status, precheck slack, segments) to a JLD2 file outside
the gitignored cluster_output directory.
"""
function load_capo_stage2_baseline(cluster_id::Int, baseline_path::AbstractString)
    if !isfile(baseline_path)
        error("Stage 2 baseline not found at $baseline_path")
    end
    
    baseline_data = load(baseline_path)
    
    return Dict{String,Any}(
        "verified" => get(baseline_data, "verified", nothing),
        "all_horizons_solved" => get(baseline_data, "all_horizons_solved", nothing),
        "precheck_slack" => get(baseline_data, "precheck_slack", nothing),
        "segments" => get(baseline_data, "segments", nothing),
        "cluster_id" => cluster_id,
        "baseline_path" => baseline_path,
    )
end

"""
    run_spaceagora_stage1(config_dict::AbstractDict, capo_baseline::AbstractDict) -> Dict{String,Any}

Run SpaceAGORA Stage 1 with cached Stage 0 coefficients.

# Arguments
- `config_dict::AbstractDict`: SpaceAGORA configuration
- `capo_baseline::AbstractDict`: CAPOConstellation baseline with coefficients

# Returns
- `Dict{String,Any}`: SpaceAGORA Stage 1 results
"""
function run_spaceagora_stage1(config_dict::AbstractDict, capo_baseline::AbstractDict)
    # Inject cached coefficients
    config_dict["optimizer_params"]["use_cached_stage0"] = true
    config_dict["optimizer_params"]["cached_stage0_data"] = capo_baseline["coefficients"]
    
    # Run Stage 1
    result = run_stage1_controllable_optimization(config_dict)
    
    return result
end

"""
    compare_arrays(a::AbstractArray, b::AbstractArray; rtol::Float64=1e-4, atol::Float64=1e-8) -> (Bool, Float64)

Compare two arrays with relative and absolute tolerance.

# Arguments
- `a::AbstractArray`: First array
- `b::AbstractArray`: Second array
- `rtol::Float64=1e-4`: Relative tolerance
- `atol::Float64=1e-8`: Absolute tolerance

# Returns
- `(Bool, Float64)`: (match status, max relative error)
"""
function compare_arrays(a::AbstractArray, b::AbstractArray; rtol::Float64=1e-4, atol::Float64=1e-8)
    if size(a) != size(b)
        return (false, Inf)
    end
    
    max_error = 0.0
    for i in eachindex(a)
        diff = abs(a[i] - b[i])
        scale = max(abs(a[i]), abs(b[i]), atol)
        rel_error = diff / scale
        max_error = max(max_error, rel_error)
    end
    
    match = max_error <= rtol
    return (match, max_error)
end

"""
    run_parity_test_stage1(cluster_id::Int, config_path::AbstractString; 
                           stage1_baseline_path::Union{AbstractString, Nothing}=nothing,
                           rtol::Float64=1e-4) -> ParityTestResult

Run Stage 1 parity test for a given cluster.

# Arguments
- `cluster_id::Int`: Cluster ID
- `config_path::AbstractString`: SpaceAGORA configuration path
- `stage1_baseline_path::Union{AbstractString, Nothing}=nothing`: Path to CAPO Stage 1 baseline JLD2
- `rtol::Float64=1e-4`: Relative tolerance

# Returns
- `ParityTestResult`: Test result
"""
function run_parity_test_stage1(cluster_id::Int, config_path::AbstractString; 
                               stage1_baseline_path::Union{AbstractString, Nothing}=nothing,
                               rtol::Float64=1e-4)
    # Load CAPOConstellation Stage 0 baseline
    capo_baseline = load_capo_baseline(cluster_id)
    
    # Load SpaceAGORA configuration
    config_dict = parse_constellation_config(config_path)
    
    # Run SpaceAGORA Stage 1
    spaceagora_result = run_spaceagora_stage1(config_dict, capo_baseline)
    
    # Load CAPO Stage 1 baseline if available
    if stage1_baseline_path !== nothing && isfile(stage1_baseline_path)
        capo_stage1 = load_capo_stage1_baseline(cluster_id, stage1_baseline_path)
        
        # Compare beta
        beta_match, beta_error = compare_arrays(spaceagora_result["beta"], capo_stage1["beta"]; rtol=rtol)
        
        # Compare z_tube
        z_tube_match, z_tube_error = compare_arrays(spaceagora_result["z_tube"], capo_stage1["z_tube"]; rtol=rtol)
        
        # Compare num_active
        num_active_match = spaceagora_result["num_active"] == capo_stage1["num_active"]
        
        passed = beta_match && z_tube_match && num_active_match
        
        return ParityTestResult(
            test_name="Stage 1 Cluster $cluster_id",
            passed=passed,
            beta_match=beta_match,
            z_tube_match=z_tube_match,
            stage2_match=true,
            beta_error=beta_error,
            z_tube_error=z_tube_error,
            stage2_error=0.0,
            details=Dict{String,Any}(
                "cluster_id" => cluster_id,
                "spaceagora_num_active" => spaceagora_result["num_active"],
                "capo_num_active" => capo_stage1["num_active"],
                "spaceagora_feasible" => spaceagora_result["feasible"],
                "num_active_match" => num_active_match,
                "baseline_path" => stage1_baseline_path,
            ),
        )
    else
        # No baseline available - just check that SpaceAGORA runs successfully
        return ParityTestResult(
            test_name="Stage 1 Cluster $cluster_id (no baseline)",
            passed=true,
            beta_match=true,
            z_tube_match=true,
            stage2_match=true,
            beta_error=0.0,
            z_tube_error=0.0,
            stage2_error=0.0,
            details=Dict{String,Any}(
                "cluster_id" => cluster_id,
                "spaceagora_num_active" => spaceagora_result["num_active"],
                "spaceagora_feasible" => spaceagora_result["feasible"],
                "note" => "No CAPO Stage 1 baseline available - test only checks SpaceAGORA runs successfully",
            ),
        )
    end
end

"""
    run_parity_test_stage2(cluster_id::Int, config_path::AbstractString; 
                           stage2_baseline_path::Union{AbstractString, Nothing}=nothing,
                           rtol::Float64=1e-3) -> ParityTestResult

Run Stage 2 parity test for a given cluster.

# Arguments
- `cluster_id::Int`: Cluster ID
- `config_path::AbstractString`: SpaceAGORA configuration path
- `stage2_baseline_path::Union{AbstractString, Nothing}=nothing`: Path to CAPO Stage 2 baseline JLD2
- `rtol::Float64=1e-3`: Relative tolerance

# Returns
- `ParityTestResult`: Test result
"""
function run_parity_test_stage2(cluster_id::Int, config_path::AbstractString; 
                               stage2_baseline_path::Union{AbstractString, Nothing}=nothing,
                               rtol::Float64=1e-3)
    # Load CAPOConstellation Stage 0 baseline
    capo_baseline = load_capo_baseline(cluster_id)
    
    # Load SpaceAGORA configuration
    config_dict = parse_constellation_config(config_path)
    
    # Run SpaceAGORA Stage 1 first
    config_dict["optimizer_params"]["use_cached_stage0"] = true
    config_dict["optimizer_params"]["cached_stage0_data"] = capo_baseline["coefficients"]
    stage1_result = run_stage1_controllable_optimization(config_dict)
    
    # Run SpaceAGORA Stage 2
    config_dict["stage1_results"] = stage1_result
    stage2_result = run_stage2_ocp_verification(config_dict)
    
    # Load CAPO Stage 2 baseline if available
    if stage2_baseline_path !== nothing && isfile(stage2_baseline_path)
        capo_stage2 = load_capo_stage2_baseline(cluster_id, stage2_baseline_path)
        
        # Compare verified status
        verified_match = stage2_result["verified"] == capo_stage2["verified"]
        
        # Compare all_horizons_solved
        all_solved_match = stage2_result["all_horizons_solved"] == capo_stage2["all_horizons_solved"]
        
        # Compare precheck_slack
        precheck_slack_match, precheck_slack_error = compare_arrays(
            [stage2_result["precheck_slack"]], 
            [capo_stage2["precheck_slack"]]; 
            rtol=rtol
        )
        
        # Compare segment counts
        segment_count_match = length(stage2_result["segments"]) == length(capo_stage2["segments"])
        
        passed = verified_match && all_solved_match && precheck_slack_match && segment_count_match
        
        return ParityTestResult(
            test_name="Stage 2 Cluster $cluster_id",
            passed=passed,
            beta_match=true,
            z_tube_match=true,
            stage2_match=passed,
            beta_error=0.0,
            z_tube_error=0.0,
            stage2_error=precheck_slack_error,
            details=Dict{String,Any}(
                "cluster_id" => cluster_id,
                "spaceagora_verified" => stage2_result["verified"],
                "capo_verified" => capo_stage2["verified"],
                "spaceagora_all_horizons_solved" => stage2_result["all_horizons_solved"],
                "capo_all_horizons_solved" => capo_stage2["all_horizons_solved"],
                "verified_match" => verified_match,
                "all_solved_match" => all_solved_match,
                "segment_count_match" => segment_count_match,
                "baseline_path" => stage2_baseline_path,
            ),
        )
    else
        # No baseline available - just check that SpaceAGORA runs successfully
        return ParityTestResult(
            test_name="Stage 2 Cluster $cluster_id (no baseline)",
            passed=true,
            beta_match=true,
            z_tube_match=true,
            stage2_match=true,
            beta_error=0.0,
            z_tube_error=0.0,
            stage2_error=0.0,
            details=Dict{String,Any}(
                "cluster_id" => cluster_id,
                "spaceagora_verified" => stage2_result["verified"],
                "spaceagora_all_horizons_solved" => stage2_result["all_horizons_solved"],
                "note" => "No CAPO Stage 2 baseline available - test only checks SpaceAGORA runs successfully",
            ),
        )
    end
end

export ParityTestResult, load_capo_baseline, load_capo_stage1_baseline, load_capo_stage2_baseline,
       run_spaceagora_stage1, compare_arrays, run_parity_test_stage1, run_parity_test_stage2

end # module ParityTestHarness
