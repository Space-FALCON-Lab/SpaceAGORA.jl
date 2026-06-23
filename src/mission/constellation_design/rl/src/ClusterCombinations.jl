module ClusterCombinations

using CSV
using DataFrames
using ..ROEBoundsCalculator

"""
    load_labeled_debris_csv(csv_path::String) -> DataFrame

Load labeled debris CSV file with cluster information.
Expected columns: NORAD_ID, cluster_id, a_km, ecc, inc_deg, raan_deg, argp_deg, true_anomaly_deg
"""
function load_labeled_debris_csv(csv_path::String)
    isfile(csv_path) || error("Labeled debris CSV not found: $csv_path")
    df = CSV.read(csv_path, DataFrame)
    return df
end

"""
    filter_cluster(df::DataFrame, cluster_ids::Vector{String}) -> DataFrame

Filter debris DataFrame to specific cluster IDs.
"""
function filter_cluster(df::DataFrame, cluster_ids::Vector{String})
    ids = Set(cluster_ids)
    return filter(row -> string(row.cluster_id) in ids, df)
end

"""
    generate_cluster_combinations(cluster_ids::Vector{String}) -> Vector{Vector{String}}

Generate all cluster combinations (individual, pairs, supersets).
"""
function generate_cluster_combinations(cluster_ids::Vector{String})
    combinations = Vector{Vector{String}}()
    
    # Individual clusters
    for id in cluster_ids
        push!(combinations, [id])
    end
    
    # Pairs
    for i in 1:length(cluster_ids)
        for j in (i+1):length(cluster_ids)
            push!(combinations, [cluster_ids[i], cluster_ids[j]])
        end
    end
    
    # Triples
    for i in 1:length(cluster_ids)
        for j in (i+1):length(cluster_ids)
            for k in (j+1):length(cluster_ids)
                push!(combinations, [cluster_ids[i], cluster_ids[j], cluster_ids[k]])
            end
        end
    end
    
    # Superset (all clusters)
    if length(cluster_ids) > 3
        push!(combinations, cluster_ids)
    end
    
    return combinations
end

"""
    generate_specific_combinations(cluster_ids::Vector{String}, config::Dict) -> Vector{Vector{String}}

Generate cluster combinations based on configuration.
"""
function generate_specific_combinations(cluster_ids::Vector{String}, config::Dict)
    combinations = Vector{Vector{String}}()
    
    # Individual clusters
    if haskey(config, "individual")
        for id in get(config, "individual", [])
            if id in cluster_ids
                push!(combinations, [id])
            end
        end
    end
    
    # Pairs
    if haskey(config, "pairs")
        for pair in get(config, "pairs", [])
            if all(id in cluster_ids for id in pair)
                push!(combinations, copy(pair))
            end
        end
    end
    
    # Supersets
    if haskey(config, "supersets")
        for superset in get(config, "supersets", [])
            if all(id in cluster_ids for id in superset)
                push!(combinations, copy(superset))
            end
        end
    end
    
    return combinations
end

"""
    convert_to_orbital_elements(df::DataFrame) -> Matrix{Float64}

Convert debris DataFrame to orbital elements matrix.
Converts from [a_km, ecc, inc_deg, raan_deg, argp_deg, true_anomaly_deg] to [a, e, inc, raan, arg_p, ta]
"""
function convert_to_orbital_elements(df::DataFrame)
    deg2rad = π / 180.0
    
    n_rows = nrow(df)
    orbitals = zeros(n_rows, 6)
    
    for i in 1:n_rows
        orbitals[i, 1] = df.a_km[i] * 1000.0  # Convert km to m
        orbitals[i, 2] = df.ecc[i]
        orbitals[i, 3] = df.inc_deg[i] * deg2rad
        orbitals[i, 4] = df.raan_deg[i] * deg2rad
        orbitals[i, 5] = df.argp_deg[i] * deg2rad
        orbitals[i, 6] = df.true_anomaly_deg[i] * deg2rad
    end
    
    return orbitals
end

"""
    write_client_csv(orbitals::Matrix{Float64}, output_path::String)

Write orbital elements to CSV file for client input.
"""
function write_client_csv(orbitals::Matrix{Float64}, output_path::String)
    mkpath(dirname(output_path))
    
    df = DataFrame(
        a = orbitals[:, 1],
        e = orbitals[:, 2],
        inc = orbitals[:, 3],
        raan = orbitals[:, 4],
        arg_p = orbitals[:, 5],
        ta = orbitals[:, 6],
    )
    
    CSV.write(output_path, df)
    return output_path
end

"""
    build_cluster_scenario(df::DataFrame, cluster_set::Vector{String}, output_dir::String) -> Dict{String,Any}

Build a training scenario for a specific cluster combination.
"""
function build_cluster_scenario(df::DataFrame, cluster_set::Vector{String}, output_dir::String)
    # Filter to cluster set
    cluster_df = filter_cluster(df, cluster_set)
    
    if nrow(cluster_df) == 0
        @warn "No debris found for cluster set: $cluster_set"
        return nothing
    end
    
    # Convert to orbital elements
    orbitals = convert_to_orbital_elements(cluster_df)
    
    # Compute orbital bounds
    bounds = compute_orbital_bounds_from_cluster(orbitals)
    
    # Write client CSV
    cluster_tag = join(cluster_set, "_")
    client_csv_path = joinpath(output_dir, "clients_$(cluster_tag).csv")
    write_client_csv(orbitals, client_csv_path)
    
    return Dict{String,Any}(
        "cluster_set" => cluster_set,
        "n_clients" => size(orbitals, 1),
        "client_orbitals" => orbitals,
        "orbital_bounds" => bounds,
        "client_csv_path" => client_csv_path,
    )
end

"""
    build_all_scenarios(labeled_csv_path::String, cluster_combinations::Vector{Vector{String}}, output_dir::String) -> Vector{Dict{String,Any}}

Build training scenarios for all cluster combinations.
"""
function build_all_scenarios(labeled_csv_path::String, cluster_combinations::Vector{Vector{String}}, output_dir::String)
    # Load labeled debris
    df = load_labeled_debris_csv(labeled_csv_path)
    
    # Get all unique cluster IDs
    all_cluster_ids = unique(string.(df.cluster_id))
    
    scenarios = Vector{Dict{String,Any}}()
    
    for cluster_set in cluster_combinations
        scenario = build_cluster_scenario(df, cluster_set, output_dir)
        if scenario !== nothing
            push!(scenarios, scenario)
        end
    end
    
    return scenarios
end

export load_labeled_debris_csv, filter_cluster, generate_cluster_combinations, generate_specific_combinations
export convert_to_orbital_elements, write_client_csv, build_cluster_scenario, build_all_scenarios

end # module ClusterCombinations
