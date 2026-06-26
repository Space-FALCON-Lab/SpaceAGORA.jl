module Scenarios

using CSV
using DataFrames
using YAML
using ..ConstellationUtils

"""
    load_training_scenarios() -> Vector{Dict{String,Any}}

Load all debris cluster configurations from data/debris_clusters/.
Returns a vector of config dictionaries for training.
"""
function load_training_scenarios()
    cluster_dir = joinpath(@__DIR__, "../../../../../../data/debris_clusters")
    
    if !isdir(cluster_dir)
        @warn "Debris cluster directory not found: $cluster_dir"
        return Dict{String,Any}[]
    end
    
    scenarios = []
    cluster_files = filter(f -> endswith(f, ".csv") && startswith(basename(f), "clients_cluster_"), 
                          readdir(cluster_dir))
    
    for file in cluster_files
        config = load_cluster_from_csv(joinpath(cluster_dir, file))
        push!(scenarios, config)
    end
    
    return scenarios
end

"""
    load_cluster_from_csv(csv_path::String) -> Dict{String,Any}

Load a single debris cluster from CSV file and return config dictionary.
"""
function load_cluster_from_csv(csv_path::String)
    df = CSV.read(csv_path, DataFrame)
    
    # Extract cluster ID from filename
    filename = basename(csv_path)
    cluster_id = replace(filename, "clients_cluster_" => "", ".csv" => "")
    
    # Parse orbital elements
    client_orbitals = Matrix{Float64}(df[:, [:a, :e, :inc, :raan, :arg_p, :ta]])
    
    # Compute orbital bounds
    a_min = minimum(df.a)
    a_max = maximum(df.a)
    e_min = minimum(df.e)
    e_max = maximum(df.e)
    inc_min = minimum(df.inc)
    inc_max = maximum(df.inc)
    raan_min = minimum(df.raan)
    raan_max = maximum(df.raan)
    arg_p_min = minimum(df.arg_p)
    arg_p_max = maximum(df.arg_p)
    ta_min = minimum(df.ta)
    ta_max = maximum(df.ta)
    
    orbital_bounds = hcat(
        [a_min, e_min, inc_min, raan_min, arg_p_min, ta_min],
        [a_max, e_max, inc_max, raan_max, arg_p_max, ta_max]
    )
    
    return Dict{String,Any}(
        "cluster_id" => cluster_id,
        "n_clients" => size(client_orbitals, 1),
        "client_orbitals" => client_orbitals,
        "orbital_bounds" => orbital_bounds,
        "csv_path" => csv_path,
    )
end

"""
    load_superset_csv(cluster_ids::Vector{Int}) -> Dict{String,Any}

Load and merge multiple debris clusters into a superset.
"""
function load_superset_csv(cluster_ids::Vector{Int})
    cluster_dir = joinpath(@__DIR__, "../../../../../../data/debris_clusters")
    
    all_orbitals = []
    for id in cluster_ids
        file = joinpath(cluster_dir, "clients_cluster_$id.csv")
        if isfile(file)
            df = CSV.read(file, DataFrame)
            orbitals = Matrix{Float64}(df[:, [:a, :e, :inc, :raan, :arg_p, :ta]])
            append!(all_orbitals, eachrow(orbitals))
        end
    end
    
    if isempty(all_orbitals)
        error("No valid cluster files found for IDs: $cluster_ids")
    end
    
    client_orbitals = vcat(all_orbitals...)
    
    # Compute combined bounds
    orbital_bounds = hcat(
        minimum(client_orbitals, dims=1)',
        maximum(client_orbitals, dims=1)'
    )
    
    return Dict{String,Any}(
        "cluster_id" => join(cluster_ids, "_"),
        "n_clients" => size(client_orbitals, 1),
        "client_orbitals" => client_orbitals,
        "orbital_bounds" => orbital_bounds,
    )
end

export load_training_scenarios, load_cluster_from_csv, load_superset_csv

end # module Scenarios
