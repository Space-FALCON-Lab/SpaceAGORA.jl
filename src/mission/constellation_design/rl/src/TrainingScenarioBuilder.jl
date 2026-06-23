module TrainingScenarioBuilder

using YAML
using ..ClusterCombinations
using ..ROEBoundsCalculator

"""
    TrainingScenario

Structure representing a single training scenario.
"""
struct TrainingScenario
    scenario_id::String
    cluster_set::Vector{String}
    cluster_tag::String
    n_clients::Int
    client_orbitals::Matrix{Float64}
    orbital_bounds::Matrix{Float64}
    client_csv_path::String
    effector_params::Dict{String,Any}
    training_config::Dict{String,Any}
    difficulty::String
end

"""
    build_training_scenarios(campaign_config::Dict{String,Any}) -> Vector{TrainingScenario}

Build all training scenarios from campaign configuration.
"""
function build_training_scenarios(campaign_config::Dict{String,Any})
    campaign = get(campaign_config, "campaign", Dict{String,Any}())
    variables = get(campaign_config, "variables", Dict{String,Any}())
    training = get(campaign_config, "training", Dict{String,Any}())
    
    # Get labeled debris CSV path
    labeled_csv = get(campaign, "labeled_debris_csv", "analysis_tools/cluster_output/labeled_debris.csv")
    
    # Get cluster combinations
    cluster_combinations_config = get(campaign, "cluster_combinations", Dict{String,Any}())
    cluster_ids = get(campaign, "cluster_ids", ["1a", "2a", "3a", "4a", "5a"])
    
    if haskey(cluster_combinations_config, "individual") || haskey(cluster_combinations_config, "pairs") || haskey(cluster_combinations_config, "supersets")
        cluster_combinations = generate_specific_combinations(cluster_ids, cluster_combinations_config)
    else
        cluster_combinations = generate_cluster_combinations(cluster_ids)
    end
    
    # Get output directory
    output_dir = get(campaign, "output_dir", "data/rl_training_scenarios")
    mkpath(output_dir)
    
    # Build cluster scenarios
    cluster_scenarios = build_all_scenarios(labeled_csv, cluster_combinations, output_dir)
    
    # Generate parameter combinations
    param_combinations = generate_parameter_combinations(variables)
    
    # Build full training scenarios
    training_scenarios = Vector{TrainingScenario}()
    
    scenario_idx = 1
    for cluster_scenario in cluster_scenarios
        for param_combo in param_combinations
            scenario_id = "scenario_$(lpad(scenario_idx, 4, '0'))"
            cluster_tag = cluster_scenario["cluster_tag"]
            
            # Create training config
            training_config = merge(
                copy(training),
                Dict{String,Any}(
                    "scenario_id" => scenario_id,
                    "cluster_set" => cluster_scenario["cluster_set"],
                    "cluster_tag" => cluster_tag,
                )
            )
            
            # Merge parameter combination
            training_config = merge(training_config, param_combo)
            
            # Compute difficulty
            difficulty = compute_difficulty(cluster_scenario, param_combo)
            
            scenario = TrainingScenario(
                scenario_id,
                cluster_scenario["cluster_set"],
                cluster_tag,
                cluster_scenario["n_clients"],
                cluster_scenario["client_orbitals"],
                cluster_scenario["orbital_bounds"],
                cluster_scenario["client_csv_path"],
                param_combo,
                training_config,
                difficulty,
            )
            
            push!(training_scenarios, scenario)
            scenario_idx += 1
        end
    end
    
    # Sort by difficulty if curriculum learning is enabled
    if Bool(get(training, "curriculum_learning", false))
        difficulty_order = get(training, "difficulty_order", "easy_to_hard")
        if difficulty_order == "easy_to_hard"
            sort!(training_scenarios, by = s -> parse_difficulty(s.difficulty))
        else
            sort!(training_scenarios, by = s -> parse_difficulty(s.difficulty), rev = true)
        end
    end
    
    return training_scenarios
end

"""
    generate_parameter_combinations(variables::Dict{String,Any}) -> Vector{Dict{String,Any}}

Generate all parameter combinations from variable definitions.
Supports full-factorial sampling.
"""
function generate_parameter_combinations(variables::Dict{String,Any})
    mode = get(variables, "mode", "full_factorial")
    
    if mode == "full_factorial"
        return generate_full_factorial(variables)
    else
        @warn "Parameter mode '$mode' not implemented, using full_factorial"
        return generate_full_factorial(variables)
    end
end

"""
    generate_full_factorial(variables::Dict{String,Any}) -> Vector{Dict{String,Any}}

Generate full-factorial parameter combinations.
"""
function generate_full_factorial(variables::Dict{String,Any})
    # Get variable names and their levels
    var_names = String[]
    var_levels = Vector{Vector{Any}}()
    
    for (var_name, var_def) in variables
        if haskey(var_def, "levels")
            push!(var_names, var_name)
            push!(var_levels, var_def["levels"])
        end
    end
    
    if isempty(var_names)
        return [Dict{String,Any}()]
    end
    
    # Generate Cartesian product of all levels
    combinations = Vector{Dict{String,Any}}()
    
    # Recursive Cartesian product
    function cartesian_product(idx::Int, current::Dict{String,Any})
        if idx > length(var_names)
            push!(combinations, copy(current))
            return
        end
        
        var_name = var_names[idx]
        for level in var_levels[idx]
            current[var_name] = level
            cartesian_product(idx + 1, current)
        end
    end
    
    cartesian_product(1, Dict{String,Any}())
    
    return combinations
end

"""
    compute_difficulty(cluster_scenario::Dict, param_combo::Dict) -> String

Compute difficulty level for a scenario based on cluster size and parameters.
"""
function compute_difficulty(cluster_scenario::Dict, param_combo::Dict)
    n_clients = cluster_scenario["n_clients"]
    
    # Difficulty based on number of clients
    if n_clients <= 5
        base_difficulty = "easy"
    elseif n_clients <= 15
        base_difficulty = "medium"
    else
        base_difficulty = "hard"
    end
    
    # Adjust based on laser range (harder with smaller range)
    if haskey(param_combo, "laser_range_m")
        laser_range = param_combo["laser_range_m"]
        if laser_range < 300000
            base_difficulty = increase_difficulty(base_difficulty)
        end
    end
    
    # Adjust based on laser power (harder with lower power)
    if haskey(param_combo, "laser_power_or_thrust_n")
        laser_power = param_combo["laser_power_or_thrust_n"]
        if laser_power < 0.1
            base_difficulty = increase_difficulty(base_difficulty)
        end
    end
    
    return base_difficulty
end

"""
    increase_difficulty(difficulty::String) -> String

Increase difficulty level by one step.
"""
function increase_difficulty(difficulty::String)
    if difficulty == "easy"
        return "medium"
    elseif difficulty == "medium"
        return "hard"
    else
        return "hard"
    end
end

"""
    parse_difficulty(difficulty::String) -> Int

Parse difficulty string to numeric value for sorting.
"""
function parse_difficulty(difficulty::String)
    if difficulty == "easy"
        return 1
    elseif difficulty == "medium"
        return 2
    elseif difficulty == "hard"
        return 3
    else
        return 2
    end
end

"""
    save_training_scenarios(scenarios::Vector{TrainingScenario}, output_path::String)

Save training scenarios to YAML file for later use.
"""
function save_training_scenarios(scenarios::Vector{TrainingScenario}, output_path::String)
    mkpath(dirname(output_path))
    
    scenarios_data = []
    for scenario in scenarios
        scenario_dict = Dict{String,Any}(
            "scenario_id" => scenario.scenario_id,
            "cluster_set" => scenario.cluster_set,
            "cluster_tag" => scenario.cluster_tag,
            "n_clients" => scenario.n_clients,
            "client_csv_path" => scenario.client_csv_path,
            "orbital_bounds" => scenario.orbital_bounds,
            "effector_params" => scenario.effector_params,
            "training_config" => scenario.training_config,
            "difficulty" => scenario.difficulty,
        )
        push!(scenarios_data, scenario_dict)
    end
    
    YAML.write_file(output_path, scenarios_data)
    return output_path
end

"""
    load_training_scenarios(input_path::String) -> Vector{TrainingScenario}

Load training scenarios from YAML file.
"""
function load_training_scenarios(input_path::String)
    scenarios_data = YAML.load_file(input_path)
    
    scenarios = Vector{TrainingScenario}()
    for scenario_dict in scenarios_data
        scenario = TrainingScenario(
            scenario_dict["scenario_id"],
            scenario_dict["cluster_set"],
            scenario_dict["cluster_tag"],
            scenario_dict["n_clients"],
            scenario_dict["client_orbitals"],
            scenario_dict["orbital_bounds"],
            scenario_dict["client_csv_path"],
            scenario_dict["effector_params"],
            scenario_dict["training_config"],
            scenario_dict["difficulty"],
        )
        push!(scenarios, scenario)
    end
    
    return scenarios
end

export TrainingScenario, build_training_scenarios, generate_parameter_combinations
export compute_difficulty, save_training_scenarios, load_training_scenarios

end # module TrainingScenarioBuilder
