using Serialization
using DataFrames
using Arrow
using Dates
using StaticArrays

#usage
# # Run the full test
# julia graph_test.jl

# # Just inspect data structure
# TEST_MODE=2 julia graph_test.jl

# Include necessary modules
include("../config.jl")
include("../utils/plot_data_multi.jl")
include("../utils/Define_mission.jl")
import .config

function load_latest_simulation_data(results_dir::String = joinpath(pwd(), "results"))
    """
    Load the most recent aerobraking_outputs.jls file from the results directory
    """
    if !isdir(results_dir)
        error("Results directory does not exist: $results_dir")
    end
    
    # Find all .jls files matching the pattern
    jls_files = filter(f -> startswith(f, "aerobraking_outputs_") && endswith(f, ".jls"), readdir(results_dir))
    
    if isempty(jls_files)
        error("No aerobraking_outputs files found in $results_dir")
    end
    
    # Sort by filename (which includes timestamp) to get the latest
    sort!(jls_files, rev=true)
    latest_file = joinpath(results_dir, jls_files[1])
    
    println("Loading simulation data from: $latest_file")
    
    # Deserialize the data
    aerobraking_outputs = open(latest_file, "r") do io
        deserialize(io)
    end
    
    return aerobraking_outputs
end

function create_mock_spacecraft_model()
    """
    Create a mock spacecraft model for testing when the original model is not available
    """
    # Create mock planet
    planet = config.Planet(
        Rp_e = 6378137.0,  # Earth radius in meters
        Rp_p = 6356752.0,
        mass = 5.972e24,
        μ = 3.986004418e14
    )
    
    # Create mock spacecraft body
    body = config.SpacecraftModel(
        joints = [],
        links = [],
        roots = [],
        n_reaction_wheels = 3,
        n_thrusters = 8,
        prop_mass = [100.0],
        laser_effector = false
    )
    
    # Create composite model
    model = (planet = planet, body = body)
    
    return model
end

function reconstruct_aerobraking_outputs(raw_data::Dict)
    """
    Reconstruct the aerobraking_outputs dictionary with proper structure for plots_multi
    """
    reconstructed = Dict{Any, Any}()
    
    for (obj_id, data) in raw_data
        println("Processing spacecraft/object ID: $obj_id")
        
        # Check if this is valid simulation data
        if !isa(data, Dict)
            println("  Warning: Invalid data format for object $obj_id")
            continue
        end
        
        # Create mock model if not present or invalid
        if !haskey(data, :m) || isnothing(data[:m])
            println("  Creating mock model for object $obj_id")
            data[:m] = create_mock_spacecraft_model()
        end
        
        # Ensure required keys exist with defaults
        reconstructed[obj_id] = Dict(
            :state => get(data, :state, 1),
            :m => data[:m],
            :name => get(data, :name, Dict(string(obj_id) => "Spacecraft_$obj_id")),
            :args => get(data, :args, Dict{Symbol, Any}()),
            :temp_name => get(data, :temp_name, tempname() * ".tmp")
        )
        
        # Ensure args has required fields
        args = reconstructed[obj_id][:args]
        if !haskey(args, :keplerian)
            args[:keplerian] = false
        end
        if !haskey(args, :closed_form)
            args[:closed_form] = 0
        end
        if !haskey(args, :body_shape)
            args[:body_shape] = "Spacecraft"
        end
        if !haskey(args, :type_of_mission)
            args[:type_of_mission] = "Orbital"
        end
        if !haskey(args, :orientation_sim)
            args[:orientation_sim] = false
        end
        
        println("  Successfully reconstructed data for object $obj_id")
    end
    
    return reconstructed
end

function create_mock_arrow_files(aerobraking_outputs::Dict, temp_dir::String = mktempdir())
    """
    Create mock Arrow files if they don't exist for testing purposes
    """
    println("Creating mock Arrow files in: $temp_dir")
    
    for (obj_id, data) in aerobraking_outputs
        # Create mock time series data
        n_points = 100
        times = collect(0.0:10.0:(n_points-1)*10.0)
        
        # Mock orbital data
        mock_data = DataFrame(
            time = times,
            pos_ii_1 = 7000e3 .+ 500e3 .* sin.(times ./ 1000),
            pos_ii_2 = 1000e3 .* cos.(times ./ 1000),
            pos_ii_3 = 500e3 .* sin.(times ./ 1500),
            vel_ii_1 = -3.0 .+ 0.5 .* cos.(times ./ 1000),
            vel_ii_2 = 7.5 .+ 0.3 .* sin.(times ./ 1000),
            vel_ii_3 = 0.1 .* cos.(times ./ 1500),
            q_1 = cos.(times ./ 2000),
            q_2 = sin.(times ./ 2000),
            q_3 = zeros(n_points),
            q_4 = zeros(n_points),
            omega_1 = 0.01 .* sin.(times ./ 500),
            omega_2 = 0.01 .* cos.(times ./ 500),
            omega_3 = 0.005 .* sin.(times ./ 1000),
            alt = 400e3 .+ 100e3 .* sin.(times ./ 1000),
            lat = 45.0 .+ 5.0 .* sin.(times ./ 2000),
            lon = -120.0 .+ 10.0 .* cos.(times ./ 2000)
        )
        
        # Create Arrow file
        filename = joinpath(temp_dir, "simulation_data_$(obj_id).arrow")
        Arrow.write(filename, mock_data)
        
        # Update temp_name to point to this directory
        aerobraking_outputs[obj_id][:temp_name] = filename
        
        println("  Created mock Arrow file: $filename")
    end
    
    return temp_dir
end

function test_plots_multi()
    """
    Main test function
    """
    println("="^50)
    println("Testing plot_data_multi functionality")
    println("="^50)
    
    try
        # Load simulation data
        println("\n1. Loading simulation data...")
        aerobraking_outputs = load_latest_simulation_data()
        println("   Loaded data for $(length(aerobraking_outputs)) objects")
        
        # Reconstruct data structure
        println("\n2. Reconstructing data structures...")
        reconstructed_data = reconstruct_aerobraking_outputs(aerobraking_outputs)
        
        # Create mock Arrow files if needed
        println("\n3. Setting up Arrow files...")
        temp_dir = create_mock_arrow_files(reconstructed_data)
        
        # Create global args for plotting
        args_global = Dict{Symbol, Any}(
            :results_dir => dirname(temp_dir),
            :plot_output => true,
            :save_plots => true
        )
        
        # Test plots_multi function
        println("\n4. Testing plots_multi function...")
        plots_multi(reconstructed_data, args_global)
        
        println("\n✅ Test completed successfully!")
        println("Temporary files created in: $temp_dir")
        
    catch e
        println("\n❌ Test failed with error:")
        println("Error type: $(typeof(e))")
        println("Error message: $e")
        println("\nStacktrace:")
        for (i, frame) in enumerate(stacktrace(catch_backtrace()))
            println("  $i. $frame")
        end
    end
end

function inspect_simulation_data(results_dir::String = joinpath(pwd(), "results"))
    """
    Inspect the structure of saved simulation data without running plots
    """
    println("Inspecting simulation data structure...")
    
    try
        data = load_latest_simulation_data(results_dir)
        
        println("\nTop-level keys: $(collect(keys(data)))")
        
        for (obj_id, obj_data) in data
            println("\nObject ID: $obj_id")
            if isa(obj_data, Dict)
                println("  Keys: $(collect(keys(obj_data)))")
                for (key, value) in obj_data
                    println("    $key: $(typeof(value))")
                end
            else
                println("  Type: $(typeof(obj_data))")
            end
        end
        
    catch e
        println("Error inspecting data: $e")
    end
end

# Run the test
if abspath(PROGRAM_FILE) == @__FILE__
    println("Choose an option:")
    println("1. Run full plot test")
    println("2. Just inspect data structure")
    
    choice = get(ENV, "TEST_MODE", "1")
    
    if choice == "2"
        inspect_simulation_data()
    else
        test_plots_multi()
    end
end