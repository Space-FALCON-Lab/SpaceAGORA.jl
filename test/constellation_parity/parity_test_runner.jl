#!/usr/bin/env julia

# parity_test_runner.jl
#
# Parity test runner for comparing SpaceAGORA constellation design outputs
# against CAPOConstellations outputs for clusters 1 and 2.
#
# Usage:
#   julia --project test/constellation_parity/parity_test_runner.jl

using Test
using YAML
using DataFrames, CSV
using Arrow
using LinearAlgebra
using Statistics

# Add parent directory to load path to access SpaceAGORA modules
push!(LOAD_PATH, joinpath(@__DIR__, "../../src"))
using SpaceAGORA.ConstellationDesign

const PARITY_TOLERANCE = 1e-2
const CLUSTER_DATA_DIR = joinpath(@__DIR__, "cluster_data")

"""
    compare_arrays(a::AbstractArray, b::AbstractArray; tol::Float64=PARITY_TOLERANCE) -> Bool

Compare two arrays with floating-point tolerance.
"""
function compare_arrays(a::AbstractArray, b::AbstractArray; tol::Float64=PARITY_TOLERANCE)
    size(a) == size(b) || return false
    max_diff = maximum(abs.(a .- b))
    return max_diff < tol
end

"""
    compare_dict_values(a::Any, b::Any; tol::Float64=PARITY_TOLERANCE) -> Bool

Compare two values with appropriate comparison method based on type.
"""
function compare_dict_values(a::Any, b::Any; tol::Float64=PARITY_TOLERANCE)
    if a isa AbstractArray && b isa AbstractArray
        return compare_arrays(a, b; tol=tol)
    elseif a isa AbstractFloat && b isa AbstractFloat
        return abs(a - b) < tol
    elseif a isa Integer && b isa Integer
        return a == b
    elseif a isa AbstractDict && b isa AbstractDict
        return compare_dicts(a, b; tol=tol)
    elseif a isa Vector && b isa Vector
        if length(a) != length(b)
            return false
        end
        return all(compare_dict_values(ai, bi; tol=tol) for (ai, bi) in zip(a, b))
    else
        return a == b
    end
end

"""
    compare_dicts(a::AbstractDict, b::AbstractDict; tol::Float64=PARITY_TOLERANCE) -> Bool

Compare two dictionaries with tolerance for floating-point values.
"""
function compare_dicts(a::AbstractDict, b::AbstractDict; tol::Float64=PARITY_TOLERANCE)
    keys_a = Set(keys(a))
    keys_b = Set(keys(b))
    
    if keys_a != keys_b
        @warn "Key mismatch: missing in a: $(setdiff(keys_b, keys_a)), missing in b: $(setdiff(keys_a, keys_b))"
        return false
    end
    
    for key in keys_a
        if !compare_dict_values(a[key], b[key]; tol=tol)
            @warn "Value mismatch for key $key: $(a[key]) vs $(b[key])"
            return false
        end
    end
    
    return true
end

"""
    test_stage0_parity(cluster_id::Int) -> Bool

Test Stage 0 stochastic seeding parity for a given cluster.
"""
function test_stage0_parity(cluster_id::Int)
    @testset "Stage 0 Seeding - Cluster $cluster_id" begin
        # Load CAPO reference results (if available)
        capo_results_path = joinpath(CLUSTER_DATA_DIR, "cluster_$(cluster_id)", "stage0_capo_results.jld2")
        
        if !isfile(capo_results_path)
            @warn "CAPO reference results not found for cluster $cluster_id, skipping parity test"
            return true
        end
        
        # Load SpaceAGORA configuration
        config_path = joinpath(CLUSTER_DATA_DIR, "cluster_$(cluster_id)", "config.yaml")
        
        if !isfile(config_path)
            @warn "Config file not found for cluster $cluster_id, skipping parity test"
            return true
        end
        
        # Run SpaceAGORA Stage 0
        spaceagora_results = ConstellationDesign.run_stage0_seeding(config_path)
        
        # Load CAPO results
        capo_results = JLD2.load(capo_results_path, "results")
        
        # Compare results
        @test compare_dicts(spaceagora_results, capo_results; tol=PARITY_TOLERANCE)
    end
    
    return true
end

"""
    test_stage1_parity(cluster_id::Int) -> Bool

Test Stage 1 constellation design parity for a given cluster.
"""
function test_stage1_parity(cluster_id::Int)
    @testset "Stage 1 Design - Cluster $cluster_id" begin
        # Load CAPO reference results (if available)
        capo_results_path = joinpath(CLUSTER_DATA_DIR, "cluster_$(cluster_id)", "stage1_capo_results.jld2")
        
        if !isfile(capo_results_path)
            @warn "CAPO reference results not found for cluster $cluster_id, skipping parity test"
            return true
        end
        
        # Load SpaceAGORA configuration
        config_path = joinpath(CLUSTER_DATA_DIR, "cluster_$(cluster_id)", "config.yaml")
        
        if !isfile(config_path)
            @warn "Config file not found for cluster $cluster_id, skipping parity test"
            return true
        end
        
        # Run SpaceAGORA Stage 1
        spaceagora_results = ConstellationDesign.run_constellation_design(config_path)
        
        # Load CAPO results
        capo_results = JLD2.load(capo_results_path, "results")
        
        # Compare results
        @test compare_dicts(spaceagora_results, capo_results; tol=PARITY_TOLERANCE)
    end
    
    return true
end

"""
    test_stage2_parity(cluster_id::Int) -> Bool

Test Stage 2 control verification parity for a given cluster.
"""
function test_stage2_parity(cluster_id::Int)
    @testset "Stage 2 Verification - Cluster $cluster_id" begin
        # Load CAPO reference results (if available)
        capo_results_path = joinpath(CLUSTER_DATA_DIR, "cluster_$(cluster_id)", "stage2_capo_results.jld2")
        
        if !isfile(capo_results_path)
            @warn "CAPO reference results not found for cluster $cluster_id, skipping parity test"
            return true
        end
        
        # Load SpaceAGORA configuration
        config_path = joinpath(CLUSTER_DATA_DIR, "cluster_$(cluster_id)", "config.yaml")
        
        if !isfile(config_path)
            @warn "Config file not found for cluster $cluster_id, skipping parity test"
            return true
        end
        
        # Run SpaceAGORA Stage 2
        spaceagora_results = ConstellationDesign.run_stage2_verification(config_path)
        
        # Load CAPO results
        capo_results = JLD2.load(capo_results_path, "results")
        
        # Compare results
        @test compare_dicts(spaceagora_results, capo_results; tol=PARITY_TOLERANCE)
    end
    
    return true
end

"""
    run_full_parity_test(cluster_ids::Vector{Int}=[1, 2]) -> Bool

Run full parity tests for all specified clusters.
"""
function run_full_parity_test(cluster_ids::Vector{Int}=[1, 2])
    @testset "Constellation Design Parity Tests" begin
        for cluster_id in cluster_ids
            @testset "Cluster $cluster_id" begin
                test_stage0_parity(cluster_id)
                test_stage1_parity(cluster_id)
                test_stage2_parity(cluster_id)
            end
        end
    end
    
    return true
end

"""
    main(args::Vector{String}=ARGS)

Main entry point for parity test runner.
"""
function main(args::Vector{String}=ARGS)
    # Parse arguments
    cluster_ids = if isempty(args)
        [1, 2]
    else
        parse.(Int, args)
    end
    
    println("Running parity tests for clusters: $(join(cluster_ids, ", "))")
    println("Tolerance: $PARITY_TOLERANCE")
    
    # Run tests
    success = run_full_parity_test(cluster_ids)
    
    if success
        println("\nAll parity tests passed!")
    else
        println("\nSome parity tests failed.")
    end
    
    return success
end

# Run main if this file is executed directly
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
