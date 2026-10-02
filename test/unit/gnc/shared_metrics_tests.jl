using Test

# This module deliberately has no planner/configuration types. Only the two
# scalar kernels are exercised here; other metrics retain their own dependencies.
module SharedMetricsOnly
    using LinearAlgebra, StaticArrays
    const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
    include(joinpath(ROOT, "src", "gnc", "shared", "rpo", "path_metrics.jl"))
end

@testset "Shared metric kernels need no HYPR configuration" begin
    K = SharedMetricsOnly
    @test !isdefined(K, :RPOPSOConfig)
    @test !isdefined(K, :RPOPSOConfigurator)
    points = [0.0 3.0; 0.0 4.0; 0.0 0.0]
    inputs = (cost_ref_distance_m=0.0, sample_ds_m=0.2, tf_s=10.0,
        mass_kg=2.0, isp_s=4.0, g0_mps2=5.0)
    @test K.rpo_path_cost_normalization_refs(points; inputs...) ===
        (straight_len=5.0, len_ref=5.0, fuel_ref=0.05)
    @test K.rpo_path_cost_normalization_refs(points;
        merge(inputs, (cost_ref_distance_m=10.0,))...) ===
        (straight_len=5.0, len_ref=10.0, fuel_ref=0.1)
    @test K.rpo_path_cost_normalization_refs(zeros(3,2);
        merge(inputs, (mass_kg=0.0,))...).fuel_ref === 1e-12
    floor = K.rpo_path_cost_normalization_refs(zeros(3,2);
        cost_ref_distance_m=0.0, sample_ds_m=0.0, tf_s=0.0,
        mass_kg=2.0, isp_s=0.0, g0_mps2=0.0)
    @test floor.straight_len === 0.0
    @test floor.len_ref === 1e-6
    @test floor.fuel_ref ≈ 2e9 rtol=4eps(Float64)
    fuel_inputs = (tf_s=2.0, mass_kg=2.0, isp_s=4.0, g0_mps2=5.0)
    quadratic = [0.0 1.0 4.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    @test K.rpo_fuel_proxy_from_samples(quadratic; fuel_inputs...) === 0.2
    @test K.rpo_fuel_proxy_from_samples(zeros(3,5); fuel_inputs...) === 0.0
    for columns in 0:2
        @test K.rpo_fuel_proxy_from_samples(zeros(3,columns); fuel_inputs...) === 0.0
    end
    # Short inputs return before the three-coordinate finite-difference loop.
    @test K.rpo_fuel_proxy_from_samples(zeros(2,2); fuel_inputs...) === 0.0
    @test_throws BoundsError K.rpo_path_cost_normalization_refs(zeros(3,0); inputs...)
    @test_throws BoundsError K.rpo_path_cost_normalization_refs(zeros(2,2); inputs...)
    @test_throws UndefKeywordError K.rpo_path_cost_normalization_refs(points)
    @test_throws UndefKeywordError K.rpo_fuel_proxy_from_samples(quadratic)
end

using SpaceAGORA

metric_bits(x::Float64) = reinterpret(UInt64, x)
metric_bits(x::NamedTuple) = map(metric_bits, x)

@testset "Existing metric access and HYPR forwarding remain compatible" begin
    G = SpaceAGORA.SimulationModel.GuidanceHooks
    M = SpaceAGORA.SimulationModel
    K = SharedMetricsOnly
    @test parentmodule(G.rpo_path_cost_normalization_refs) === G
    @test parentmodule(G.rpo_fuel_proxy_from_samples) === G
    configs = [M.RPOPSOConfig(),
        M.RPOPSOConfig(cost_ref_distance_m=0.0, sample_ds_m=0.17, tf_s=7.3,
            mass_kg=13.2, isp_s=44.5, g0_mps2=9.78),
        M.RPOPSOConfig(cost_ref_distance_m=23.4, sample_ds_m=0.03, tf_s=18.7,
            mass_kg=6.1, isp_s=61.7, g0_mps2=3.71),
        M.RPOPSOConfig(cost_ref_distance_m=0.0, sample_ds_m=0.0, tf_s=0.0,
            mass_kg=0.0, isp_s=0.0, g0_mps2=0.0)]
    paths = [reshape([2.0,-1.0,0.5],3,1),
        [2.0 2.0; -1.0 -1.0; 0.5 0.5],
        [0.0 3.0; 0.0 4.0; 0.0 0.0],
        [0.0 1.0 4.0; 0.0 0.0 2.0; 0.0 1.0 1.0],
        [0.0 0.0 1.0 4.0 9.0; 2.0 2.0 -1.0 3.0 3.0; 1.0 1.0 0.0 2.0 -2.0]]
    for cfg in configs
        normalization_inputs = (cost_ref_distance_m=cfg.cost_ref_distance_m,
            sample_ds_m=cfg.sample_ds_m, tf_s=cfg.tf_s, mass_kg=cfg.mass_kg,
            isp_s=cfg.isp_s, g0_mps2=cfg.g0_mps2)
        fuel_inputs = (tf_s=cfg.tf_s, mass_kg=cfg.mass_kg,
            isp_s=cfg.isp_s, g0_mps2=cfg.g0_mps2)
        for points in paths
            expected = K.rpo_path_cost_normalization_refs(points; normalization_inputs...)
            @test metric_bits(G.rpo_path_cost_normalization_refs(points,cfg)) == metric_bits(expected)
            @test metric_bits(G.rpo_path_cost_normalization_refs(points; normalization_inputs...)) == metric_bits(expected)
            fuel = K.rpo_fuel_proxy_from_samples(points; fuel_inputs...)
            @test metric_bits(G.rpo_fuel_proxy_from_samples(points,cfg)) == metric_bits(fuel)
            @test metric_bits(G.rpo_fuel_proxy_from_samples(points; fuel_inputs...)) == metric_bits(fuel)
        end
        @test G.rpo_fuel_proxy_from_samples(zeros(3,0),cfg) === 0.0
        @test_throws BoundsError G.rpo_path_cost_normalization_refs(zeros(3,0),cfg)
    end
    @test_throws MethodError G.rpo_fuel_proxy_from_samples(["x" "y";"z" "q";"a" "b"], configs[1])
    @test_throws MethodError G.rpo_path_cost_normalization_refs([1.0,2.0,3.0], configs[1])
end
