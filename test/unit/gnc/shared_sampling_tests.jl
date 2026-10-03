using Test

# Load only shared geometry and sampling into a module with no planner types.
# This fixture is deliberately evaluated before this file imports SpaceAGORA.
module SharedSamplingOnly
    using LinearAlgebra, StaticArrays
    const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
    include(joinpath(ROOT, "src", "gnc", "shared", "path_geometry.jl"))
    include(joinpath(ROOT, "src", "gnc", "shared", "rpo", "path_geometry.jl"))
    include(joinpath(ROOT, "src", "gnc", "shared", "rpo", "sampling_settings.jl"))
    include(joinpath(ROOT, "src", "gnc", "shared", "rpo", "path_sampling.jl"))
    # Analytic one-point station used only to exercise the standalone sampler.
    rpo_clearance_distance_to_station(p, g) = norm(p) - g.station.keepout_radius_m - maximum(g.chaser.half_extents_body)
end

@testset "Shared sampling loads and runs without HYPR configuration" begin
    K = SharedSamplingOnly
    @test !isdefined(K, :RPOPSOConfig)
    @test !isdefined(K, :RPOPSOConfigurator)
    geometry = (station=(keepout_radius_m=0.25,), chaser=(half_extents_body=(0.05,0.05,0.05),))
    settings = K.RPOAdaptiveSamplingSettings()
    @test K.rpo_adaptive_sampling_min_ds_m(0.4, geometry, settings; safe_distance_m=0.2) == 0.1
    points = [-1.0 1.0; 1.0 1.0; 0.0 0.0]
    expected = [-1.0 -0.5 0.0 0.5 1.0; 1.0 1.0 1.0 1.0 1.0; 0.0 0.0 0.0 0.0 0.0]
    fixed = K.RPOAdaptiveSamplingSettings(enabled=false)
    for curve in (:bezier, :polyline)
        samples, params, clearances = K.rpo_sample_path_with_params(points, fixed, geometry; base_ds_m=0.5, curve_type=curve)
        @test samples == expected
        @test all(isnan, clearances)
        @test curve == :bezier ? params == [0.0,0.25,0.5,0.75,1.0] : isempty(params)
        adaptive = K.rpo_sample_path(points, settings, geometry; base_ds_m=0.2, curve_type=curve)
        @test adaptive[:,1] == points[:,1]
        @test adaptive[:,end] == points[:,end]
        @test all(isfinite, adaptive)
    end
end

using SpaceAGORA, StaticArrays

sampling_bits(x::AbstractArray{<:AbstractFloat}) = (size(x), reinterpret(UInt64, vec(x)))
sampling_bits(x::Tuple) = sampling_bits.(x)
sampling_bits(x) = x

@testset "Existing HYPR sampling calls preserve defaults and overrides" begin
    G = SpaceAGORA.SimulationModel.GuidanceHooks
    M = SpaceAGORA.SimulationModel
    geometry = M.RPOReferenceGeometry(M.RPOStationGeometry([0.0 0.0;0.0 2.0;0.0 0.0];keepout_radius_m=0.2); chaser=M.RPOCubeSatGeometry(dims_m=(0.1,0.1,0.1)))
    paths = [reshape([2.0,1.0,0.0],3,1), [2.0 2.0;1.0 1.0;0.0 0.0], [-2.0 2.0;1.0 1.0;0.0 0.0], [-2.0 0.0 2.0;1.0 0.1 1.0;0.0 0.5 0.0], [-2.0 -2.0 0.0 2.0;1.0 1.0 0.4 1.0;0.0 0.0 0.5 0.0]]
    @test parentmodule(G.RPOAdaptiveSamplingSettings) === G
    @test fieldnames(G.RPOAdaptiveSamplingSettings) == (:enabled, :max_ds_m, :far_clearance_m, :power, :safe_distance_fraction, :obstacle_guard_fraction)
    for enabled in (false,true), curve in (:bezier,:polyline)
        cfg = G.RPOPSOConfig(adaptive_sampling_enable=enabled, curve_type=curve, sample_ds_m=0.17, safe_distance_m=0.23, adaptive_sampling_max_ds_m=0.41, adaptive_sampling_far_clearance_m=1.2, adaptive_sampling_power=0.7, adaptive_sampling_safe_distance_fraction=0.4, adaptive_sampling_obstacle_guard_fraction=0.6)
        # Construct the expected shared inputs independently of the adapter.
        settings = G.RPOAdaptiveSamplingSettings(enabled=enabled, max_ds_m=0.41, far_clearance_m=1.2, power=0.7, safe_distance_fraction=0.4, obstacle_guard_fraction=0.6)
        @test (@inferred G._rpo_sampling_settings(cfg)) === settings
        for points in paths, override in (false,true)
            oldopts = override ? (safe_distance_m=0.11,base_ds_m=0.31,curve_type=(curve==:bezier ? :polyline : :bezier)) : NamedTuple()
            inputs = merge((safe_distance_m=0.23,base_ds_m=0.17,curve_type=curve),oldopts)
            for sample in (G.rpo_sample_path,G.rpo_sample_path_with_params)
                @test isequal(sampling_bits(sample(points,cfg,geometry;oldopts...)), sampling_bits(sample(points,settings,geometry;inputs...)))
            end
            old_simple = override ? (safe_distance_m=0.11,base_ds_m=0.31) : NamedTuple()
            simple = merge((safe_distance_m=0.0,base_ds_m=0.17),old_simple)
            for sample in (G.rpo_sample_path_polyline_adaptive,G.rpo_sample_path_bezier_adaptive,G.rpo_sample_path_bezier_adaptive_with_params)
                @test isequal(sampling_bits(sample(points,geometry,cfg;old_simple...)), sampling_bits(sample(points,geometry,settings;simple...)))
            end
        end
        for ds in (-1.0,0.0,0.2), safe in (-0.1,0.0,0.3)
            @test G.rpo_adaptive_sampling_min_ds_m(ds,geometry,cfg;safe_distance_m=safe) === G.rpo_adaptive_sampling_min_ds_m(ds,geometry,settings;safe_distance_m=safe)
        end
        for sample in (G.rpo_sample_path,G.rpo_sample_path_with_params)
            @test_throws ArgumentError sample(paths[3],cfg,geometry;curve_type=:unknown)
            @test_throws ArgumentError sample(paths[3],settings,geometry;base_ds_m=0.17,curve_type=:unknown)
        end
    end
end
