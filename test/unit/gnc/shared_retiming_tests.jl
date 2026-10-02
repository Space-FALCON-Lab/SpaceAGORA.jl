using Test, LinearAlgebra, StaticArrays, Logging

module SharedRetimingOnly
    using LinearAlgebra, StaticArrays
    const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
    include(joinpath(ROOT, "src/gnc/shared/path_geometry.jl"))
    include(joinpath(ROOT, "src/gnc/shared/rpo/path_geometry.jl"))
    include(joinpath(ROOT, "src/gnc/shared/rpo/profile_evaluation.jl"))
    include(joinpath(ROOT, "src/gnc/shared/rpo/path_retiming.jl"))
    # Constant-clearance fixture supplies the existing shared geometry queries.
    rpo_clearance_to_station(p, geometry) = (clearance=geometry.clearance, distance=geometry.clearance)
    rpo_clearance_distance_to_station(p, geometry) = geometry.clearance
end

@testset "Shared retiming executes without HYPR policy or configuration" begin
    R = SharedRetimingOnly
    for name in (:RPOPSOConfig, :rpo_retime_available_distance,
                 :rpo_retime_pointwise_speed, :rpo_retime_sampling_ds_m)
        @test !isdefined(R, name)
    end
    geometry = (station=(keepout_radius_m=0.0,), chaser=(half_extents_body=[0.0,0.0,0.0],), clearance=100.0)
    samples = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    input = copy(samples)
    distances = Float64[]
    function available(c, d, safe)
        push!(distances, c)
        return d - safe
    end
    inputs = (max_speed_mps=1.0, min_speed_mps=0.0, dt_s=0.25, max_steps=20,
        available_distance=available, pointwise_speed=(d,k)->1.0)
    path, s, v = R.rpo_retime_samples(samples, geometry; inputs...)
    @test path[1,:] == [0.0,0.25,0.5,0.75,1.0]
    @test path[2:3,:] == zeros(2,5)
    @test s == [0.0,0.25,0.5,0.75,1.0]
    @test v == ones(5)
    @test samples == input
    @test distances == [100.0,100.0]
    slower = R.rpo_retime_samples(samples, geometry; merge(inputs,(pointwise_speed=(d,k)->0.5,))...)
    @test size(slower[1],2) == 9
    @test slower[1] != path
    duplicate = hcat(samples[:,1], samples)
    @test R.rpo_retime_samples(duplicate, geometry; inputs...) == (path,s,v)
    @test_throws ErrorException R.rpo_retime_samples(samples, geometry;
        merge(inputs,(pointwise_speed=(d,k)->error("policy failure"),))...)
    zero, zs, zv = with_logger(NullLogger()) do
        R.rpo_retime_samples(zeros(3,2), geometry; inputs...)
    end
    @test zero == zeros(3,1)
    @test zs == [0.0] && zv == [0.0]
    limited = with_logger(NullLogger()) do
        R.rpo_retime_samples(samples, geometry; merge(inputs,(max_steps=1,))...)
    end
    @test limited[1][:,end] == samples[:,end]
    @test length(limited[2]) == 3
    @test limited[3][end] == 1e-3

    # Analytic rest-to-rest motion: four metres at 1 m/s² gives a 4 s trip.
    line = [0.0 4.0; 0.0 0.0; 0.0 0.0]
    curve = R.RPORetimeCurve(line,:polyline)
    clearances = [NaN,NaN]
    profile_inputs = (max_speed_mps=Inf,min_speed_mps=0.0,initial_speed_mps=0.0,
        a_max_mps2=1.0,available_distance=(c,d,safe)->d,
        pointwise_speed=(d,k)->100.0,warn=false)
    profile = R.rpo_retime_profile(curve,line,Float64[],clearances,geometry; profile_inputs...)
    @test profile.s == [0.0,2.0,4.0]
    @test profile.v == [0.0,2.0,0.0]
    @test profile.t == [0.0,2.0,4.0]
    @test profile.a_seg == [1.0,-1.0]
    @test profile.duration_s == 4.0
    @test profile.length_m == 4.0
    @test profile.fallback_count == 0
    @test isequal(clearances,[NaN,NaN])
    @test line == [0.0 4.0;0.0 0.0;0.0 0.0]
    reference = R.rpo_retimed_reference_from_profile(profile,1.0)
    @test reference.r_rtn[1,:] == [0.0,0.5,2.0,3.5,4.0]
    @test reference.v_rtn[1,:] == [0.0,1.0,2.0,1.0,0.0]
    @test reference.t_s == collect(0.0:4.0)
    @test all(iszero,reference.v_rtn[:,end])
    @test_throws ErrorException R.rpo_retime_profile(curve,line,Float64[],clearances,geometry;
        merge(profile_inputs,(available_distance=(c,d,safe)->error("policy failure"),))...)
end

@testset "Shared retiming rejects definition-time HYPR coupling" begin
    R = SharedRetimingOnly
    source = read(joinpath(R.ROOT,"src/gnc/shared/rpo/path_retiming.jl"),String)
    target = Module(gensym(:SharedRetimeWitness))
    Core.eval(target, :(using LinearAlgebra, StaticArrays))
    Base.include(target,joinpath(R.ROOT,"src/gnc/shared/rpo/path_geometry.jl"))
    Base.include_string(target,source)
    @test !isdefined(target,:RPOPSOConfig)
    rejected = try
        Base.include_string(target,"forbidden_retime(cfg::RPOPSOConfig) = cfg")
        nothing
    catch e
        e
    end
    @test rejected isa LoadError && rejected.error isa UndefVarError && rejected.error.var === :RPOPSOConfig
end
