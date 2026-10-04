using Test, LinearAlgebra, StaticArrays, Logging

module SharedRetimingOnly
    using LinearAlgebra, StaticArrays
    const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
    include(joinpath(ROOT, "src/gnc/shared/path_geometry.jl"))
    include(joinpath(ROOT, "src/gnc/shared/rpo/path_geometry.jl"))
    include(joinpath(ROOT, "src/gnc/shared/rpo/profile_evaluation.jl"))
    include(joinpath(ROOT, "src/gnc/shared/rpo/path_retiming.jl"))
    # Constant-clearance fixture supplies the existing shared geometry queries.
    rpo_clearance_to_station(p, geometry) = (clearance=geometry.clearance, distance=geometry.clearance + geometry.station.keepout_radius_m + maximum(geometry.chaser.half_extents_body))
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

@testset "Shared retiming forwards geometry and preserves speed-policy semantics" begin
    R = SharedRetimingOnly
    # Clearance excludes the 3 m body/station margin; distance includes it.
    geometry = (station=(keepout_radius_m=1.0,), chaser=(half_extents_body=[2.0,1.0,1.0],), clearance=2.0)
    samples = [0.0 1.0 2.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    calls = Tuple{Float64,Float64,Float64}[]
    available(c,d,safe) = (push!(calls,(c,d,safe)); d-c-safe)
    inputs = (max_speed_mps=1.0, min_speed_mps=0.0, dt_s=0.4, max_steps=40,
        available_distance=available, pointwise_speed=(d,k)->d, safe_distance_m=0.5)
    path = R.rpo_retime_samples(samples,geometry; inputs...)
    @test calls == fill((2.0,5.0,0.5),3)
    # The policy returns 2.5 m/s. max_speed_mps only caps fallback, not that result.
    @test path[2] == [0.0,1.0,2.0]
    @test path[3] == fill(2.5,3)

    empty!(calls)
    duplicate = samples[:,[1,2,2,3]]
    @test R.rpo_retime_samples(duplicate,geometry; inputs...) == path
    @test calls == fill((2.0,5.0,0.5),3)
    floor_path = R.rpo_retime_samples(samples,geometry;
        merge(inputs,(pointwise_speed=(d,k)->0.25,min_speed_mps=0.5,dt_s=1.0))...)
    @test floor_path[3] == fill(0.5,5)

    logger = Test.TestLogger()
    fallback_path = with_logger(logger) do
        R.rpo_retime_samples(samples,geometry;
            merge(inputs,(pointwise_speed=(d,k)->0.0,max_speed_mps=0.125,
                          fallback_speed_mps=0.5,dt_s=1.0))...)
    end
    @test fallback_path[3] == fill(0.125,17)
    @test length(logger.logs) == 1
    warning = only(logger.logs)
    @test warning.level == Logging.Warn
    @test occursin("infeasible zero-speed samples",warning.message)
    # Existing tuple-style logging has a source-derived key; pin its values.
    fields = only(values(warning.kwargs))
    @test fields.count == 3
    @test fields.first_idx == 1
    @test fields.fallback_speed_mps == 0.125

    line = [0.0 2.0 4.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    curve = R.RPORetimeCurve(line,:polyline)
    inputs_p = (max_speed_mps=1.0,min_speed_mps=0.0,initial_speed_mps=1.0,
        a_max_mps2=1.0,available_distance=available,pointwise_speed=(d,k)->d,
        safe_distance_m=0.5,warn=false)
    empty!(calls)
    profile = R.rpo_retime_profile(curve,line,Float64[],fill(NaN,3),geometry; inputs_p...)
    @test calls == fill((2.0,5.0,0.5),3)
    @test profile.v_point == fill(2.5,3)
    @test profile.v == [1.0,2.0,0.0]
    @test profile.clearance == fill(2.0,3)

    fallback_inputs = merge(inputs_p,(pointwise_speed=(d,k)->0.0,
        initial_speed_mps=0.0,max_speed_mps=0.125,fallback_speed_mps=0.5))
    fallback = R.rpo_retime_profile(curve,line,Float64[],fill(2.0,3),geometry; fallback_inputs...)
    @test fallback.v_point == [0.0,0.125,0.0]
    @test fallback.v == [0.0,0.125,0.0]
    @test fallback.fallback_count == 1
    floored = R.rpo_retime_profile(curve,line,Float64[],fill(2.0,3),geometry;
        merge(fallback_inputs,(min_speed_mps=0.5,))...)
    @test floored.v_point == [0.0,0.5,0.0]
    @test floored.v == [0.0,0.5,0.0]
    @test floored.fallback_count == 1

    # A right-angle polyline has discrete curvature 2sqrt(2) at the interior.
    elbow = [0.0 1.0 1.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
    curvatures = Float64[]
    curved = R.rpo_retime_profile(R.RPORetimeCurve(elbow,:polyline),elbow,
        Float64[],fill(2.0,3),geometry;
        merge(inputs_p,(pointwise_speed=(d,k)->(push!(curvatures,k);1.0),))...)
    @test curved.curvature ≈ fill(2sqrt(2.0),3) rtol=8eps(Float64)
    @test curvatures == curved.curvature

    # A straight Bezier has exact 4 m length, even through its quadrature route.
    endpoints = line[:,[1,3]]
    bezier = R.rpo_retime_profile(R.RPORetimeCurve(endpoints,:bezier),endpoints,
        [0.0,1.0],[2.0,2.0],geometry;
        merge(inputs_p,(initial_speed_mps=0.0,pointwise_speed=(d,k)->100.0))...)
    @test bezier.bezier
    @test bezier.params == [0.0,0.5,1.0]
    @test bezier.s ≈ [0.0,2.0,4.0] rtol=8eps(Float64)
    @test bezier.length_m ≈ 4.0 rtol=8eps(Float64)
    @test bezier.duration_s ≈ 4.0 rtol=8eps(Float64)
    @test bezier.v ≈ [0.0,2.0,0.0] rtol=8eps(Float64)
    @test bezier.curvature == zeros(3)
end
