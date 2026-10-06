using Test
using Random

include(joinpath(@__DIR__, "..", "..", "examples", "Earth_RPO_CubeSat_MPC_Replanning.jl"))

@testset "Replanning obstacle endpoint clearance" begin
    rng = MersenneTwister(740)
    for separation in (0.01, 0.1, 1.5, 5.0), required in (0.25, 0.75), _ in 1:25
        start = SVector{3, Float64}(randn(rng, 3))
        goal = start + separation .* normalize(SVector{3, Float64}(randn(rng, 3)))
        demo = (
            initial_relative_state_rtn=vcat(start, zeros(3)),
            goal_rtn=goal,
            geometry=(station=(keepout_radius_m=0.25,), chaser=(half_extents_body=SVector(0.1, 0.1, 0.15),)),
        )
        # Overlapping exclusion regions and a collinear center exercise the old
        # alternating-correction failure, including deliberately short paths.
        center = 0.5 .* (start + goal)
        adjusted = _rpo_replanning_endpoint_safe_center(demo, center, 0.22, required)
        for endpoint in (start, goal)
            @test _rpo_replanning_sphere_endpoint_clearance(demo, adjusted, 0.22, endpoint) >= required
        end
        @test _rpo_replanning_endpoint_safe_center(demo, adjusted, 0.22, required) ≈ adjusted
        @test_throws ArgumentError rpo_replanning_case_config(
            demo, (id=:baseline_static_map,); safe_distance_m=0.25, obstacle_endpoint_clearance_m=0.1,
        )
    end
end
