using Test, LinearAlgebra, SpaceAGORA

@testset "RPO corner-stop retiming fallback" begin
    guidance = SpaceAGORA.SimulationModel.GuidanceHooks
    points = [0.0 3.0 3.0; 0.0 0.0 4.0; 0.0 0.0 0.0]
    distances = guidance.rpo_arc_length_params(points)
    speed_limit = fill(0.1, 3)
    retime(; kwargs...) = guidance.rpo_hcw_acceleration_constrained_profile(
        points, distances, speed_limit, 0.1, 0.0011, 0.00625, 20_000;
        kwargs...,
    )
    @test retime(corner_stop_fallback=false) === nothing
    profile = retime()
    @test profile !== nothing
    if profile !== nothing
        @test guidance.rpo_hcw_feedforward_max_acceleration(profile, 0.0011) <= 0.00625
        reference, _, _ = guidance.rpo_retime_path_from_profile(profile)
        @test reference[:, 1] ≈ points[:, 1]
        @test reference[:, end] ≈ points[:, end]
        @test all(diff(profile.s_ref) .>= 0.0)
        @test all(p -> abs(p[2]) < 1.0e-12 || abs(p[1] - 3.0) < 1.0e-12, eachcol(reference))
        @test maximum(profile.speed_ref) <= maximum(speed_limit)
    end
    # A fallback must still respect the configured resource budget.
    @test guidance.rpo_hcw_acceleration_constrained_profile(
        points, distances, speed_limit, 0.1, 0.0011, 0.00625, 10,
    ) === nothing
    straight = points[:, 1:2]
    args = (straight, distances[1:2], speed_limit[1:2], 0.1, 0.0011, 0.00625, 20_000)
    @test guidance.rpo_hcw_acceleration_constrained_profile(args...).s_ref ==
        guidance.rpo_hcw_acceleration_constrained_profile(args...; corner_stop_fallback=false).s_ref
end

@testset "RPO rest-to-rest reference endpoints" begin
    guidance = SpaceAGORA.SimulationModel.GuidanceHooks
    points = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    profile = guidance.RPORetimingProfile(
        points, [0.0, 1.0], [0.0, 0.25, 0.5, 0.75, 1.0], fill(0.5, 5), 0.5,
    )
    times, positions, velocities = guidance.rpo_reference_from_profile(profile)
    @test iszero(velocities[:, 1])
    @test iszero(velocities[:, end])
    @test velocities[:, 2:end-1] ≈ repeat([0.5, 0.0, 0.0], 1, 3)
    @test guidance.rpo_hcw_feedforward_max_acceleration(profile, 0.0) ≈ 1.0
    @test times[end] == 2.0
    @test positions[:, end] == points[:, end]

    # Even a sub-timestep path has separate acceleration and braking intervals.
    short = guidance.rpo_profile_from_speed_envelope(
        points, [0.0, 1.0], [10.0, 10.0], 1.0, 10,
    )
    @test length(short.s_ref) == 3
    _, _, short_velocity = guidance.rpo_reference_from_profile(short)
    @test iszero(short_velocity[:, 1])
    @test iszero(short_velocity[:, end])
    @test short_velocity[1, 2] > 0.0
    @test guidance.rpo_hcw_feedforward_max_acceleration(short, 0.0) > 0.0
    @test guidance.rpo_profile_from_speed_envelope(
        points, [0.0, 1.0], [10.0, 10.0], 1.0, 1,
    ) === nothing
end
