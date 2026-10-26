module InterlinkScoringTests

using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel.SchedulingPolicies: matching_workspace!

function circular_target_fixture(; shared_target_terminal=false, shared_partner=false,
    eccentricity=0.0, second_power=299_792_458.0)
    spacecraft = [(n_terminal=2,) for _ in 1:3]
    state = (sc=[
        (pos=SVector(1.0, 0.0, 0.0), vel=SVector(0.0, sqrt(1 + eccentricity), 0.0), mass=1.0),
        (pos=SVector(1.0, -1.0, 0.0), vel=zeros(SVector{3, Float64}), mass=1.0),
        (pos=SVector(1.0, 1.0, 0.0), vel=zeros(SVector{3, Float64}), mass=1.0),
    ],)
    model = ConstellationInterLinkModel()
    register_allowed_link!(model, spacecraft, (1, 1), (2, 1);
        parameters=InterLinkParameters(P=299_792_458.0, B=1.0, range=10.0))
    register_allowed_link!(model, spacecraft, (1, shared_target_terminal ? 1 : 2),
        (shared_partner ? 2 : 3, 1);
        parameters=InterLinkParameters(P=second_power, B=1.0, range=10.0))
    update_availability!(model, spacecraft, state)
    return model, state
end

@testset "Near-circular eccentricity objective admission" begin
    model, state = circular_target_fixture()
    target = state.sc[1]
    first_force = force_on_endpoint(target, state.sc[2], model.graph.edge_properties[1])
    second_force = force_on_endpoint(target, state.sc[3], model.graph.edge_properties[2])
    # With mu=mass=radius=speed=1, a tangential force gives e_dot=(2F_y,0,0).
    # The separate positive norms sum to four, while the combined force is zero.
    @test first_force == SVector(0.0, 1.0, 0.0)
    @test second_force == -first_force
    @test norm(2first_force) + norm(2second_force) == 4.0
    @test iszero(first_force + second_force)
    for targets in (1, :both)
        policy = SchedulingPolicyModel(EccentricityRatePolicy(; targets))
        @test_throws ArgumentError score_candidates!(model, policy, state, 1.0;
            candidate_force=force_on_endpoint)
    end
    near_model, near_state = circular_target_fixture(eccentricity=1e-10)
    @test_throws ArgumentError score_candidates!(near_model, SchedulingPolicyModel(EccentricityRatePolicy()),
        near_state, 1.0; candidate_force=force_on_endpoint)

    # Distinct terminals are harmless when terminal matching forbids concurrent
    # links, a link is unavailable, or its contribution is zero.
    for configuration in ((shared_target_terminal=true,), (shared_partner=true,), (second_power=0.0,))
        supported, supported_state = circular_target_fixture(; configuration...)
        policy = SchedulingPolicyModel(EccentricityRatePolicy())
        @test isnothing(score_candidates!(supported, policy, supported_state, 1.0;
            candidate_force=force_on_endpoint))
        @test length(select_interlinks!(supported, policy)) == 1
    end
    model.graph.available[2] = false
    policy = SchedulingPolicyModel(EccentricityRatePolicy())
    @test isnothing(score_candidates!(model, policy, state, 1.0; candidate_force=force_on_endpoint))
    @test length(select_interlinks!(model, policy)) == 1

    # The projected eccentricity rate away from circularity remains additive.
    eccentric_model, eccentric_state = circular_target_fixture(eccentricity=1e-4)
    @test isnothing(score_candidates!(eccentric_model, policy, eccentric_state, 1.0;
        candidate_force=force_on_endpoint))
    @test matching_workspace!(eccentric_model, policy).score[1] > 0
    @test matching_workspace!(eccentric_model, policy).score[2] < 0

    # SMA policies keep their existing linear objective for the same graph.
    model, state = circular_target_fixture()
    for objective in (SemimajorAxisRatePolicy(), SemimajorAxisRatePolicy(targets=:lower_index))
        sma_policy = SchedulingPolicyModel(objective)
        @test isnothing(score_candidates!(model, sma_policy, state, 1.0;
            candidate_force=force_on_endpoint))
        @test matching_workspace!(model, sma_policy).score == [2.0, -2.0]
        @test length(select_interlinks!(model, sma_policy)) == 1
    end
end

end # module InterlinkScoringTests
