# Run by the optional-interlink-installation job after explicitly loading PythonCall.
module InterlinkPythonTests
using Test, SpaceAGORA, Random
using ..InterlinkTests: matching_workspace!, exhaustive_matching_score, gve_sma_configuration
const PythonCall = Base.get_extension(SpaceAGORA, :SpaceAGORAPythonCallExt).PythonCall
const SM = SpaceAGORA.SimulationModel

@testset "Explicit Python interlink backend" begin
    spacecraft = [SpacecraftModel(n_terminal=2) for _ in 1:5]
    triangle = ConstellationInterLinkModel()
    for (i, j) in ((1, 2), (2, 3), (1, 3))
        register_allowed_link!(triangle, spacecraft, (i, 1), (j, 1))
    end
    triangle.graph.available .= true
    triangle_policy = SchedulingPolicyModel(SemimajorAxisRatePolicy())
    matching_workspace!(triangle, triangle_policy).score .= [3.0, 2.0, 1.0]
        @test select_interlinks!(triangle, triangle_policy) == [((1, 1), (2, 1))]
        @test any(triangle_policy.matching.cyclic)
        PythonCall.GIL.@unlock begin
            @test_throws ArgumentError select_interlinks!(triangle, triangle_policy)
        end
        @test select_interlinks!(triangle, triangle_policy) == [((1, 1), (2, 1))]
        rng = MersenneTwister(718)
        pairs = [(i, j) for i in 1:5 for j in (i + 1):6]
        six = [SpacecraftModel() for _ in 1:6]
        for _ in 1:40
            random_model = ConstellationInterLinkModel()
            for (i, j) in shuffle(rng, pairs)[1:8]
                register_allowed_link!(random_model, six, (i, 1), (j, 1))
            end
            random_model.graph.available .= rand(rng, Bool, 8)
            random_policy = SchedulingPolicyModel(SemimajorAxisRatePolicy())
            scores = matching_workspace!(random_model, random_policy).score
            scores .= rand(rng, [-1.0, 0.0, 1e-20, 1.0, nextfloat(1.0), 3.0], 8)
            expected = exhaustive_matching_score(random_model, random_policy)
            selected = select_interlinks!(random_model, random_policy)
            actual = sum((Rational{BigInt}(scores[random_model.graph.edge_index_lookup[key]])
                for key in selected); init=zero(Rational{BigInt}))
            @test actual == expected
            @test length(Set(endpoint for key in selected for endpoint in key)) == 2 * length(selected)
            @test selected == select_interlinks!(random_model, random_policy)
        end
        # Directed links in a cycle keep their transmitter-first keys, and of a
        # two-way pair only the better direction is a candidate.
        three = [SpacecraftModel(n_terminal=1) for _ in 1:3]
        directed = InterlinkModel()
        beam = PowerBeamLink(power_w=1.0, range_m=1.0)
        for (i, j) in ((2, 1), (3, 2), (1, 3), (1, 2))
            register_link!(directed, three, (i, 1), (j, 1), beam)
        end
        directed.graph.available .= true
        directed_policy = SchedulingPolicyModel(PowerDeliveryPolicy())
        matching_workspace!(directed, directed_policy).score .= [1.0, 2.0, 1.0, 3.0]
        @test select_interlinks!(directed, directed_policy) == [((1, 1), (2, 1))]
        @test any(directed_policy.matching.cyclic)
        matching_workspace!(directed, directed_policy).score .= [1.0, 4.0, 1.0, 0.5]
        @test select_interlinks!(directed, directed_policy) == [((3, 1), (2, 1))]
        matching_workspace!(directed, directed_policy).score .= [3.0, 1.0, 5.0, 1.0]
        @test select_interlinks!(directed, directed_policy) == [((1, 1), (3, 1))]
        @test directed.graph.active == [false, false, true, false]
        matching_workspace!(directed, directed_policy).score .= [5.0, 1.0, 1.0, 2.0]   # the first of the pair wins
        @test select_interlinks!(directed, directed_policy) == [((2, 1), (1, 1))]
        multi_target = gve_sma_configuration(duration=20.0)
        register_allowed_link!(multi_target.interlink_model, multi_target.dynamics_model.spacecraft,
            (2, 1), (3, 1); parameters=InterLinkParameters(range=500e3))
        multi_target = SM.SimConfig._with_configuration(multi_target;
            scheduling_policy_model=SchedulingPolicyModel(SemimajorAxisRatePolicy(targets=[3, 1, 2])))
        multi_solution = run_simulation(multi_target; return_solution=true)
        multi_ran = multi_solution.prob.p.args
        @test multi_solution.t[end] == 20.0
        @test all(isfinite, multi_solution.u[end])
        @test any(multi_ran.scheduling_policy_model.matching.cyclic)
        @test length(multi_ran.interlink_model.history) == multi_solution.stats.naccept + 1
        @test all(length(sample.active) == 1 for sample in multi_ran.interlink_model.history)

end
end
