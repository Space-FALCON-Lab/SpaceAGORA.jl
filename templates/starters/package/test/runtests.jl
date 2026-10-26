using Test, SpaceAGORA, StarterName
using SpaceAGORA.Conformance: test_algorithm, rpo_conformance_requests, rpo_trusted_queries, rpo_conformance_simulation

@testset "StarterName" begin
    # Interface conformance with the RPO family: the family's representative
    # requests, its trusted validation queries and a short simulation. A pass
    # certifies the interface, not the method.
    test_algorithm(StraightLinePlanner(), rpo_conformance_requests(); trusted=rpo_trusted_queries(),
        simulation=rpo_conformance_simulation())

    include("method_tests.jl")
end
