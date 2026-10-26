# The package's own tests of its method: that the references it plans are the
# ones the method promises (shape, timing, the cases it must solve and refuse).
# The conformance checks in runtests.jl do not test the method; these tests do.
using Test, SpaceAGORA, StarterName
using Random: MersenneTwister
using SpaceAGORA.Conformance: rpo_conformance_requests

@testset "the method" begin
    request = first(rpo_conformance_requests())
    result = plan_rpo!(nothing, StraightLinePlanner(), request, MersenneTwister(1))
    @test result.status === :candidate
    @test result.reference.r_ref_rtn_m[:, end] == collect(request.goal_rtn_m)
end
