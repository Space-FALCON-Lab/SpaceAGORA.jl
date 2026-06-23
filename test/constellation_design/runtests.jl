using Test
using SpaceAGORA.ConstellationDesign

"""
Run all constellation design tests.
"""
@testset "Constellation Design Tests" begin
    include("test_stage1_cached_regression.jl")
    include("test_stage2_cached_regression.jl")
    include("test_e2e_cached_regression.jl")
    include("test_cache_validation.jl")
    include("test_invariants.jl")
end
