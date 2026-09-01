using Test

@testset "Unit Test Placeholder" begin
    @test true
end

include("rpo_port_tests.jl")
include(joinpath(@__DIR__, "parallel", "cost_work_counts_tests.jl"))
include(joinpath(@__DIR__, "parallel", "cost_robust_timing_tests.jl"))
include(joinpath(@__DIR__, "parallel", "cost_machine_calibration_tests.jl"))
include(joinpath(@__DIR__, "parallel", "outer_route_persistence_tests.jl"))
include(joinpath(@__DIR__, "simulation", "jac_prototype_tests.jl"))
include(joinpath(@__DIR__, "robotics", "runtests.jl"))
