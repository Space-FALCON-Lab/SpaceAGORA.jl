using Test

@testset "Unit Test Placeholder" begin
    @test true
end

include("state_access_tests.jl")
include("rpo_port_tests.jl")
include("rpo_nonlinear_plant_tests.jl")
include("rpo_capsule_tests.jl")
include("rpo_corner_retiming_tests.jl")
include(joinpath(@__DIR__, "robotics", "runtests.jl"))
