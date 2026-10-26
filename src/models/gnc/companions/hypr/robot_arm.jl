# Core configuration and compatibility contracts. Execution is optional.
include(joinpath(@__DIR__, "robot_arm", "config.jl"))
include(joinpath(@__DIR__, "robot_arm", "rrt_types.jl"))
include(joinpath(@__DIR__, "robot_arm", "execution_contracts.jl"))
