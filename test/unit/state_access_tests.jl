using Test

import SpaceAGORA
using RecursiveArrayTools: ArrayPartition
using StaticArrays

const _SE = SpaceAGORA.SimulationEngine

@testset "state access accepts flat solver states" begin
    u = @SVector [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    @test _SE._state_position_ii(u, 1) == SVector{3, Float64}(1.0, 2.0, 3.0)
    @test _SE._state_velocity_ii(u, 1) == SVector{3, Float64}(4.0, 5.0, 6.0)
    @test_throws BoundsError _SE._state_position_ii(u, 2)

    second_order_u = ArrayPartition(
        SVector{3, Float64}(4.0, 5.0, 6.0),
        SVector{3, Float64}(1.0, 2.0, 3.0),
    )
    @test _SE._state_position_ii(second_order_u, 1) == SVector{3, Float64}(1.0, 2.0, 3.0)
    @test _SE._state_velocity_ii(second_order_u, 1) == SVector{3, Float64}(4.0, 5.0, 6.0)
end
