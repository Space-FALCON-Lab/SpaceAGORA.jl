using Test
using SpaceAGORA
using ComponentArrays
using RecursiveArrayTools: ArrayPartition
using StaticArrays

@testset "Callback state contract preserves both solver layouts" begin
    service = SpaceAGORA.CallbackServices
    engine = SpaceAGORA.SimulationEngine
    args = (dynamics_model=(spacecraft=[(dry_mass=8.0, prop_mass=2.0, links=[1, 2])],),)
    u = ComponentVector(sc=[(pos=[1.0, 2.0, 3.0], vel=[4.0, 5.0, 6.0],
        q=[0.0, 0.0, 0.0, 1.0], mass=9.0, heat_loads=[0.2, 0.3])])
    split = ArrayPartition(
        ComponentVector(sc=[(vel=[4.0, 5.0, 6.0],)]),
        ComponentVector(sc=[(pos=[1.0, 2.0, 3.0],)]))
    for state in (u, split)
        @test service._state_position_ii(state, 1) == SVector(1.0, 2.0, 3.0)
        @test service._state_velocity_ii(state, 1) == SVector(4.0, 5.0, 6.0)
        @test service._state_mass_kg(state, args, 1) == engine._state_mass_kg(state, args, 1)
        @test service._state_heat_loads(state, args, 1) == engine._state_heat_loads(state, args, 1)
        @test service._state_quaternion(state, 1) == engine._state_quaternion(state, 1)
        @test service._is_gravity_backbone_state(state) == (state === split)
    end
    @test service._state_mass_kg(u, args, 1) == 9.0
    @test service._state_mass_kg(split, args, 1) == 10.0
    @test service._state_quaternion(split, 1) === nothing
    changed = deepcopy(u)
    @test service._set_state_position_velocity_ii!(changed, 1,
        SVector(8.0, 9.0, 10.0), SVector(11.0, 12.0, 13.0)) === nothing
    @test changed.sc[1].pos == [8.0, 9.0, 10.0]
    @test changed.sc[1].vel == [11.0, 12.0, 13.0]
    @test changed.sc[1].mass == u.sc[1].mass
    @test changed.sc[1].heat_loads == u.sc[1].heat_loads
    @test_throws ArgumentError service._set_state_position_velocity_ii!(split, 1,
        SVector(8.0, 9.0, 10.0), SVector(11.0, 12.0, 13.0))
end

@testset "Callback declarations load without engine or physical models" begin
    root = dirname(dirname(pathof(SpaceAGORA)))
    isolated = Module(gensym(:CallbackContractFixture))
    Base.include(isolated, joinpath(root, "src", "orchestration", "simulation", "callback_services.jl"))
    contract = getproperty(isolated, :CallbackServices)
    @test !isdefined(isolated, :SimulationEngine)
    @test !isdefined(isolated, :SimulationModel)
    @test_throws MethodError contract.sample_environment(nothing)
    # The installed provider handles unsupported calls through the existing
    # engine method errors rather than a mutable service registry.
    @test_throws MethodError SpaceAGORA.CallbackServices.sample_environment(nothing)
end
