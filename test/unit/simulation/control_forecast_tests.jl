module ControlForecastTests
using Test, SpaceAGORA, Distributed, Serialization
const fixture = normpath(joinpath(@__DIR__, "..", "..", "fixtures", "control_forecast_fixture.jl"))
Base.include(Main, fixture)
const F = Main.ControlForecastFixture
const SE = SpaceAGORA.SimulationEngine

struct RemoteBackedArray <: AbstractVector{Float64}
    channel::RemoteChannel
end
Base.size(::RemoteBackedArray) = (1,)
Base.getindex(::RemoteBackedArray, ::Int) = 1.0

@testset "isolated pre-command forecasts" begin
    @test_throws ArgumentError SE._check_forecast_data(RemoteBackedArray(RemoteChannel(()->Channel{Int}(1))))
    @test_throws ArgumentError SE._check_forecast_data(view([1.0], :))
    @test SE._check_forecast_data((data=[1.0, 2.0], mask=trues(2))) === nothing
    mktempdir() do root
        directory = joinpath(root,"run")
        F.parent(directory)
        snapshot = deserialize(joinpath(directory,"snapshot.bin"))
        before = Set(Distributed.workers())
        Base.ScopedValues.with(SE._CONTROL_FORECAST_RESTART=>(;)) do
            @test_throws ArgumentError capture_control_forecast(nothing,nothing,0.,nothing)
            @test_throws ArgumentError with_control_forecast_worker(identity,snapshot)
        end
        @test Set(Distributed.workers()) == before
        with_control_forecast_worker(snapshot) do evaluate
            worker = only(setdiff(Set(Distributed.workers()),before))
            # The fixture is test code, loaded explicitly into the disposable
            # worker. Production companion types use the `packages` argument.
            remotecall_fetch(Base.include,worker,Main,fixture)
            @test remotecall_fetch(F.child,worker,directory) === nothing
            @test remotecall_fetch(F.kernel_replay,worker) === nothing
            candidate = only(snapshot.configuration.control_model.control_effectors)
            @test_throws RemoteException evaluate(candidate;duration_s=-1.,sample_period_s=0.1)
        end
        @test Set(Distributed.workers()) == before
    end
end
end
