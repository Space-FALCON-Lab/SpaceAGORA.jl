using Test
include(joinpath(@__DIR__, "..", "II_examples", "targeted_40sat_sma_interlinks.jl"))

@testset "40-satellite animation frame data" begin
    args = targeted_40sat_configuration(save_results=false)
    model = args.interlink_model
    state = SimulationEngine.build_initial_conditions(args)
    distant = deepcopy(state)
    distant.sc[1].pos[1] += 1e7
    key = ((1, 1), (21, 1))
    append!(model.history, [
        InterLinkHistory(0.0, InterLinkKey[]),
        InterLinkHistory(1.0, [key]),
        InterLinkHistory(2.0, InterLinkKey[]),
    ])
    solution = (prob=(p=(args=args,),), t=[0.0, 1.0, 2.0], u=[state, state, distant])
    frames = interlink_animation_frames(solution; frame_count=3)
    @test frames.times == [0.0, 1.0, 2.0]
    @test length(frames.endpoints) == 400
    @test Set(Iterators.flatten(frames.endpoints)) == Set(1:40)
    @test all(1 <= a <= 20 && 21 <= b <= 40 for (a, b) in frames.endpoints)
    @test all(length(status) == 400 for status in frames.status)
    edge = model.graph.edge_index_lookup[key]
    @test [status[edge] for status in frames.status] == [1, 2, 0]
    @test all(status -> all(code -> code in 0:2, status), frames.status)
    @test !any(model.graph.available) # Postprocessing uses its own graph.
    @test !any(model.graph.active)
    @test length(model.history) == 3
    @test args.mission_configuration.mission_time ≈ 20frames.period
    @test interlink_animation_frames(solution; frame_count=10).times == [0.0, 1.0, 2.0]
    @test interlink_animation_frames(solution; frame_count=2).times == [0.0, 2.0]
    @test_throws ArgumentError interlink_animation_frames(solution; frame_count=1)
    missing_states = merge(solution, (; t=[0.0, 2.0], u=[state, distant]))
    @test_throws ErrorException interlink_animation_frames(missing_states; frame_count=3)
    @test_throws ArgumentError run_targeted_40sat_example(fps=0)
    @test_throws ArgumentError run_targeted_40sat_example(frame_count=1)
    mktempdir() do output
        path = save_interlink_animation(solution; output, frame_count=3, fps=5)
        @test basename(path) == "interlinks.gif"
        @test readdir(output) == ["interlinks.gif"]
        @test open(io -> String(read(io, 6)), path) == "GIF89a"
    end
end
