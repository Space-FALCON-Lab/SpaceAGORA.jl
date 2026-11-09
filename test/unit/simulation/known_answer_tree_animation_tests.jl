module KnownAnswerTreeAnimationTests

using Test
using CSV
using DataFrames
include(joinpath(@__DIR__, "..", "..", "..", "examples", "animate_known_answer_tree_interlinks.jl"))
const A = KnownAnswerTreeAnimation

function fixture(directory)
    spacecraft = DataFrame(index=1:4, role=["helper", "helper", "target", "target"],
        mass_kg=[180.0, 260.0, 200.0, 300.0], altitude_m=[1e6, 1e6, 1.1e6, 1.1e6],
        terminal_capacity=ones(Int, 4))
    candidates = DataFrame(pair=["1-2", "1-3", "1-4", "2-3", "2-4", "3-4"],
        a_spacecraft=[1, 1, 1, 2, 2, 3], b_spacecraft=[2, 3, 4, 3, 4, 4],
        a_terminal=ones(Int, 6), b_terminal=ones(Int, 6),
        allowed=[false, true, false, true, true, false], range_m=fill(200e3, 6))
    samples = DataFrame(time_s=[0.0, 2.0, 4.0],
        available=["1-3;2-3;2-4", "2-3", "empty"], active=["1-3;2-4", "2-3", "empty"],
        selected_total_m_per_s=[0.034, 0.027, 0.0],
        exhaustive_best_m_per_s=[0.034, 0.027, 0.0],
        largest_edge_m_per_s=[0.027, 0.027, 0.0], maximum_terminal_use=[1, 1, 0],
        excluded_1_4_distance_m=[103e3, 102e3, 201e3],
        excluded_1_4_geometrically_feasible=[true, true, false], excluded_1_4_active=falses(3),
        geometrically_feasible_1_2=trues(3), geometrically_feasible_1_3=[true, false, false],
        geometrically_feasible_1_4=[true, true, false], geometrically_feasible_2_3=[true, true, false],
        geometrically_feasible_2_4=[true, false, false], geometrically_feasible_3_4=trues(3))
    CSV.write(joinpath(directory, "spacecraft_before_run.csv"), spacecraft)
    CSV.write(joinpath(directory, "candidates_before_run.csv"), candidates)
    CSV.write(joinpath(directory, "accepted_step_comparison.csv"), samples)
    return samples
end

@testset "Known-answer graph animation" begin
    mktempdir() do directory
        @test_throws ErrorException A.read_recording(directory)
        samples = fixture(directory)
        recording = A.read_recording(directory)
        @test recording.candidates.pair == ["1-2", "1-3", "1-4", "2-3", "2-4", "3-4"]
        frames = A.animation_frames(recording; frame_count=5)
        @test getproperty.(frames, :time_s) == [0, 1, 2, 3, 4]
        @test getproperty.(frames, :sample_time_s) == [0, 0, 2, 2, 4]
        @test getproperty.(frames, :status) == [
            [-1, 2, -1, 1, 2, -1], [-1, 2, -1, 1, 2, -1],
            [-1, 0, -1, 2, 0, -1], [-1, 0, -1, 2, 0, -1], [-1, 0, 0, 0, 0, -1]]
        @test [count(==(-1), frame.status) for frame in frames] == [3, 3, 3, 3, 2]
        @test A.NODE_POSITIONS[1][1] == A.NODE_POSITIONS[4][1] > A.NODE_POSITIONS[2][1] == A.NODE_POSITIONS[3][1]
        @test A.NODE_POSITIONS[1][2] == A.NODE_POSITIONS[2][2] < A.NODE_POSITIONS[3][2] == A.NODE_POSITIONS[4][2]
        figure = A.graph_frame(recording, first(frames))
        @test length(figure.series_list) == 16 # Six edges, four nodes, three legend swatches, timeline.
        for (index, pair) in enumerate(A.ALL_PAIRS)
            edge = only(eachrow(recording.candidates[recording.candidates.pair .== pair, :]))
            a, b = A.NODE_POSITIONS[edge.a_spacecraft], A.NODE_POSITIONS[edge.b_spacecraft]
            @test figure.series_list[index][:x] == [a[1], b[1]]
            @test figure.series_list[index][:y] == [a[2], b[2]]
            code = first(frames).status[index]
            color = code == -1 ? A.EXCLUDED_COLOR : A.COLORS[code + 1]
            @test figure.series_list[index][:linecolor] == A.Plots.plot_color(color)
            @test figure.series_list[index][:linestyle] == :solid
        end
        final_figure = A.graph_frame(recording, last(frames))
        @test final_figure.series_list[3][:linestyle] == :dash
        @test final_figure.series_list[3][:linecolor] == final_figure.series_list[2][:linecolor]
        annotations = final_figure.subplots[1][:annotations]
        @test any(annotation -> annotation[3].str == "Not allowed / feasible (2)", annotations)
        @test any(annotation -> annotation[3].str == "Unavailable (4)", annotations)
        for (index, progress) in ((1, 0.0), (3, 0.5), (5, 1.0))
            timed = index == 1 ? figure : A.graph_frame(recording, frames[index])
            @test timed.series_list[end - 2][:x] == [140, 1060]
            @test timed.series_list[end - 1][:x] == [140, 140 + 920progress]
            @test timed.series_list[end][:x] == [140 + 920progress]
            @test timed.series_list[end][:y] == [738]
            text = "Simulation time: $(Float64(frames[index].time_s)) / 4.0 s"
            @test any(annotation -> annotation[3].str == text, timed.subplots[1][:annotations])
        end
        @test getproperty.(frames, :selected_total_m_per_s) == [0.034, 0.034, 0.027, 0.027, 0.0]
        @test length(A.animation_frames(recording; frame_count=2)) == 2
        @test_throws ArgumentError A.animation_frames(recording; frame_count=1)
        @test_throws ArgumentError A.save_animation(input=directory, fps=0)
        @test_throws ArgumentError A.matching_pairs("1-3;1-3")

        rendered = A.save_animation(input=directory, output=joinpath(directory, "rendered"),
            frame_count=3, fps=2)
        @test rendered.frame_count == 3 && rendered.fps == 2
        @test open(io -> String(read(io, 6)), rendered.gif_path) == "GIF89a"
        @test open(io -> read(io, 8), rendered.preview_path) == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
        @test filesize(rendered.gif_path) > 1000

        for pair in ("1-2", "3-4")
            altered = copy(samples)
            altered[2, Symbol("geometrically_feasible_", replace(pair, "-" => "_"))] = false
            CSV.write(joinpath(directory, "accepted_step_comparison.csv"), altered)
            changing = A.read_recording(directory)
            changing_frames = A.animation_frames(changing; frame_count=3)
            index = findfirst(==(pair), A.ALL_PAIRS)
            @test first(changing_frames).status[index] == -1
            @test changing_frames[2].status[index] == 0
            @test last(changing_frames).status[index] == -1
            unavailable = A.graph_frame(changing, changing_frames[2])
            @test unavailable.series_list[index][:linestyle] == :dash
            @test unavailable.series_list[index][:linecolor] == A.Plots.plot_color(A.COLORS[1])
        end
        for (field, value) in ((:active, "1-4"), (:available, "1-4"),
            (:selected_total_m_per_s, NaN), (:time_s, 0.0),
            (:maximum_terminal_use, 2), (:excluded_1_4_active, true),
            (:excluded_1_4_geometrically_feasible, false),
            (:geometrically_feasible_1_3, true), (:geometrically_feasible_1_4, false))
            invalid = copy(samples)
            invalid[2, field] = value
            CSV.write(joinpath(directory, "accepted_step_comparison.csv"), invalid)
            @test_throws ArgumentError A.read_recording(directory)
        end
        CSV.write(joinpath(directory, "accepted_step_comparison.csv"), select(samples, Not(A.GEOMETRY_COLUMNS)))
        @test_throws ArgumentError A.read_recording(directory)
        fixture(directory)
        candidates = copy(recording.candidates)
        deleteat!(candidates, findfirst(==("1-2"), candidates.pair))
        CSV.write(joinpath(directory, "candidates_before_run.csv"), candidates)
        @test_throws ArgumentError A.read_recording(directory)
    end
end

end # module KnownAnswerTreeAnimationTests
