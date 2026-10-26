module RawQuaternionHistoryTests
using Test
using SpaceAGORA: SimulationModel
const SimulationConfiguration = SimulationModel.SimulationConfiguration
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "examples", "support", "aerobraking_result_tables.jl"))

# Read only the named function definition. Plotting and legacy dispatcher imports
# are deliberately excluded; the production function bodies remain unmodified.
function source_function(path, name)
    definitions = Expr[]
    for form in Meta.parseall(read(path, String)).args
        form isa Expr || continue
        if form.head === :macrocall && form.args[1] === GlobalRef(Core, Symbol("@doc"))
            form = form.args[end]
        end
        form isa Expr && form.head === :function || continue
        signature = form.args[1]
        signature isa Expr && signature.head === :call && signature.args[1] === name || continue
        push!(definitions, form)
    end
    length(definitions) == 1 || error("Expected exactly one definition of $(name) in $(path).")
    return only(definitions)
end

module RenderBoundary
using ..RawQuaternionHistoryTests: _read_simulation_results, _quaternion_history_samples
const FIGURES = Any[]
function plot(; kwargs...)
    figure = (attributes=(; kwargs...), traces=Any[], saves=String[])
    push!(FIGURES, figure)
    return figure
end
function plot!(figure, x, y; kwargs...)
    push!(figure.traces, (; x, y, kwargs...))
    return figure
end
savefig(figure, path) = push!(figure.saves, path)
display(args...) = error("The explicit saver must not open a plot window.")
end
Core.eval(RenderBoundary, source_function(
    joinpath(ROOT, "examples", "support", "aerobraking_mission_plot_utils.jl"),
    :_save_quaternion_history_plot))

module LegacyBoundary
scatter(; kwargs...) = (; kwargs...)
attr(; kwargs...) = (; kwargs...)
Layout(; kwargs...) = (; kwargs...)
plot(traces, layout) = (; traces, layout)
display(figure) = nothing
savefig(figure, args...; kwargs...) = figure
end
Core.eval(LegacyBoundary, source_function(
    joinpath(ROOT, "scripts", "plotting", "plot_data.jl"), :quaternion_plot))

function fixture()
    # Deliberately scramble column order, times and components. Include identity,
    # sign-flipped identity, a non-unit quaternion, and a zero quaternion.
    DataFrame(sc2_q_4=[8.0, 9.0, 10.0, 11.0],
              sc1_q_2=[0.0, 0.0, -3.0, 0.0],
              time=[120.0, 0.0, 60.0, 60.0],
              sc1_q_4=[1.0, -1.0, -5.0, 0.0],
              sc2_q_1=[1.0, -2.0, 3.0, -4.0],
              sc1_q_1=[0.0, -0.0, 2.0, 0.0],
              sc2_q_3=[5.0, -6.0, 7.0, -8.0],
              sc1_q_3=[0.0, 0.0, 4.0, 0.0],
              sc2_q_2=[-9.0, 10.0, -11.0, 12.0])
end

function caught_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

@testset "raw quaternion preparation preserves saved data" begin
    df = fixture()
    before = deepcopy(df)
    samples = _quaternion_history_samples(df; spacecraft_index=1)
    @test samples.time_min == [2.0, 0.0, 1.0, 1.0]
    @test samples.spacecraft_index == 1
    @test samples.labels == ("q1", "q2", "q3", "q4")
    @test samples.colors == (:red, :green, :blue, :orange)
    @test isequal(samples.components, ([0.0, -0.0, 2.0, 0.0], [0.0, 0.0, -3.0, 0.0],
                                      [0.0, 0.0, 4.0, 0.0], [1.0, -1.0, -5.0, 0.0]))
    @test isequal(df, before)
    for i in 1:4
        @test samples.components[i] !== df[!, Symbol("sc1_q_$(i)")]
    end
    samples.components[1][1] = 99.0
    @test isequal(df, before)
    other = _quaternion_history_samples(df; spacecraft_index=2)
    @test other.components == ([1.0, -2.0, 3.0, -4.0], [-9.0, 10.0, -11.0, 12.0],
                              [5.0, -6.0, 7.0, -8.0], [8.0, 9.0, 10.0, 11.0])
    @test other.time_min == [2.0, 0.0, 1.0, 1.0]
    df.sc2_q_1 = fill(missing, 4)
    @test _quaternion_history_samples(df; spacecraft_index=1).components[4] == before.sc1_q_4
    singleton = _quaternion_history_samples(before[3:3, :]; spacecraft_index=1)
    @test singleton.time_min == [1.0]
    @test singleton.components == ([2.0], [-3.0], [4.0], [-5.0])

    # Float64 coercion would silently change these valid saved component values.
    exact = DataFrame(time=[0, 60], sc1_q_1=[2^53 + 1, -(2^53 + 1)],
                      sc1_q_2=BigInt[big(2)^100 + 1, -(big(2)^100 + 1)],
                      sc1_q_3=[1//3, -2//3], sc1_q_4=BigFloat[1, -1])
    raw = _quaternion_history_samples(exact; spacecraft_index=1)
    for i in 1:4
        column = exact[!, Symbol("sc1_q_$(i)")]
        @test isequal(raw.components[i], column)
        @test eltype(raw.components[i]) === eltype(column)
    end
end

@testset "selected quaternion columns fail clearly" begin
    for index in (0, -1, true, 1.5, "1", missing)
        err = caught_error(() -> _quaternion_history_samples(fixture(); spacecraft_index=index))
        @test err isa ArgumentError
        @test occursin("spacecraft_index", sprint(showerror, err))
    end
    @test_throws UndefKeywordError _quaternion_history_samples(fixture())
    @test_throws ArgumentError _quaternion_history_samples(fixture(); spacecraft_index=3)
    for column in (:time, :sc1_q_1, :sc1_q_2, :sc1_q_3, :sc1_q_4)
        df = select(fixture(), Not(column))
        err = caught_error(() -> _quaternion_history_samples(df; spacecraft_index=1))
        @test err isa ArgumentError
        @test occursin(string(column), sprint(showerror, err))
        for bad in (missing, "0.5", true, 1 + 0im, NaN, Inf, -Inf)
            df = fixture()
            df[!, column] = Any[df[!, column]...]
            df[2, column] = bad
            err = caught_error(() -> _quaternion_history_samples(df; spacecraft_index=1))
            @test err isa ArgumentError
            @test occursin("$(column), row 2", sprint(showerror, err))
        end
    end
    err = caught_error(() -> _quaternion_history_samples(fixture()[1:0, :]; spacecraft_index=1))
    @test err isa ArgumentError
    @test occursin("at least one saved row", sprint(showerror, err))
end

@testset "CSV facade and exact legacy trace comparison" begin
    mktempdir() do directory
        csv_path = joinpath(directory, "simulation_results.csv")
        output_path = joinpath(directory, "raw_quaternion.png")
        df = fixture()
        CSV.write(csv_path, df)
        input_bytes = read(csv_path)
        result = Base.invokelatest(() -> RenderBoundary._save_quaternion_history_plot(
            csv_path; spacecraft_index=2, output_path))
        @test length(result.plot.traces) == 4
        @test result.plot.saves == [output_path]
        @test result.output_path == output_path
        @test result.spacecraft_index == 2
        @test result.plot.attributes.xlabel == "Time [min]"
        @test result.plot.attributes.ylabel == "Quaternion"
        @test read(csv_path) == input_bytes
        @test !isfile(output_path) # The recording renderer never writes a file.

        # Execute the actual retained legacy function with a schema-only adapter.
        legacy_table = DataFrame(time=df.time, q_1=df.sc2_q_1, q_2=df.sc2_q_2,
                                 q_3=df.sc2_q_3, q_4=df.sc2_q_4)
        legacy = Base.invokelatest(() -> LegacyBoundary.quaternion_plot("unused", legacy_table))
        for i in 1:4
            trace = result.plot.traces[i]
            @test isequal(trace.x, legacy.traces[i].x)
            @test isequal(trace.y, legacy.traces[i].y)
            @test trace.label == legacy.traces[i].name
            @test string(trace.color) == legacy.traces[i].line.color
            @test trace.y === result.components[i]
        end
        @test result.plot.attributes.xlabel == legacy.layout.xaxis_title
        @test result.plot.attributes.ylabel == legacy.layout.yaxis_title

        # Invalid requests must fail before constructing a figure or saving over
        # an existing destination. Parsing is the same real CSV boundary.
        write(output_path, "existing destination")
        figure_count = length(RenderBoundary.FIGURES)
        invalid_path = joinpath(directory, "partial.csv")
        CSV.write(invalid_path, select(df, Not(:sc2_q_3)))
        empty_path = joinpath(directory, "empty.csv")
        CSV.write(empty_path, df[1:0, :])
        malformed_path = joinpath(directory, "malformed.csv")
        malformed = copy(df)
        malformed.sc2_q_1 = ["1", "bad", "3", "4"]
        CSV.write(malformed_path, malformed)
        for path in (invalid_path, empty_path, malformed_path, joinpath(directory, "absent.csv"))
            @test_throws ArgumentError Base.invokelatest(() -> RenderBoundary._save_quaternion_history_plot(
                path; spacecraft_index=2, output_path))
        end
        for destination in ("", "  ", joinpath(directory, "extension_missing"), csv_path)
            @test_throws ArgumentError Base.invokelatest(() -> RenderBoundary._save_quaternion_history_plot(
                csv_path; spacecraft_index=2, output_path=destination))
        end
        # Existing links to the source are still the input file, even when the
        # destination has a supported image suffix. Never hand them to savefig.
        hardlink_path = joinpath(directory, "hardlink.png")
        hardlink(csv_path, hardlink_path)
        alias_paths = [hardlink_path]
        if !Sys.iswindows() # Windows symlink creation may require privileges.
            symlink_path = joinpath(directory, "symlink.png")
            symlink(csv_path, symlink_path)
            push!(alias_paths, symlink_path)
            # Plots expands a leading `~` before writing. The alias guard must
            # inspect the same destination even when its literal path is absent.
            push!(alias_paths, joinpath("~", relpath(symlink_path, homedir())))
        end
        for alias_path in alias_paths
            err = caught_error(() -> Base.invokelatest(() -> RenderBoundary._save_quaternion_history_plot(
                csv_path; spacecraft_index=2, output_path=alias_path)))
            @test err isa ArgumentError
            @test occursin("filesystem aliases", sprint(showerror, err))
            @test read(csv_path) == input_bytes
        end
        @test_throws UndefKeywordError Base.invokelatest(() -> RenderBoundary._save_quaternion_history_plot(
            csv_path; output_path))
        @test_throws UndefKeywordError Base.invokelatest(() -> RenderBoundary._save_quaternion_history_plot(
            csv_path; spacecraft_index=1))
        @test length(RenderBoundary.FIGURES) == figure_count
        @test read(output_path, String) == "existing destination"
        @test read(csv_path) == input_bytes
    end
end
end # module RawQuaternionHistoryTests
