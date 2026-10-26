module BodyRateHistoryTests
using Test
using StaticArrays
using SpaceAGORA: SimulationModel
const SimulationConfiguration = SimulationModel.SimulationConfiguration
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "examples", "support", "aerobraking_result_tables.jl"))

# Inspect actual named definitions without importing a campaign or plot backend.
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
    length(definitions) == 1 || error("Expected exactly one definition of $(name).")
    return only(definitions)
end

module RenderBoundary
using ..BodyRateHistoryTests: _read_simulation_results, _body_rate_history_samples
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
Core.eval(RenderBoundary, source_function(joinpath(ROOT, "examples", "support", "aerobraking_mission_plot_utils.jl"),
                                          :_save_body_rate_history_plot))

module LegacyBoundary
scatter(; kwargs...) = (; kwargs...)
attr(; kwargs...) = (; kwargs...)
Layout(; kwargs...) = (; kwargs...)
plot(traces, layout) = (; traces, layout)
display(figure) = nothing
savefig(figure, args...; kwargs...) = figure
end
Core.eval(LegacyBoundary, source_function(joinpath(ROOT, "scripts", "plotting", "plot_data.jl"),
                                         :angular_velocity_plot))

module SavedProducer
using StaticArrays
using ..BodyRateHistoryTests: SimulationModel
const SC = SimulationModel.SimulationCallbacks
end
const SLEW = source_function(joinpath(ROOT, "test", "fixtures", "research_contracts", "save_field_recipes.jl"), :slew_save_fields)
# Read the attributed getter and SaveField expression unchanged. Other demo
# channels, telemetry loading and propagation are neither imported nor run.
const RATE_GETTER = only(filter(SLEW.args[2].args) do expr
    expr isa Expr && expr.head === :(=) && expr.args[1] isa Expr &&
        expr.args[1].head === :call && expr.args[1].args[1] === :ω_of
end)
function body_rate_field_calls(form)
    form isa Expr || return Expr[]
    if form.head === :call && form.args[1] == :(SC.SaveField) &&
       any(arg -> arg isa QuoteNode && arg.value === :body_rate_rad_s, form.args)
        return [form]
    end
    return reduce(vcat, (body_rate_field_calls(arg) for arg in form.args); init=Expr[])
end
Core.eval(SavedProducer, RATE_GETTER)
const RECORDED_FIELD = Core.eval(SavedProducer, only(body_rate_field_calls(SLEW)))

function fixture()
    DataFrame(sc2_body_rate_3=[-0.02, 0.0, 0.04, 0.05],
              sc1_body_rate_2=[0.01, -0.02, 0.03, -0.04], time=[120.0, 0.0, 60.0, 60.0],
              sc2_body_rate_1=[-0.0, 0.01, -0.01, 0.02],
              sc1_body_rate_1=[0.0, -0.0, 0.02, -0.03],
              sc1_body_rate_3=[-0.04, 0.05, 0.0, 0.01],
              sc2_body_rate_2=[0.03, -0.04, 0.05, -0.06])
end
function caught_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

@testset "existing body-rate producer and actual result-table writer" begin
    @test RECORDED_FIELD.name === :body_rate_rad_s
    @test RECORDED_FIELD.column_prefix == "body_rate"
    @test RECORDED_FIELD.per_satellite
    times = [60.0, 0.0]
    rates = [SVector(-0.0, 0.01, -0.02), SVector(0.03, -0.04, 0.05)]
    states = [(sc=[(ω=r,)],) for r in rates]
    snapshots = [Dict(RECORDED_FIELD.name => Base.invokelatest(RECORDED_FIELD.getter, u, t, nothing))
                 for (u, t) in zip(states, times)]
    args = (dynamics_model=(spacecraft=[(id=41,)],),)
    df = SimulationModel.IOOutputs._build_results_dataframe(times, snapshots, [RECORDED_FIELD], args)
    @test names(df) == ["time", "sc1_body_rate_1", "sc1_body_rate_2", "sc1_body_rate_3"]
    samples = _body_rate_history_samples(df; spacecraft_index=1)
    @test samples.time_min == [1.0, 0.0]
    for i in 1:3
        @test isequal(samples.rates_rad_s[i], [r[i] for r in rates])
    end
    @test_throws ArgumentError _body_rate_history_samples(df; spacecraft_index=41)

    # Generalized explicit save recipe: positional spacecraft indices, not IDs.
    field = SimulationModel.SaveField(:body_rate_rad_s,
        (u, t, integrator) -> [SVector{3, Float64}(sc.ω) for sc in u.sc];
        per_satellite=true, column_prefix="body_rate")
    multi = [(sc=[(ω=r,), (ω=-r,)],) for r in rates]
    data = [Dict(field.name => field.getter(u, t, nothing)) for (u, t) in zip(multi, times)]
    args2 = (dynamics_model=(spacecraft=[(id=41,), (id=97,)],),)
    table = SimulationModel.IOOutputs._build_results_dataframe(times, data, [field], args2)
    @test _body_rate_history_samples(table; spacecraft_index=2).rates_rad_s ==
          ntuple(i -> [-r[i] for r in rates], 3)
    @test_throws ArgumentError _body_rate_history_samples(table; spacecraft_index=97)
end

@testset "body-rate preparation preserves raw data and legacy conversion" begin
    df = fixture()
    before = deepcopy(df)
    for index in (1, 2)
        samples = _body_rate_history_samples(df; spacecraft_index=index)
        @test samples.time_min == [2.0, 0.0, 1.0, 1.0]
        @test samples.spacecraft_index == index
        @test samples.labels == ("ω1", "ω2", "ω3")
        @test samples.colors == (:red, :green, :blue)
        for i in 1:3
            stored = df[!, Symbol("sc$(index)_body_rate_$(i)")]
            @test isequal(samples.rates_rad_s[i], stored)
            @test samples.rates_rad_s[i] !== stored
            @test isequal(samples.components[i], rad2deg.(stored) .* 60)
        end
    end
    @test isequal(df, before)
    samples = _body_rate_history_samples(df; spacecraft_index=1)
    samples.rates_rad_s[1][1] = 1.0
    samples.components[2][1] = 1.0
    @test isequal(df, before)
    df.sc2_body_rate_1 = fill(missing, nrow(df))
    @test _body_rate_history_samples(df; spacecraft_index=1).rates_rad_s[3] == before.sc1_body_rate_3
    singleton = _body_rate_history_samples(before[3:3, :]; spacecraft_index=1)
    @test singleton.time_min == [1.0]
    @test singleton.rates_rad_s == ([0.02], [0.03], [0.0])
    exact = DataFrame(time=[0, 60], sc1_body_rate_1=[2^53 + 1, -(2^53 + 1)],
                      sc1_body_rate_2=BigInt[big(2)^100 + 1, -(big(2)^100 + 1)],
                      sc1_body_rate_3=[1//3, -2//3])
    raw = _body_rate_history_samples(exact; spacecraft_index=1)
    for i in 1:3
        stored = exact[!, Symbol("sc1_body_rate_$(i)")]
        @test isequal(raw.rates_rad_s[i], stored)
        @test eltype(raw.rates_rad_s[i]) === eltype(stored)
        @test isequal(raw.components[i], rad2deg.(stored) .* 60)
    end
end

@testset "selected body-rate schema fails clearly" begin
    for index in (0, -1, true, 1.5, "1", missing)
        err = caught_error(() -> _body_rate_history_samples(fixture(); spacecraft_index=index))
        @test err isa ArgumentError
        @test occursin("spacecraft_index", sprint(showerror, err))
    end
    @test_throws UndefKeywordError _body_rate_history_samples(fixture())
    @test_throws ArgumentError _body_rate_history_samples(fixture(); spacecraft_index=3)
    for column in (:time, :sc1_body_rate_1, :sc1_body_rate_2, :sc1_body_rate_3)
        err = caught_error(() -> _body_rate_history_samples(select(fixture(), Not(column)); spacecraft_index=1))
        @test err isa ArgumentError
        @test occursin(string(column), sprint(showerror, err))
        for bad in (missing, "0.5", true, 1 + 0im, NaN, Inf, -Inf)
            df = fixture()
            df[!, column] = Any[df[!, column]...]
            df[2, column] = bad
            err = caught_error(() -> _body_rate_history_samples(df; spacecraft_index=1))
            @test err isa ArgumentError
            @test occursin("$(column), row 2", sprint(showerror, err))
        end
    end
    @test_throws ArgumentError _body_rate_history_samples(fixture()[1:0, :]; spacecraft_index=1)
    huge_time = fixture(); huge_time.time[1] = floatmax(Float64)
    @test isfinite(_body_rate_history_samples(huge_time; spacecraft_index=1).time_min[1])
    alias = DataFrame(time=[0.0], sc1_omega_1=[1.0], sc1_omega_2=[2.0], sc1_omega_3=[3.0])
    @test_throws ArgumentError _body_rate_history_samples(alias; spacecraft_index=1)
    for i in 1:3
        df = fixture()
        column = Symbol("sc1_body_rate_$(i)")
        df[2, column] = floatmax(Float64)
        err = caught_error(() -> _body_rate_history_samples(df; spacecraft_index=1))
        @test err isa ArgumentError
        @test occursin("$(column), row 2", sprint(showerror, err))
        @test occursin("after conversion", sprint(showerror, err))
    end
end

@testset "body-rate facade matches actual legacy inputs and protects CSV" begin
    mktempdir() do directory
        csv_path = joinpath(directory, "simulation_results.csv")
        output_path = joinpath(directory, "body_rate.png")
        df = fixture()
        CSV.write(csv_path, df)
        input_bytes = read(csv_path)
        result = Base.invokelatest(() -> RenderBoundary._save_body_rate_history_plot(
            csv_path; spacecraft_index=2, output_path))
        @test length(result.plot.traces) == 3
        @test result.plot.saves == [abspath(output_path)]
        @test result.output_path == output_path
        @test result.spacecraft_index == 2
        @test read(csv_path) == input_bytes
        @test !isfile(output_path)
        legacy = Base.invokelatest(() -> LegacyBoundary.angular_velocity_plot("unused",
            DataFrame(time=df.time, omega_1=df.sc2_body_rate_1,
                      omega_2=df.sc2_body_rate_2, omega_3=df.sc2_body_rate_3)))
        for i in 1:3
            trace = result.plot.traces[i]
            @test isequal(trace.x, legacy.traces[i].x)
            @test isequal(trace.y, legacy.traces[i].y)
            @test trace.label == legacy.traces[i].name
            @test string(trace.color) == legacy.traces[i].line.color
            @test trace.y === result.components[i]
        end
        @test result.plot.attributes.xlabel == legacy.layout.xaxis_title
        @test result.plot.attributes.ylabel == legacy.layout.yaxis_title
        write(output_path, "existing destination")
        figure_count = length(RenderBoundary.FIGURES)
        invalid = joinpath(directory, "partial.csv")
        CSV.write(invalid, select(df, Not(:sc2_body_rate_3)))
        empty = joinpath(directory, "empty.csv")
        CSV.write(empty, df[1:0, :])
        overflow = joinpath(directory, "overflow.csv")
        bad = deepcopy(df); bad.sc2_body_rate_1[1] = floatmax(Float64)
        CSV.write(overflow, bad)
        for input in (invalid, empty, overflow, joinpath(directory, "absent.csv"))
            @test_throws ArgumentError Base.invokelatest(() -> RenderBoundary._save_body_rate_history_plot(
                input; spacecraft_index=2, output_path))
        end
        for dest in ("", " ", joinpath(directory, "no_extension"), csv_path)
            @test_throws ArgumentError Base.invokelatest(() -> RenderBoundary._save_body_rate_history_plot(
                csv_path; spacecraft_index=2, output_path=dest))
        end
        link = joinpath(directory, "hardlink.png")
        hardlink(csv_path, link)
        aliases = [link]
        if !Sys.iswindows()
            symlink_path = joinpath(directory, "symlink.png")
            symlink(csv_path, symlink_path)
            push!(aliases, symlink_path, joinpath("~", relpath(symlink_path, homedir())))
        end
        for dest in aliases
            err = caught_error(() -> Base.invokelatest(() -> RenderBoundary._save_body_rate_history_plot(
                csv_path; spacecraft_index=2, output_path=dest)))
            @test err isa ArgumentError
            @test occursin("filesystem aliases", sprint(showerror, err))
            @test read(csv_path) == input_bytes
        end
        @test_throws UndefKeywordError Base.invokelatest(() -> RenderBoundary._save_body_rate_history_plot(csv_path; output_path))
        @test_throws UndefKeywordError Base.invokelatest(() -> RenderBoundary._save_body_rate_history_plot(csv_path; spacecraft_index=1))
        @test length(RenderBoundary.FIGURES) == figure_count
        @test read(output_path, String) == "existing destination"
        @test read(csv_path) == input_bytes
    end
end
end # module BodyRateHistoryTests
