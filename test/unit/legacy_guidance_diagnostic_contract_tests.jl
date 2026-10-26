module LegacyGuidanceDiagnosticContractTests
using Test
using SpaceAGORA: SimulationModel

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const Guidance = SimulationModel.GuidanceHooks

# Load single definitions from retained sources without running a simulation.
function source_function(path, name)
    definitions = Expr[]
    for form in Meta.parseall(read(joinpath(ROOT, path), String)).args
        form isa Expr || continue
        if form.head === :macrocall && form.args[1] === GlobalRef(Core, Symbol("@doc"))
            form = form.args[end]
        end
        form isa Expr && form.head === :function || continue
        signature = form.args[1]
        signature isa Expr && signature.head === :call && signature.args[1] === name || continue
        push!(definitions, form)
    end
    return only(definitions)
end

module SavedSwitchBoundary
using ..LegacyGuidanceDiagnosticContractTests: SimulationModel
const SaveField = SimulationModel.SaveField
end
Core.eval(SavedSwitchBoundary, source_function("examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl", :_edg_runtime_control))
const SWITCH_FIELD_SOURCE = only(filter(Meta.parseall(read(joinpath(ROOT,
    "examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl"), String)).args) do x
    x isa Expr && x.head === :(=) && x.args[1] === :targeting_switch_save_field
end)
Core.eval(SavedSwitchBoundary, SWITCH_FIELD_SOURCE)

@testset "explicit saved switch telemetry and its reset lifetime" begin
    state = Guidance.AerobrakingEnergyDepletionState(num_sats=2)
    field = SavedSwitchBoundary.targeting_switch_save_field
    @test field.name === :edg_targeting_switch_s
    @test field.per_satellite
    @test field.column_prefix == "edg_targeting_switch_s"
    @test state.targeting_switch_s == [Inf, Inf]
    @test !hasproperty(state, :lamv_list)
    @test !hasproperty(state, :lambda_switch_list)
    state.targeting_switch_s .= [140.0, 160.0]
    state.heat_load_switches_s .= [(125.0, 130.0), (145.0, 150.0)]
    state.last_switch_solve_t .= [123.0, 124.0]
    state.energy_bracketing_count .= [2, 3]
    state.last_alpha_rad .= [0.1, 1.2]
    args = (dynamics_model=(spacecraft=[(id=41,), (id=97,)],),
        control_model=(control_effectors=[(state=state,)],))
    integrator = (p=(args=args,),)
    snapshot1 = field.getter(nothing, 123.0, integrator)
    Guidance._edg_reset_pass!(state, 1)
    snapshot2 = field.getter(nothing, 130.0, integrator)
    @test snapshot1 == [140.0, 160.0]
    @test snapshot2 == [Inf, 160.0]
    @test state.heat_load_switches_s == [(Inf, Inf), (145.0, 150.0)]
    @test state.last_switch_solve_t == [-Inf, 124.0]
    @test state.energy_bracketing_count == [2, 3]
    @test state.last_alpha_rad == [0.1, 1.2]
    table = SimulationModel.IOOutputs._build_results_dataframe([123.0, 130.0],
        [Dict(field.name => snapshot1), Dict(field.name => snapshot2)], [field], args)
    @test names(table) == ["time", "sc1_edg_targeting_switch_s", "sc2_edg_targeting_switch_s"]
    @test table.time == [123.0, 130.0]
    @test table.sc1_edg_targeting_switch_s == [140.0, Inf]
    @test table.sc2_edg_targeting_switch_s == [160.0, 160.0]
    Guidance._edg_initialize_state!(state)
    @test state.targeting_switch_s == [Inf, Inf]
    @test state.energy_bracketing_count == [0, 0]
    @test all(isnan, state.last_alpha_rad)
end

module LegacyLambdaBoundary
scatter(; kwargs...) = (; kwargs...)
attr(; kwargs...) = (; kwargs...)
Layout(; kwargs...) = (; kwargs...)
plot(traces, layout) = (; traces, layout)
const DISPLAY_CALLS = Ref(0)
display(figure) = (DISPLAY_CALLS[] += 1; nothing)
savefig(figure, path; format) = (; figure, path, format)
end
Core.eval(LegacyLambdaBoundary, source_function("scripts/plotting/plot_data.jl", :lambda_plot))

@testset "legacy lambda renderer is a raw two-pair consumer" begin
    # Synthetic consumer input characterizes rendering only. It does not assert
    # that a tracked producer supplies these costate histories.
    input = (time_switch_list=[12.0, 11.0, 11.0], lambda_switch_list=[3.0, -2.0, -0.0],
        time_list=[25.0, 24.0], lamv_list=[-4.0, 5.0])
    result = LegacyLambdaBoundary.lambda_plot("diagnostic", input, nothing)
    threshold, costate = result.figure.traces
    @test threshold.x === input.time_switch_list
    @test threshold.y === input.lambda_switch_list
    @test costate.x === input.time_list
    @test costate.y === input.lamv_list
    @test threshold.name == "λ_switch"
    @test costate.name == "λᵥ"
    @test threshold.line.color == "black"
    @test costate.line.color == "red"
    @test threshold.mode == costate.mode == "lines"
    @test result.figure.layout.xaxis_title == "Time [s]"
    @test result.figure.layout.yaxis_title == "λᵥ"
    @test result.path == "diagnostic_lambda_v.pdf"
    @test result.format == "pdf"
    @test LegacyLambdaBoundary.DISPLAY_CALLS[] == 1

    empty_input = (time_switch_list=Float64[], lambda_switch_list=Float64[],
        time_list=Float64[], lamv_list=Float64[])
    empty_result = LegacyLambdaBoundary.lambda_plot("empty", empty_input, nothing)
    @test length(empty_result.figure.traces) == 2
    @test all(trace -> isempty(trace.x) && isempty(trace.y), empty_result.figure.traces)
    # Even a mismatched x/y pair reaches the backend unchanged. The retained
    # renderer establishes neither data validation nor a valid costate history.
    mismatched = merge(input, (time_list=[1.0], lamv_list=[2.0, 3.0]))
    mismatch_result = LegacyLambdaBoundary.lambda_plot("mismatch", mismatched, nothing)
    @test mismatch_result.figure.traces[2].x == [1.0]
    @test mismatch_result.figure.traces[2].y == [2.0, 3.0]
end
end # module
