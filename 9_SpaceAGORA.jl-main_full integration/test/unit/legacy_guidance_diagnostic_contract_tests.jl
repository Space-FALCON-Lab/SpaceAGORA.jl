module LegacyGuidanceDiagnosticContractTests
using Test
using LinearAlgebra
using SpaceAGORA: SimulationModel

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const Cnf = SimulationModel.ConfigTypes.Cnf
const EDG = SimulationModel.EDGAlgorithms
const Guidance = SimulationModel.GuidanceHooks

# Characterize the retained storage boundary without running a shooting solve,
# loading SPICE kernels, or claiming a physical qualification of its trajectory.
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

const TARGETING = source_function("src/models/gnc/guidance/aerobraking/t_edg/eom_predictor.jl", :asim_ctrl_targeting_plot)
const TARGETING_BODY = filter(x -> !(x isa LineNumberNode), TARGETING.args[2].args)
const TAIL_START = only(findall(x -> x isa Expr && x.head === :(=) && x.args[1] === :v_ii_mag, TARGETING_BODY))
# Execute the original speed extraction, switching formula and two append calls.
# The final return is also retained; the solver result is a synthetic boundary input.
Core.eval(@__MODULE__, Expr(:function,
    :(record_targeting!(cnf_state, sol, m, CD_slope)),
    Expr(:block, TARGETING_BODY[TAIL_START:end]...)))
const TIME_INITIALIZATION = only(filter(TARGETING_BODY) do x
    x isa Expr && x.head === :(=) && x.args[1] == :(initial_time, final_time)
end)
Core.eval(@__MODULE__, Expr(:function, :(targeting_time_span(time_0)),
    Expr(:block, TIME_INITIALIZATION, :(return (initial_time, final_time)))))

const CONSTRAINT = source_function("src/models/gnc/control/aerobraking/constraint_tracking.jl", :asim_ctrl_plot)
const CONSTRAINT_BODY = filter(x -> !(x isa LineNumberNode), CONSTRAINT.args[2].args)
function original_push(field)
    return only(filter(CONSTRAINT_BODY) do x
        x isa Expr && x.head === :call && x.args[1] === :push! &&
            x.args[2] == Expr(:., :cnf_state, QuoteNode(field))
    end)
end
Core.eval(@__MODULE__, Expr(:function, :(constraint_value_append!(cnf_state, lambda_switch_list)),
    Expr(:block, original_push(:lambda_switch_list))))
Core.eval(@__MODULE__, Expr(:function, :(constraint_time_append!(cnf_state, sol)),
    Expr(:block, original_push(:time_list2))))

struct SyntheticSolution
    t::Vector{Float64}
    values::Matrix{Float64}
end
Base.getindex(sol::SyntheticSolution, rows, column) = sol.values[rows, column]
function solver_samples(times, velocities)
    values = zeros(10, length(times))
    for (i, velocity) in enumerate(velocities)
        values[4:6, i] .= velocity
    end
    return SyntheticSolution(Float64.(times), values)
end

@testset "legacy runtime arrays and T-EDG recording boundary" begin
    first_runtime, second_runtime = Cnf(), Cnf()
    for field in (:lambda_switch_list, :time_switch_list, :time_list, :lamv_list)
        @test fieldtype(Cnf, field) === Vector{Float64}
        @test isempty(getfield(first_runtime, field))
        @test getfield(first_runtime, field) !== getfield(second_runtime, field)
    end
    @test !hasproperty(first_runtime, :time_list2)
    @test targeting_time_span(123.0) == (123.0, 1623.0)

    # Repeated/nonmonotone solver times deliberately demonstrate no sorting or
    # local-origin subtraction in this recording boundary, not solver behavior.
    sol = solver_samples([123, 125, 125, 124],
        [(3.0, -4.0, 0.0), (-5.0, 0.0, 12.0), (0.0, 0.0, 0.0), (0.0, 8.0, -15.0)])
    model = (body=(mass=3.0, area_tot=2.0),)
    @test record_targeting!(first_runtime, sol, model, 2.0) === sol
    expected = 1.5 .* [5.0, 13.0, 0.0, 17.0] ./ pi
    @test first_runtime.lambda_switch_list ≈ expected
    @test first_runtime.time_switch_list == sol.t
    @test first_runtime.time_switch_list !== sol.t
    @test isempty(first_runtime.time_list)
    @test isempty(first_runtime.lamv_list)

    # Same caller-owned Cnf accumulates another invocation without segmentation.
    record_targeting!(first_runtime, sol, model, -2.0)
    @test first_runtime.lambda_switch_list ≈ vcat(expected, -expected)
    @test first_runtime.time_switch_list == vcat(sol.t, sol.t)
    @test isempty(second_runtime.lambda_switch_list)
    @test isempty(second_runtime.time_switch_list)

    # Record actual defect behavior rather than silently fixing the producer.
    @test_throws MethodError constraint_value_append!(second_runtime, expected)
    @test isempty(second_runtime.lambda_switch_list)
    @test_throws ErrorException constraint_time_append!(second_runtime, sol)
end

@testset "typed costates remain separate call-local prediction arrays" begin
    track = (time=[0.0, 2.0], speed=[4.0, 5.0], h=[0.0, 0.0],
        gamma=[-0.2, 0.0], rho=[0.0, 0.0])
    seed = [0.4, 0.7]
    coeffs = (cd_slope=2.0, cd_low=1.0, cl_low=0.0)
    planet = (μ=100.0, Rp_e=10.0, g_ref=1.0)
    before = deepcopy(track)
    threshold, costate = EDG._edg_heat_load_lambdas(track, seed, coeffs, 3.0, 2.0, planet, 10.0, 4.0)
    @test threshold ≈ 6.0 .* track.speed ./ pi
    @test costate ≈ [4.6, 5.0]
    @test track == before
    @test seed == [0.4, 0.7]
    other_threshold, other_costate = EDG._edg_heat_load_lambdas(track, seed, coeffs, 3.0, 2.0, planet, 10.0, 4.0)
    @test threshold !== other_threshold
    @test costate !== other_costate
    negative_slope, _ = EDG._edg_heat_load_lambdas(track, seed,
        merge(coeffs, (cd_slope=-2.0,)), 3.0, 2.0, planet, 10.0, 4.0)
    @test negative_slope == threshold # Typed denominator uses abs; legacy does not.
    negative_gain, _ = EDG._edg_heat_load_lambdas(track, seed, coeffs, 3.0, 2.0, planet, 10.0, -4.0)
    @test negative_gain == -threshold

    # Exact equality selects the maximum angle in this algorithm.
    one_point = (time=[0.0], speed=[5.0], h=[0.0], gamma=[0.0], rho=[0.0])
    config = (min_alpha_rad=0.1, max_alpha_rad=1.2)
    profile, boundary, lambda = EDG._edg_heat_load_alpha_profile(config, one_point,
        (cd_slope=2/pi, cd_low=1.0, cl_low=0.0), 1.0, 1.0, planet, 10.0, 1.0, [0.7])
    @test boundary == lambda == [5.0]
    @test profile == [config.max_alpha_rad]
    profile_low, _, _ = EDG._edg_heat_load_alpha_profile(config, one_point,
        (cd_slope=2/pi, cd_low=1.0, cl_low=0.0), 1.0, 1.0, planet, 10.0, 2.0, [0.7])
    @test profile_low == [config.min_alpha_rad]

    # Relative prediction samples become elapsed simulation switch timestamps.
    time_track = (time=[0.0, 2.0, 7.0, 10.0],)
    @test EDG._edg_low_alpha_switch_window(config, 123.0, time_track, [1.2, 0.1, 0.1, 1.2]) == (125.0, 130.0)
    @test EDG._edg_low_alpha_switch_window(config, 123.0, time_track, fill(1.2, 4)) == (Inf, Inf)
    @test EDG._edg_low_alpha_switch_window(config, 123.0, time_track, [1.2, 0.1, 1.2, 1.2]) == (Inf, Inf)
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
