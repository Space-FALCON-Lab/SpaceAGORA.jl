module LegacyWheelContractTests

using Test
using LinearAlgebra
using StaticArrays
using DataFrames
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SC = SM.SimulationCallbacks
const SE = SpaceAGORA.SimulationEngine
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Evaluate only named source boundaries. Importing the demo would load telemetry
# and launch a campaign; rendering the legacy functions would open plot windows.
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

function save_field_call(form, field_name)
    form isa Expr || return Expr[]
    if form.head === :call && form.args[1] == :(SC.SaveField) &&
       any(arg -> arg isa QuoteNode && arg.value === field_name, form.args)
        return [form]
    end
    return reduce(vcat, (save_field_call(arg, field_name) for arg in form.args); init=Expr[])
end

module PrescribedBoundary
using LinearAlgebra
using StaticArrays
using ..LegacyWheelContractTests: SM, SC
const AXES = [1.0 0.0 0.0 0.5; 0.0 2.0 0.0 -0.25; 0.0 0.0 1.0 0.75]
const SPEED0 = [4.0, -3.0, 2.0, -1.0]
const RATE = [0.5, 1.0, -0.25, 2.0]
const WHEEL_I = 0.125
const OFFSET = 0.25
const KNOTS = [0.0, 1.0, 2.0]
const wheels = SM.ReactionWheelMomentumModel(KNOTS,
    [SPEED0[k] + RATE[k] * t for t in KNOTS, k in 1:4], AXES, WHEEL_I;
    time_offset_s=OFFSET)
ω_of(u) = SVector{3, Float64}(u.sc[1].ω)
end

const SLEW = source_function(joinpath(ROOT, "scripts", "dev", "viewer_demos", "cygnss_slew.jl"),
                             :slew_save_fields)
const MOMENTUM_FIELD = Core.eval(PrescribedBoundary, only(save_field_call(SLEW, :wheel_momentum_nms)))
const TORQUE_MAGNITUDE_FIELD = Core.eval(PrescribedBoundary, only(save_field_call(SLEW, :wheel_torque_nm)))

module LegacyBoundary
scatter(; kwargs...) = (; kwargs...)
attr(; kwargs...) = (; kwargs...)
Layout(; kwargs...) = (; kwargs...)
plot(traces, layout) = (; traces, layout)
display(figure) = nothing
savefig(figure, args...; kwargs...) = figure
end
for name in (:reaction_wheel_h_plot, :reaction_wheel_torque_plot,
             :total_reaction_wheel_torque_plot, :torque_plot, :inertia_plot)
    Core.eval(LegacyBoundary,
        source_function(joinpath(ROOT, "scripts", "plotting", "plot_data.jl"), name))
end

@testset "commanded root wheel momenta and actual saved output boundary" begin
    # Four wheel coordinates cannot be mistaken for three body components.
    axes = PrescribedBoundary.AXES
    h = [0.25, -0.125, 0.5, -0.375]
    wheel = SM.ReactionWheelAssembly(n_wheels=4, J_rw=axes, h_wheels=copy(h))
    stale_hdot = [11.0, -13.0, 17.0, -19.0]
    stale_torque = SVector(23.0, -29.0, 31.0)
    wheel.h_dot_wheels .= stale_hdot
    wheel.tau_body_net .= stale_torque
    omega = SVector(0.01, -0.02, 0.03)
    inertia = Diagonal([2.0, 3.0, 5.0])
    body_torque = SVector(0.125, -0.25, 0.5)
    state = (q=[0.0, 0.0, 0.0, 1.0], ω=collect(omega), h_wheels=copy(h))
    derivative = (q=zeros(4), ω=zeros(3), h_wheels=zeros(4))
    SE._assign_orientation_rhs!(derivative, state, inertia, body_torque;
        propagate_quaternion=false, include_gyroscopic=true,
        rw_assembly=wheel, rw_torque_body=body_torque)
    @test derivative.h_wheels ≈ pinv(axes) * (-body_torque)
    @test axes * derivative.h_wheels ≈ -body_torque
    @test derivative.ω ≈ inertia \ (body_torque - cross(omega, inertia * omega + axes * h))
    @test state.h_wheels == h
    # These assembly diagnostics are not populated by this producer.
    @test wheel.h_dot_wheels == stale_hdot
    @test wheel.tau_body_net == stale_torque

    # Explicit recipe: copy integrated coordinates in root assembly column order.
    # Different spacecraft can have different wheel counts; N is list position.
    field = SM.SaveField(:wheel_h_nms,
        (u, t, integrator) -> [copy(u.sc[i].h_wheels)
            for i in eachindex(integrator.p.args.dynamics_model.spacecraft)];
        per_satellite=true, column_prefix="wheel_h_nms")
    times = [60.0, 0.0, 60.0]
    states = [(sc=[(h_wheels=copy(h .+ t),), (h_wheels=[-t, t + 0.5],)],) for t in times]
    args = (dynamics_model=(spacecraft=[(id=41,), (id=97,)],),)
    integrator = (p=(args=args,),)
    snapshots = [SC._save_snapshot([field], u, t, integrator) for (u, t) in zip(states, times)]
    states[1].sc[1].h_wheels .= 999.0
    table = SM.IOOutputs._build_results_dataframe(times, snapshots, [field], args)
    @test names(table) == ["time"; ["sc1_wheel_h_nms_$(i)" for i in 1:4];
                                    ["sc2_wheel_h_nms_$(i)" for i in 1:2]]
    @test table.time == times
    for i in 1:4
        @test table[!, "sc1_wheel_h_nms_$(i)"] == h[i] .+ times
    end
    @test table.sc2_wheel_h_nms_1 == -times
    @test !any(contains("sc41"), names(table))
    @test !any(contains("torque"), names(table))
    @test !any(contains("inertia"), names(table))

    # Legacy nested history is wheel-major, not a body vector.
    histories = [table[!, "sc1_wheel_h_nms_$(i)"] for i in 1:4]
    solution = (orientation=(time=times,), physical_properties=(rw_h=histories,))
    figure = LegacyBoundary.reaction_wheel_h_plot("unused", solution, table)
    @test [trace.y for trace in figure.traces] == histories
    @test all(trace -> trace.x == times, figure.traces)
    @test figure.layout.yaxis_title == "Angular Momentum (kg*m²/s)"
    # Retained legacy renderer limitation: i % 9 selects Julia index zero at
    # wheel 9. This is separate from the saver, which has no nine-wheel limit.
    nine_histories = [fill(Float64(i), length(times)) for i in 1:9]
    nine_solution = (orientation=(time=times,), physical_properties=(rw_h=nine_histories,))
    @test_throws BoundsError LegacyBoundary.reaction_wheel_h_plot("unused", nine_solution, table)
    nine_state = (sc=[(h_wheels=collect(1.0:9.0),)],)
    nine_args = (dynamics_model=(spacecraft=[(id=41,)],),)
    nine_snapshot = SC._save_snapshot([field], nine_state, 0.0, (p=(args=nine_args,),))
    nine_table = SM.IOOutputs._build_results_dataframe([0.0], [nine_snapshot], [field], nine_args)
    @test nine_table.sc1_wheel_h_nms_9 == [9.0]
end

@testset "prescribed wheel saver records aggregate body momentum and torque norm" begin
    P = PrescribedBoundary
    times = [0.75, 0.0, 0.75]
    omega = SVector(0.1, -0.2, 0.3)
    state = (sc=[(ω=omega,)],)
    fields = [MOMENTUM_FIELD, TORQUE_MAGNITUDE_FIELD]
    # Last-stage diagnostic state is deliberately stale. The real demo getters
    # query the prescribed model at each requested save time instead.
    fill!(P.wheels.state.momentum_body, SVector(999.0, 999.0, 999.0))
    fill!(P.wheels.state.torque_body, SVector(888.0, 888.0, 888.0))
    snapshots = [SC._save_snapshot(fields, state, t, nothing) for t in times]
    args = (dynamics_model=(spacecraft=[(id=41,)],),)
    table = SM.IOOutputs._build_results_dataframe(times, snapshots, fields, args)
    @test names(table) == ["time", "sc1_wheel_momentum_nms_1", "sc1_wheel_momentum_nms_2",
                          "sc1_wheel_momentum_nms_3", "sc1_wheel_torque_nm"]
    for (row, t) in enumerate(times)
        speed = P.SPEED0 + P.RATE * (t + P.OFFSET)
        momentum = P.AXES * (P.WHEEL_I * speed)
        momentum_rate = P.AXES * (P.WHEEL_I * P.RATE)
        torque = -momentum_rate - cross(omega, momentum)
        @test SM.wheel_speeds_rad_s(P.wheels, t) ≈ speed
        @test SM.wheel_momentum_rate_body(P.wheels, t) ≈ momentum_rate
        @test SM.wheel_reaction_torque(SVector{3}(momentum), SVector{3}(momentum_rate), omega) ≈ torque
        @test [table[row, "sc1_wheel_momentum_nms_$(k)"] for k in 1:3] ≈ momentum
        @test table.sc1_wheel_torque_nm[row] ≈ norm(torque)
        @test length(speed) == 4
        @test length(momentum) == 3
    end
    @test all(==(SVector(999.0, 999.0, 999.0)), P.wheels.state.momentum_body)
    @test all(==(SVector(888.0, 888.0, 888.0)), P.wheels.state.torque_body)
end

@testset "legacy torque and inertia require their distinct input histories" begin
    table = DataFrame(time=[2.0, 0.0, 2.0])
    for (prefix, count) in (("rw_tau", 4), ("rw_torque_ii", 3), ("tau_ii", 3), ("J_ii", 9))
        for i in 1:count
            table[!, "$(prefix)_$(i)"] = [-i, i + 0.5, -i - 0.25]
        end
    end
    per_wheel = LegacyBoundary.reaction_wheel_torque_plot("unused", (body=(n_reaction_wheels=4,),), table)
    reaction = LegacyBoundary.total_reaction_wheel_torque_plot("unused", table)
    total = LegacyBoundary.torque_plot("unused", table)
    inertia = LegacyBoundary.inertia_plot("unused", table)
    for (figure, prefix, count) in ((per_wheel, "rw_tau", 4), (reaction, "rw_torque_ii", 3),
                                    (total, "tau_ii", 3), (inertia, "J_ii", 9))
        @test length(figure.traces) == count
        for i in 1:count
            @test figure.traces[i].x == table.time
            @test figure.traces[i].y == table[!, "$(prefix)_$(i)"]
        end
    end
    @test inertia.layout.yaxis_title == "Inertia Tensor (kg⋅m²)"
    for i in 5:9
        table[!, "rw_tau_$(i)"] = fill(Float64(i), nrow(table))
    end
    @test_throws BoundsError LegacyBoundary.reaction_wheel_torque_plot(
        "unused", (body=(n_reaction_wheels=9,),), table)
end

end # module
