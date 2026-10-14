module LegacyEulerContractTests

using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const QM = SM.QuaternionMath
const ENGINE = SpaceAGORA.SimulationEngine
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const QI = SVector(0.0, 0.0, 0.0, 1.0)

# Extract actual definitions, excluding plotting imports and the dispatcher.
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
    return only(definitions)
end

module LegacyBoundary
using ..LegacyEulerContractTests: QM
const qToEulerAngles = QM.qToEulerAngles
const DISPLAYS = Any[]
const SAVES = Any[]
scatter(; kwargs...) = (; kwargs...)
attr(; kwargs...) = (; kwargs...)
Layout(; kwargs...) = (; kwargs...)
plot(traces, layout) = (; traces, layout)
display(figure) = push!(DISPLAYS, figure)
function savefig(figure, path; kwargs...)
    push!(SAVES, (; path, kwargs...))
    return figure
end
end
Core.eval(LegacyBoundary, source_function(
    joinpath(ROOT, "scripts", "plotting", "plot_data.jl"), :attitude_plot))

# Evaluate the real default quaternion SaveField expression and its actual
# getter. No model construction, propagation, or external assets are needed.
function quaternion_field_calls(form)
    form isa Expr || return Expr[]
    if form.head === :call && form.args[1] === :SaveField &&
       any(arg -> arg isa QuoteNode && arg.value === :quaternion, form.args)
        return [form]
    end
    return reduce(vcat, (quaternion_field_calls(arg) for arg in form.args); init=Expr[])
end
module SavedProducer
using ..LegacyEulerContractTests: SM
const SaveField = SM.SaveField
const _save_quaternion = SM.SimulationCallbacks._save_quaternion
const num_sats = 2
end
const RECORDED_FIELD = Core.eval(SavedProducer, only(quaternion_field_calls(source_function(
    joinpath(ROOT, "src", "orchestration", "simulation", "callbacks", "save_fields.jl"), :default_save_fields))))

axis_quaternion(axis, angle) = SVector{4, Float64}(
    ntuple(i -> i == 4 ? cos(angle / 2) : (i == axis ? sin(angle / 2) : 0.0), 4))

@testset "legacy Euler arithmetic is a compatibility contract" begin
    # The old Roll/Pitch/Yaw labels do not imply conventional Euler angles.
    # In particular, atan's arguments are reversed for roll and yaw.
    @test QM.qToEulerAngles(QI) == [pi / 2, 0.0, pi / 2]
    for angle in (-0.7, 0.4), axis in 1:3
        expected = [pi / 2, 0.0, pi / 2]
        expected[axis] = axis == 2 ? -angle : pi / 2 - angle
        @test QM.qToEulerAngles(axis_quaternion(axis, angle)) ≈ expected atol=1e-15
    end

    roll, pitch, yaw = 0.3, -0.4, 0.6
    q = QM.quat_mult(axis_quaternion(3, yaw),
                     QM.quat_mult(axis_quaternion(2, pitch), axis_quaternion(1, roll)))
    # Independent matrix composition pins the scalar-last inertial-to-body
    # convention before characterizing the nonstandard legacy angle output.
    rx = [1.0 0.0 0.0; 0.0 cos(roll) -sin(roll); 0.0 sin(roll) cos(roll)]
    ry = [cos(pitch) 0.0 sin(pitch); 0.0 1.0 0.0; -sin(pitch) 0.0 cos(pitch)]
    rz = [cos(yaw) -sin(yaw) 0.0; sin(yaw) cos(yaw) 0.0; 0.0 0.0 1.0]
    @test QM.rot(q) ≈ transpose(rz * ry * rx) atol=1e-15
    expected = [pi / 2 - roll, -pitch, pi / 2 - yaw]
    @test QM.qToEulerAngles(q) ≈ expected atol=1e-15
    for scale in (-7.0, -1.0, 0.25, 8.0, 1e-300, 1e300)
        raw = collect(scale * q)
        before = copy(raw)
        @test QM.qToEulerAngles(raw) ≈ expected atol=2e-15
        @test raw == before
    end

    # Unlike project_unit_quaternion, normalize has no near-zero identity
    # policy: even a quaternion far below the canonical cutoff retains attitude.
    tiny = [1e-300, 0.0, 0.0, 1e-300]
    @test QM.project_unit_quaternion(tiny) == QI
    @test QM.qToEulerAngles(tiny) ≈ [0.0, 0.0, pi / 2] atol=1e-15
    @test all(isnan, QM.qToEulerAngles(zeros(4)))
    for value in (NaN, Inf, -Inf), component in 1:4
        qbad = collect(QI)
        qbad[component] = value
        @test all(isnan, QM.qToEulerAngles(qbad))
    end
end

@testset "legacy branch and gimbal behavior is not repaired by plotting" begin
    delta = 1e-6
    for axis in (1, 3)
        before = QM.qToEulerAngles(axis_quaternion(axis, -pi / 2 - delta))[axis]
        after = QM.qToEulerAngles(axis_quaternion(axis, -pi / 2 + delta))[axis]
        @test before ≈ -pi + delta atol=1e-15
        @test after ≈ pi - delta atol=1e-15
        @test after - before > 6.0
    end
    # Across pitch lock, roll and yaw can change branches while pitch remains
    # close to -pi/2. There is no chosen unique roll/yaw pair at the lock.
    for sign in (-1.0, 1.0)
        lower = QM.qToEulerAngles(axis_quaternion(2, sign * (pi / 2 - delta)))
        upper = QM.qToEulerAngles(axis_quaternion(2, sign * (pi / 2 + delta)))
        @test lower ≈ [pi / 2, -sign * (pi / 2 - delta), pi / 2] atol=3e-10
        @test upper ≈ [-pi / 2, -sign * (pi / 2 - delta), -pi / 2] atol=3e-10
    end
    # Exact pitch lock can overshoot the sqrt domain after normalization. Norm
    # rounding can vary across platforms, so characterize the actual arithmetic
    # without requiring every platform to produce the same overshoot.
    lock_q = [1.0, 5.0, -1.0, 5.0]
    unit = normalize(lock_q)
    pitch_intermediate = 2 * (unit[4] * unit[2] - unit[1] * unit[3])
    if abs(pitch_intermediate) > 1.0
        @test_throws DomainError QM.qToEulerAngles(lock_q)
    else
        lock_angles = QM.qToEulerAngles(lock_q)
        @test all(isfinite, lock_angles)
        @test lock_angles[2] ≈ -pi / 2 atol=4sqrt(eps(Float64))
    end
end

@testset "actual saved quaternion producer, frame and column order" begin
    @test RECORDED_FIELD.name === :quaternion
    @test RECORDED_FIELD.column_prefix == "q"
    @test RECORDED_FIELD.per_satellite
    times = [120.0, 0.0, 60.0, 60.0]
    quaternions = [collect(QI), [2.0, -3.0, 4.0, -5.0], -collect(QI), zeros(4)]
    z3 = SVector(0.0, 0.0, 0.0)
    states = [(sc=[(q=q, pos_ii=z3, vel_ii=z3, mass_kg=1.0, ω=z3),
                   (q=-q, pos_ii=z3, vel_ii=z3, mass_kg=1.0, ω=z3)],) for q in quaternions]
    before = deepcopy(states)
    snapshots = [Dict(RECORDED_FIELD.name =>
        Base.invokelatest(RECORDED_FIELD.getter, u, t, nothing)) for (u, t) in zip(states, times)]
    args = (dynamics_model=(spacecraft=[(id=41,), (id=97,)],),)
    table = SM.IOOutputs._build_results_dataframe(times, snapshots, [RECORDED_FIELD], args)
    @test propertynames(table) == [:time; [Symbol("sc$(sc)_q_$(i)") for sc in 1:2 for i in 1:4]]
    @test table.time == times
    for spacecraft in 1:2, component in 1:4
        @test isequal(table[!, Symbol("sc$(spacecraft)_q_$(component)")],
            [u.sc[spacecraft].q[component] for u in states])
    end
    for u in states
        sample = ENGINE.build_state_sample(u.sc[1], nothing, true)
        @test sample.q_ib == SVector{4, Float64}(u.sc[1].q)
    end
    @test isequal(states, before)
    @test_throws ArgumentError SM.SimulationCallbacks._save_quaternion(1, (sc=[(ω=z3,)],), 0.0, nothing)

    # The legacy plot consumes only its unprefixed schema. The explicit adapter
    # is test-local; no default reader silently changes current saved tables.
    @test_throws ArgumentError Base.invokelatest(LegacyBoundary.attitude_plot, "unused", nothing, table)
    legacy_table = (; time=table.time, (Symbol("q_$(i)") => table[!, Symbol("sc1_q_$(i)")] for i in 1:4)...)
    figure = Base.invokelatest(LegacyBoundary.attitude_plot, "contract", nothing, legacy_table)
    @test getproperty.(figure.traces, :name) == ["Roll", "Pitch", "Yaw"]
    @test [trace.line.color for trace in figure.traces] == ["red", "green", "blue"]
    @test figure.layout.xaxis_title == "Time [min]"
    @test figure.layout.yaxis_title == "Angle (Degrees)"
    for component in 1:3
        trace = figure.traces[component]
        @test trace.x == [2.0, 0.0, 1.0, 1.0]
        @test isequal(trace.y, [rad2deg(QM.qToEulerAngles(q)[component]) for q in quaternions])
    end
    @test [trace.y[1] for trace in figure.traces] == [90.0, 0.0, 90.0]
    @test all(isnan, [trace.y[end] for trace in figure.traces])
    @test last(LegacyBoundary.DISPLAYS) === figure
    @test last(LegacyBoundary.SAVES) == (path="contract_attitude.pdf", format="pdf")
    @test isequal(states, before)

    branch_qs = [axis_quaternion(1, -pi / 2 - 1e-6), axis_quaternion(1, -pi / 2 + 1e-6)]
    branch_table = (; time=[0.0, 1.0], (Symbol("q_$(i)") => [q[i] for q in branch_qs] for i in 1:4)...)
    branch_figure = Base.invokelatest(LegacyBoundary.attitude_plot, "branch", nothing, branch_table)
    @test branch_figure.traces[1].y[2] - branch_figure.traces[1].y[1] > 359.0

    singular_table = (time=[0.0], q_1=[1.0], q_2=[5.0], q_3=[-1.0], q_4=[5.0])
    unit = normalize([1.0, 5.0, -1.0, 5.0])
    pitch_intermediate = 2 * (unit[4] * unit[2] - unit[1] * unit[3])
    if abs(pitch_intermediate) > 1.0
        @test_throws DomainError Base.invokelatest(LegacyBoundary.attitude_plot, "singular", nothing, singular_table)
    else
        singular = Base.invokelatest(LegacyBoundary.attitude_plot, "singular", nothing, singular_table)
        @test all(trace -> all(isfinite, trace.y), singular.traces)
        @test only(singular.traces[2].y) ≈ -90.0 atol=rad2deg(4sqrt(eps(Float64)))
    end
end

end # module LegacyEulerContractTests
