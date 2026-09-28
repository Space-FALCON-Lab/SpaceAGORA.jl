# Contract probe for the actual legacy callers. This is not a physical predictor
# validation: only expensive prediction dependencies use deterministic fixtures.
# Run with --depwarn=error and a real Roots installation. Supply the unchanged
# source/candidate repository root as ARGS[1] or SPACEAGORA_LEGACY_SOURCE_ROOT.
module LegacyRootsCallerProbes
using Test, Roots, LinearAlgebra, SHA

const SOURCE_ROOT = abspath(isempty(ARGS) ?
    get(ENV, "SPACEAGORA_LEGACY_SOURCE_ROOT", joinpath(@__DIR__, "..", "..")) : ARGS[1])
const CALLER_FILES = [
    "src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl",
    "src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl",
    "src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl",
]
const SOURCE_HASHES = Dict(path => bytes2hex(sha256(read(joinpath(SOURCE_ROOT, path))))
    for path in CALLER_FILES)

# Include entire actual caller files unchanged. Undefined dependencies of unused
# functions are harmless. Fixture-specific methods below replace predictors,
# never Roots or the four caller functions being tested.
for path in CALLER_FILES
    include(joinpath(SOURCE_ROOT, path))
end

struct PredictorFailure <: Exception end
mutable struct FixtureTrace
    residual_inputs::Vector{Float64}
    predictor_evals::Int
    post_roots::Vector{Float64}
end
FixtureTrace() = FixtureTrace(Float64[], 0, Float64[])
struct FixtureMission
    planet::NamedTuple
    aerodynamics::NamedTuple
    body::Nothing
    kind::Symbol
    mode::Symbol
    root::Float64
    trace::FixtureTrace
end
function mission(kind; mode=:root, root=110.0, heat_limit=100.0)
    FixtureMission((T=1.0, R=1.0, γ=1.4, μ=1.0, Rp_e=1.0, name="earth"),
        (α=1.0, thermal_accomodation_factor=1.0, heat_load_limit=heat_limit),
        nothing, kind, mode, Float64(root), FixtureTrace())
end
const SAVED = (time_switch_1=5.0, time_switch_2=50.0, α_list=[1.0], heat_load_past=[0.0])
args(; heat_load_sol=1, verbose=false) = Dict{Symbol,Any}(
    :heat_load_sol => heat_load_sol, :multiplicative_factor_heatload => 1.0,
    :control_mode => 0, :max_heat_rate => 100.0, :verbose => verbose)
_bridge_get_cnf(::AbstractDict; cnf=nothing) = cnf
_bridge_verbose_enabled(a::AbstractDict) = a[:verbose]

function _control_asim_ctrl(ip, m::FixtureMission, t, current_position, a, k,
        heat_rate_control, initial_guess, gram_atmosphere, ts, reevaluation_mode; cnf=nothing)
    push!(m.trace.residual_inputs, ts)
    n = length(m.trace.residual_inputs)
    if m.mode === :throw_precheck || (m.mode === :throw_root && n > 2)
        throw(PredictorFailure())
    end
    residual = m.mode === :same_sign ? 2.0 + ts / 1000.0 : m.root - ts
    return fill(m.aerodynamics.heat_load_limit + residual, 1, 1)
end

function closed_form(a, m::FixtureMission, position, T, calculate, alpha, profile=nothing)
    if profile === nothing
        m.trace.predictor_evals += 1
        if m.mode === :throw_precheck || (m.mode === :throw_root &&
                m.kind === :second && m.trace.predictor_evals > 2)
            throw(PredictorFailure())
        end
    end
    times = m.kind === :second ? collect(20.0:20.0:220.0) : [0.0, 40.0, 80.0, 160.0]
    return times, ones(length(times)), zeros(length(times)), ones(length(times))
end
density_exp(heights, planet::NamedTuple) = (ones(length(heights)),)
heat_rate_calc(factor, density, T, Tw, R, gamma, speed, angles) = Float64.(angles .> 0.0)
aerodynamic_coefficient_fM(alpha, ::Nothing, T, speed, aero, flag) = (0.0, 1.0 + alpha)
aoa(m::FixtureMission, k, times, heights, gamma, velocity, coeff, nu) = (ones(length(times)), nothing)
lambdas(m::FixtureMission, profile, k, times, heights, gamma, velocity, coeff, nu) =
    (0.0, [1.0, -1.0, -1.0, 1.0], nothing)

function func(k, m::FixtureMission, a, coeff, position, heat_rate_control, approx_sol,
        profile, initial_guess=false, approx_calc=false)
    if approx_calc
        push!(m.trace.post_roots, k)
        times, heights, gamma, velocity = approx_sol
        return times, velocity, gamma, heights
    end
    push!(m.trace.residual_inputs, k)
    n = length(m.trace.residual_inputs)
    m.mode === :throw_root && n > 2 && throw(PredictorFailure())
    m.mode === :below && return -1.0
    m.mode === :above && return 1.0
    # Preserve opposite-sign prechecks, then exercise an invalid solver bracket.
    m.mode === :invalid_after_precheck && n > 2 && return 1.0
    return k - m.root
end

# More-specific dispatch replaces the expensive iterative func_e predictor in
# the included file. The actual targeting caller and its known aoa_cf bug stay.
function func_e(nu, m::FixtureMission, a, coeff, position, heat_rate_control,
        approx_sol, energy_target, initial_guess=false, approx_calc=false)
    if approx_calc
        push!(m.trace.post_roots, nu)
        times, heights, gamma, velocity = approx_sol
        return times, velocity, gamma, heights
    end
    push!(m.trace.residual_inputs, nu)
    m.mode === :throw_root && length(m.trace.residual_inputs) > 2 && throw(PredictorFailure())
    return m.mode === :same_sign ? nu + 1.0 : nu - m.root
end

integrated(m; heat_load_sol=1) = second_time_switch_recalc_with_integration(
    nothing, m, nothing, args(; heat_load_sol), 10.0, false, 0; cnf=SAVED)
sampled(m) = second_time_switch_recalc(
    nothing, m, nothing, args(), 10.0, false; cnf=SAVED)
window(m) = switch_calculation(nothing, m, nothing, args(), 0.0, false, 0)
targeting(m) = control_solarpanels_targeting_closed_form(
    0.0, nothing, m, nothing, args(), 0.0, false, 0)

struct FixtureParams
    mission::FixtureMission
    args::Dict{Symbol,Any}
    ip::Nothing
    time_0::Float64
    gram_atmosphere::Nothing
    cnf::NamedTuple
end
_bridge_get_cnf(p::FixtureParams) = p.cnf
function energy_state(x, m::FixtureMission)
    push!(m.trace.residual_inputs, x)
    # Radius=1 and μ=1 give E=1e6*(1000+x), before the caller's /1e6 scaling.
    return reshape([1.0, 0.0, 0.0, sqrt(2.0 * (1e6 * (1000.0 + x) + 1.0)), 0.0, 0.0], 6, 1)
end
asim_ctrl_targeting(ts, p::FixtureParams, t0, initial; cnf=nothing) = energy_state(ts, p.mission)
_control_asim_ctrl_rf(ip, m::FixtureMission, t0, OE, a, nu, scale, flag, gram; cnf=nothing) =
    (energy_state(nu, m), nothing)
struct CapturedDisplay <: AbstractDisplay
    values::Vector{Any}
end
Base.display(d::CapturedDisplay, value) = (push!(d.values, value); nothing)

@testset "Legacy Roots caller migration contracts" begin
    @test Base.JLOptions().depwarn == 2 # required: deprecated fzero cannot hide in a catch
    expected_version = get(ENV, "SPACEAGORA_EXPECT_ROOTS", "")
    if !isempty(expected_version)
        @test Base.pkgversion(Roots) == VersionNumber(expected_version)
    end
    println("Roots version: ", Base.pkgversion(Roots), "; source root: ", SOURCE_ROOT)
    for path in CALLER_FILES
        println("caller_sha256 ", SOURCE_HASHES[path], " ", path)
    end

    @testset "Integrated second switch: interior roots and unchanged fallback" begin
        for (mode, upper, root) in ((1, 210.0, 110.0), (3, 1510.0, 710.0))
            m = mission(:integrated; root)
            result = integrated(m; heat_load_sol=mode)
            @test result[1] == SAVED.time_switch_1
            @test isapprox(result[2], root; atol=1e-10, rtol=0)
            @test result[2] != SAVED.time_switch_2
            @test m.trace.residual_inputs[1:2] == [50.0, 10.0]
            @test upper in m.trace.residual_inputs[3:end]
            @test any(x -> 10.0 < x < upper, m.trace.residual_inputs[3:end])
            @test abs(m.root - result[2]) <= 1e-10
        end
        for mode in (:same_sign, :throw_root)
            m = mission(:integrated; mode)
            @test integrated(m) == (5.0, 50.0)
            @test m.trace.residual_inputs[1:2] == [50.0, 10.0]
            @test length(m.trace.residual_inputs) > 2
        end
        @test_throws PredictorFailure integrated(mission(:integrated; mode=:throw_precheck))
    end

    @testset "Sampled second switch: zero plateau and 10 percent adjustment" begin
        m = mission(:second)
        result = sampled(m)
        root = result[2] / 0.9
        # Masks make this synthetic residual piecewise constant. Its exact zero
        # plateau is [100,120), not a uniquely defined switch at 110 seconds.
        @test result[1] == 5.0
        @test 100.0 <= root < 120.0
        @test result[2] != 0.9 * SAVED.time_switch_2
        remaining = filter(t -> t > root, collect(20.0:20.0:220.0))
        @test last(remaining) - first(remaining) - m.aerodynamics.heat_load_limit == 0.0
        @test m.trace.predictor_evals > 2
        no_bracket = mission(:second; heat_limit=1000.0)
        @test sampled(no_bracket) == (5.0, 45.0)
        @test no_bracket.trace.predictor_evals > 2
        failing = mission(:second; mode=:throw_root)
        @test sampled(failing) == (5.0, 45.0)
        @test failing.trace.predictor_evals == 3 # two prechecks succeeded
        @test_throws PredictorFailure sampled(mission(:second; mode=:throw_precheck))
    end

    @testset "Switch window: interior root, endpoint policies, uncaught failures" begin
        m = mission(:window; root=0.04)
        @test window(m) == [40.0, 72.0]
        @test m.trace.residual_inputs[1:2] == [0.1, 0.0]
        @test length(m.trace.post_roots) == 1
        if length(m.trace.post_roots) == 1
            @test isapprox(only(m.trace.post_roots), 0.04; atol=1e-12, rtol=0)
        end
        @test any(x -> 0.0 < x < 0.1, m.trace.residual_inputs[3:end])
        @test window(mission(:window; mode=:below)) == [0.0, 0.0]
        @test window(mission(:window; mode=:above)) == [0.0, 80.0]
        for (mode, error_type) in ((:invalid_after_precheck, ArgumentError), (:throw_root, PredictorFailure))
            m = mission(:window; mode, root=0.04)
            @test_throws error_type window(m)
            @test m.trace.residual_inputs[1:2] == [0.1, 0.0]
            @test length(m.trace.residual_inputs) > 2
            @test isempty(m.trace.post_roots)
        end
    end

    @testset "Closed-form targeting: root reached, known post-root defect explicit" begin
        m = mission(:targeting; root=37.0)
        err = try
            targeting(m)
            nothing
        catch caught
            caught
        end
        @test err isa UndefVarError
        if err isa UndefVarError
            @test err.var === :aoa_cf
        end
        @test any(x -> 1.0 < x < 100.0, m.trace.residual_inputs)
        @test any(x -> isapprox(x, 37.0; atol=1e-10, rtol=0), m.trace.residual_inputs)
        @test length(m.trace.post_roots) == 1 # caller advanced beyond Roots
        @test !isdefined(@__MODULE__, :aoa_cf) # do not mask the existing bug
        @test_throws ArgumentError targeting(mission(:targeting; mode=:same_sign))
        @test_throws PredictorFailure targeting(mission(:targeting; mode=:throw_root, root=37.0))
    end

    @testset "Two targeting trace displays are explicitly controlled" begin
        for (kind, root) in ((:num_int, 250.0), (:heatload, 333.0)), enabled in (false, true)
            m = mission(:logging; root)
            p = FixtureParams(m, args(; verbose=enabled), nothing, 0.0, nothing, SAVED)
            target_energy = 1e6 * (1000.0 + root)
            output = CapturedDisplay(Any[])
            pushdisplay(output)
            recovered = try
                kind === :num_int ? control_solarpanels_targeting_num_int(target_energy, p, 0.0, nothing) :
                    control_solarpanels_targeting_heatload(target_energy, p, nothing)
            finally
                popdisplay(output)
            end
            @test isapprox(recovered, root; atol=1e-6, rtol=0)
            @test length(m.trace.residual_inputs) > 2
            @test length(output.values) == Int(enabled)
            if enabled && length(output.values) == 1
                @test only(output.values) isa Roots.Tracks
                @test occursin("Converged", sprint(show, only(output.values)))
            end
        end
    end

    @testset "Included caller files remain unchanged" begin
        for path in CALLER_FILES
            @test bytes2hex(sha256(read(joinpath(SOURCE_ROOT, path)))) == SOURCE_HASHES[path]
        end
    end
end
println("legacy_roots_callers_ok; known aoa_cf post-root failure remains unresolved")
end
