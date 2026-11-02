# Work-precision of the J2-only KS splitting propagator against Tsit5 and DP8.
#
# Earth, 300 x 3000 km altitude, i = 51.6 deg, starting at perigee, ten orbits.
# The error is the maximum Cartesian position error over the trajectory against
# a Cartesian Vern9 reference: at every step for the splitting runs, and on a
# dense grid for the Runge-Kutta runs. Timings and evaluation counts come from
# separate runs without saved output.
#
# Usage:
#   julia --project=. --threads=1 benchmarks/studies/ks_splitting_work_precision/work_precision.jl

using CSV
using DataFrames
using Dates
using LinearAlgebra
using OrdinaryDiffEq
using OrdinaryDiffEqHighOrderRK: DP8
using SpaceAGORA
using SpaceAGORA.SimulationModel.DynamicsKS
using StaticArrays
using Statistics

const PLANET = SpaceAGORA.SimulationModel.Earth()
const PARAMS = KSPropagationParams(Re=PLANET.Rp_e, μ=PLANET.μ, J2=PLANET.J2)
const ORBIT_COUNT = 10
const ACCURACY_SAMPLES = 20_001
const REFERENCE_TOLERANCE = 1.0e-13
const SPLITTING_SWEEPS = (
    2 => [16, 32, 64, 128, 256, 512, 1024, 2048, 4096],
    4 => [8, 16, 32, 64, 128, 256, 512, 1024],
)
const RK_TOLERANCES = [1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10, 1.0e-11, 1.0e-12]
const RESULTS_DIR = joinpath(@__DIR__, "results")

function cartesian_j2_rhs(x, params, _)
    r = SVector(x[1], x[2], x[3])
    a = (-params.μ / norm(r)^3) * r + ks_j2_acceleration_si(r, params)
    return SVector(x[4], x[5], x[6], a[1], a[2], a[3])
end

ks_j2_rhs(x, params, _) = nonlinear_ks_rhs(x, params)

function maximum_position_error(times, positions, reference)
    return maximum(zip(times, positions)) do (t, position)
        norm(position - SVector{3}(view(reference(t), 1:3)))
    end
end

# Median after one warm-up call, which excludes compilation.
function median_runtime(operation; target_total_s=0.5)
    first_time = @elapsed operation()
    runs = clamp(ceil(Int, target_total_s / max(first_time, 1.0e-9)), 5, 500)
    samples = [@elapsed(operation()) for _ in 1:runs]
    return median(samples), std(samples)
end

row(method, knob_name, knob_value, runtime, steps, evaluations, error) = (
    method=method, knob_name=knob_name, knob_value=Float64(knob_value),
    runtime_median_s=runtime[1], runtime_std_s=runtime[2], steps=steps,
    force_evaluations=evaluations, maximum_position_error_m=error)

function splitting_rows!(rows, initial, final_time_s, reference)
    for (order, sweep) in SPLITTING_SWEEPS, steps_per_orbit in sweep
        runtime = median_runtime(() -> ks_splitting_propagate(initial, PARAMS,
            final_time_s, steps_per_orbit; order, save_history=false))
        states, work = ks_splitting_propagate(initial, PARAMS, final_time_s,
            steps_per_orbit; order)
        cartesian = ks_state_to_cartesian.(states)
        error = maximum_position_error(getproperty.(cartesian, :elapsed_time_s),
            getproperty.(cartesian, :position_ii_m), reference)
        push!(rows, row("KS splitting order $order", "steps_per_orbit",
            steps_per_orbit, runtime, work.steps, work.evaluations, error))
    end
end

# Stepped until the physical time reaches final_time_s, like the splitting runs;
# a terminal ContinuousCallback made each solve about ten times slower.
function solve_ks_tsit5(problem, tolerance, final_time_s; dense::Bool)
    integrator = init(problem, Tsit5(); abstol=tolerance, reltol=tolerance,
        save_everystep=dense, dense=dense)
    while integrator.u[10] < final_time_s
        step!(integrator)
    end
    return integrator.sol, integrator.stats
end

function ks_tsit5_rows!(rows, initial, final_time_s, reference)
    s_span = 1.5 * ORBIT_COUNT * ks_splitting_step_size(initial, 1)
    problem = ODEProblem(ks_j2_rhs, initial, (0.0, s_span), PARAMS)
    for tolerance in RK_TOLERANCES
        runtime = median_runtime(() -> solve_ks_tsit5(problem, tolerance, final_time_s; dense=false))
        _, stats = solve_ks_tsit5(problem, tolerance, final_time_s; dense=false)
        solution, _ = solve_ks_tsit5(problem, tolerance, final_time_s; dense=true)
        cartesian = ks_state_to_cartesian.(solution.(
            range(solution.t[1], solution.t[end]; length=ACCURACY_SAMPLES)))
        error = maximum_position_error(getproperty.(cartesian, :elapsed_time_s),
            getproperty.(cartesian, :position_ii_m), reference)
        push!(rows, row("KS nonlinear Tsit5", "tolerance", tolerance, runtime,
            stats.naccept, stats.nf, error))
    end
end

function cartesian_rows!(rows, method, algorithm, x0, final_time_s, reference)
    problem = ODEProblem(cartesian_j2_rhs, x0, (0.0, final_time_s), PARAMS)
    times = range(0.0, final_time_s; length=ACCURACY_SAMPLES)
    for tolerance in RK_TOLERANCES
        work_solve = () -> solve(problem, algorithm; abstol=tolerance,
            reltol=tolerance, save_everystep=false, dense=false)
        runtime = median_runtime(work_solve)
        stats = work_solve().stats
        solution = solve(problem, algorithm; abstol=tolerance, reltol=tolerance)
        positions = [SVector{3}(view(x, 1:3)) for x in solution.(times)]
        push!(rows, row(method, "tolerance", tolerance, runtime, stats.naccept,
            stats.nf, maximum_position_error(times, positions, reference)))
    end
end

function main()
    orbit = InitialCondition(ra=PARAMS.Re + 3000.0e3, rp=PARAMS.Re + 300.0e3,
        i=51.6, ν=0.0)
    position, velocity = orbital_elements_to_cartesian(orbit, PLANET)
    x0 = SVector{6,Float64}(position..., velocity...)
    period_s = 2pi * sqrt(orbit.a^3 / PARAMS.μ)
    final_time_s = ORBIT_COUNT * period_s

    # The splitting and KS Tsit5 runs end up to one step past final_time_s.
    reference_problem = ODEProblem(cartesian_j2_rhs, x0, (0.0, final_time_s + period_s), PARAMS)
    reference = solve(reference_problem, Vern9(); abstol=REFERENCE_TOLERANCE,
        reltol=REFERENCE_TOLERANCE, maxiters=10^8)
    tighter = solve(reference_problem, Vern9(); abstol=REFERENCE_TOLERANCE / 10,
        reltol=REFERENCE_TOLERANCE / 10, maxiters=10^8)
    check_times = range(0.0, final_time_s + period_s; length=ACCURACY_SAMPLES)
    reference_error = maximum_position_error(check_times,
        [SVector{3}(view(x, 1:3)) for x in tighter.(check_times)], reference)

    initial = SVector{10,Float64}(cartesian_to_ks_state(position, velocity, PARAMS))
    rows = NamedTuple[]
    splitting_rows!(rows, initial, final_time_s, reference)
    ks_tsit5_rows!(rows, initial, final_time_s, reference)
    cartesian_rows!(rows, "Cartesian Tsit5", Tsit5(), x0, final_time_s, reference)
    cartesian_rows!(rows, "Cartesian DP8", DP8(), x0, final_time_s, reference)

    mkpath(RESULTS_DIR)
    raw = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "raw_work_precision.csv"), raw)
    CSV.write(joinpath(RESULTS_DIR, "metadata.csv"), DataFrame(
        quantity=["generated_at", "orbit_count", "period_s", "reference_tolerance",
            "reference_vs_tighter_reference_max_position_m"],
        value=string.([now(), ORBIT_COUNT, period_s, REFERENCE_TOLERANCE, reference_error])))

    println("Reference vs tighter reference: $(reference_error) m")
    show(stdout, MIME("text/plain"), raw; allrows=true)
    println()
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
