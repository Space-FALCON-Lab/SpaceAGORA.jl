##
# Empirically checks whether `SimulationSettings.normalize` (a legacy flag
# from an earlier nondimensional-units era of the engine) still changes
# propagation accuracy or runtime in the current typed `run_simulation` path.
#
# Source inspection (src/simulation/engine/execution.jl,
# src/simulation/engine/reporting.jl) shows `normalize=true` either throws
# an ArgumentError, or — with SPACEAGORA_ALLOW_TYPED_NORMALIZE=1 set — just
# prints a transition warning and then propagates in SI units exactly as
# `normalize=false` does. There is no surviving code path that actually
# nondimensionalizes state. This script verifies that claim against real
# runs rather than static reading of the source: (1) confirms the reject
# case actually throws, (2) runs the override case against the baseline and
# compares final trajectory state and wall-clock time.
#
# The quaternion unit-norm re-projection callback (a separate, unrelated,
# currently-active feature keyed to `orientation_sim`) is deliberately left
# untouched here -- this scenario uses orientation_sim=false throughout, so
# that callback never registers for either arm.
#
# Usage: julia --project=. benchmarks/studies/normalize_flag_noop_check.jl
##

include(joinpath(@__DIR__, "..", "..", "examples", "common.jl"))

using Statistics
using Printf
using CSV
using DataFrames

const N_REPEATS = 5
const RESULTS_DIR = joinpath(REPO_ROOT, "output", "performance", "normalize_flag_noop_check")
mkpath(RESULTS_DIR)

function _build_config(; normalize::Bool, results_directory::String)
    planet = make_no_gram_planet(:earth)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.05, 2.05, 2.8),
        panel_dims=(0.01, 5.7 / 2.0, 1.0),
        bus_mass=620.0,
        panel_mass_each=10.0,
        panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
        ic=InitialCondition(
            ra=56_378.7978559e3,
            rp=planet.Rp_e + 200_590.0,
            i=89.876,
            ω=75.505,
            Ω=104.115,
            ν=175.0
        ),
        prop_mass=200.0,
        id=1
    )

    args = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=12.0 * 3600.0,
        initial_time=InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        orientation_sim=false,
        keplerian=true,
        EI_km=120.0,
        verbose=false,
        results_directory=results_directory
    )

    # make_example_config always sets normalize=false; override explicitly.
    return SimulationModel.SimulationConfiguration(
        simulation_settings=SimulationModel.SimulationSettings(
            results=args.simulation_settings.results,
            verbose=args.simulation_settings.verbose,
            generate_plots=false,
            results_directory=args.simulation_settings.results_directory,
            normalize=normalize
        ),
        mission_configuration=args.mission_configuration,
        environment_model=args.environment_model,
        dynamics_model=args.dynamics_model,
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances,
        solver_config=args.solver_config
    )
end

function _run_timed(args::SimulationModel.SimulationConfiguration)
    t = @elapsed run_simulation(args)
    csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
    df = CSV.read(csv_path, DataFrame)
    return (time_s=t, df=df)
end

function _final_state_vector(df::DataFrame)
    row = df[end, :]
    cols = [:sc1_pos_1, :sc1_pos_2, :sc1_pos_3, :sc1_vel_1, :sc1_vel_2, :sc1_vel_3]
    return Float64[row[c] for c in cols]
end

println("=" ^ 70)
println("normalize flag no-op check")
println("=" ^ 70)

# --- Part 1: does normalize=true still throw without the override flag? ---
delete!(ENV, "SPACEAGORA_ALLOW_TYPED_NORMALIZE")
reject_dir = joinpath(RESULTS_DIR, "reject_case")
mkpath(reject_dir)
reject_args = _build_config(normalize=true, results_directory=reject_dir)

reject_threw, reject_message = let reject_threw = false, reject_message = ""
    try
        run_simulation(reject_args)
    catch err
        reject_threw = err isa ArgumentError
        reject_message = sprint(showerror, err)
    end
    (reject_threw, reject_message)
end

println()
println("Part 1: normalize=true, SPACEAGORA_ALLOW_TYPED_NORMALIZE unset")
println("  Threw ArgumentError: $reject_threw")
if reject_threw
    println("  Message: $reject_message")
else
    println("  UNEXPECTED: did not throw. Message/state: $reject_message")
end

# --- Part 2: with the override, is propagation actually different? ---
ENV["SPACEAGORA_ALLOW_TYPED_NORMALIZE"] = "1"

baseline_times = Float64[]
baseline_final = Vector{Float64}[]
for i in 1:N_REPEATS
    dir = joinpath(RESULTS_DIR, "baseline_run$(i)")
    mkpath(dir)
    args = _build_config(normalize=false, results_directory=dir)
    result = _run_timed(args)
    push!(baseline_times, result.time_s)
    push!(baseline_final, _final_state_vector(result.df))
end

normalized_times = Float64[]
normalized_final = Vector{Float64}[]
for i in 1:N_REPEATS
    dir = joinpath(RESULTS_DIR, "normalized_run$(i)")
    mkpath(dir)
    args = _build_config(normalize=true, results_directory=dir)
    result = _run_timed(args)
    push!(normalized_times, result.time_s)
    push!(normalized_final, _final_state_vector(result.df))
end

delete!(ENV, "SPACEAGORA_ALLOW_TYPED_NORMALIZE")

println()
println("Part 2: normalize=false vs. normalize=true (with override), $N_REPEATS repeats each")
println("-" ^ 70)
@printf("  normalize=false median time: %.4f s  (min %.4f, max %.4f)\n",
    median(baseline_times), minimum(baseline_times), maximum(baseline_times))
@printf("  normalize=true  median time: %.4f s  (min %.4f, max %.4f)\n",
    median(normalized_times), minimum(normalized_times), maximum(normalized_times))

# Final-state agreement: compare each normalize=true run's final position/
# velocity vector against the baseline median run, componentwise.
baseline_ref = baseline_final[1]
max_pos_diff_m = 0.0
max_vel_diff_ms = 0.0
for (bf, nf) in zip(baseline_final, normalized_final)
    pos_diff = maximum(abs.(bf[1:3] .- nf[1:3]))
    vel_diff = maximum(abs.(bf[4:6] .- nf[4:6]))
    global max_pos_diff_m = max(max_pos_diff_m, pos_diff)
    global max_vel_diff_ms = max(max_vel_diff_ms, vel_diff)
end
# Also compare across baseline's own repeats, to give the position/velocity
# diffs above a run-to-run noise floor to be judged against.
baseline_self_pos_diff = 0.0
baseline_self_vel_diff = 0.0
for bf in baseline_final
    pos_diff = maximum(abs.(bf[1:3] .- baseline_ref[1:3]))
    vel_diff = maximum(abs.(bf[4:6] .- baseline_ref[4:6]))
    global baseline_self_pos_diff = max(baseline_self_pos_diff, pos_diff)
    global baseline_self_vel_diff = max(baseline_self_vel_diff, vel_diff)
end

println()
@printf("  max |Δposition| (normalize=true vs. false), any component: %.3e m\n", max_pos_diff_m)
@printf("  max |Δvelocity| (normalize=true vs. false), any component: %.3e m/s\n", max_vel_diff_ms)
@printf("  baseline run-to-run noise floor, position: %.3e m\n", baseline_self_pos_diff)
@printf("  baseline run-to-run noise floor, velocity: %.3e m/s\n", baseline_self_vel_diff)

println()
println("=" ^ 70)
println("Summary")
println("=" ^ 70)
println("  normalize=true without the override flag: ",
    reject_threw ? "rejected (ArgumentError), simulation does not run" : "did NOT reject (unexpected)")
indistinguishable = max_pos_diff_m <= 10 * baseline_self_pos_diff && max_vel_diff_ms <= 10 * baseline_self_vel_diff
println("  normalize=true with the override vs. normalize=false: ",
    indistinguishable ? "indistinguishable from normalize=false within run-to-run noise" :
    "measurably different -- see diffs above, does not match source-level expectation")
println("  Runtime difference is not attributable to the flag beyond run-to-run noise: ",
    abs(median(normalized_times) - median(baseline_times)) <= max(std(baseline_times), std(normalized_times)))

summary_path = joinpath(RESULTS_DIR, "summary.txt")
open(summary_path, "w") do io
    println(io, "normalize=true without override threw ArgumentError: $reject_threw")
    if reject_threw
        println(io, "Rejection message: $reject_message")
    end
    println(io, "normalize=false median time (s): $(median(baseline_times))")
    println(io, "normalize=true  median time (s): $(median(normalized_times))")
    println(io, "max |Δposition| (m): $max_pos_diff_m")
    println(io, "max |Δvelocity| (m/s): $max_vel_diff_ms")
    println(io, "baseline self noise floor position (m): $baseline_self_pos_diff")
    println(io, "baseline self noise floor velocity (m/s): $baseline_self_vel_diff")
end
println()
println("Summary written to $(abspath(summary_path))")
