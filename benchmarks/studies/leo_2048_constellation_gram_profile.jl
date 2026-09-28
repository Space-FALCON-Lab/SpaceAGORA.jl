# Re-profile of the 1024-satellite LEO GRAM constellation scenario (see
# LEO_GRAM_CONSTELLATION_HANDOFF.md) AFTER the routing/sort/hot-path fixes landed.
# The pre-fix profile (in the handoff doc) showed 51% utilization, ~90% of samples
# waiting on the persistent-worker dispatch channel -- caused by the 2-effector
# flat-queue routing bug, since fixed. This script re-runs the same scenario
# (8 threads, parallel/flat mode) under Profile to see whether utilization improved
# and what the next-largest remaining cost is.
#
# Root cause found in that re-profile: `auto_thread_min_budget(:density_callback)`
# (src/parallel/policy/env_config.jl) hardcodes a 16-thread minimum before
# SPACEAGORA_DENSITY_CALLBACK_PARALLEL=auto will ever enable threading, so at
# --threads=8 the entire per-satellite GRAM atmosphere sample loop
# (_prefill_environment_samples!, dynamics_rhs.jl) runs serially on the calling
# thread regardless of satellite count. Pass "on" as the first arg to validate the
# fix (SPACEAGORA_DENSITY_CALLBACK_PARALLEL=on bypasses that gate); omit/pass "auto"
# to reproduce the original (broken) baseline.
#
#   julia --project=. --threads=8 benchmarks/studies/leo_2048_constellation_gram_profile.jl auto
#   julia --project=. --threads=8 benchmarks/studies/leo_2048_constellation_gram_profile.jl on

density_mode = length(ARGS) >= 1 ? ARGS[1] : "auto"
density_mode in ("auto", "on") || error("Usage: ... leo_2048_constellation_gram_profile.jl <auto|on>")

include(joinpath(@__DIR__, "..", "..", "examples", "common.jl"))
ensure_gramsuite_loaded!()
const GRAMAtmosphereModel = SimulationModel.GRAMAtmosphereModel

using LinearAlgebra
using Profile

const N_SATS = 1024
const ALT_M = 150e3
const MISSION_TIME_S = 3600.0

function build_constellation_config()
    planet = Earth("", SPICE_PATH)
    harmonics = GravitationalHarmonicsModel(
        20, 20, joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv"), planet
    )

    spacecraft = SpacecraftModel[]
    for i in 1:N_SATS
        root = Link(root=true, m=500.0, ref_area=12.0)
        jitter = 50.0 * (i - 1) / N_SATS
        ic = InitialCondition(
            ra=planet.Rp_e + ALT_M + jitter,
            rp=planet.Rp_e + ALT_M + jitter,
            i=53.0,
            ω=0.0,
            Ω=10.0,
            ν=360.0 * (i - 1) / N_SATS
        )
        push!(spacecraft, SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, i))
    end

    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=MISSION_TIME_S,
            orientation_sim=false,
            num_steps_to_save=20
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=300.0,
            density_model=GRAMAtmosphereModel(planet_name="earth"),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel(spacecraft, (harmonics, AerodynamicCoefficientfM())),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=2.0)
    )
end

const ENV_PAIRS = [
    "SPACEAGORA_RHS_EXECUTION_MODE" => "flat",
    "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "1",
    "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
    "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => density_mode,
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
    "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
    "SPACEAGORA_RHS_CALIBRATE" => "off",
    "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
    "SPACEAGORA_VACUUM_GRAM_CACHE" => "1",
    "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "20",
    "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => string(MISSION_TIME_S + 500.0),
    "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1e8",
]

function run_once(args)
    result = SpaceAGORA.run_simulation(args; isolate_state=false, return_solution=true, return_solver_metadata=true)
    return result.solution
end

args = build_constellation_config()

println("threads=$(Threads.nthreads())  n_sat=$(N_SATS)  alt_km=$(ALT_M/1e3)  mission_s=$(MISSION_TIME_S)  density_mode=$(density_mode)")
flush(stdout)

println("Warmup (unprofiled, absorbs JIT compilation)...")
flush(stdout)
t_warmup = @elapsed begin
    sol_warmup = withenv(ENV_PAIRS...) do
        run_once(args)
    end
end
println("warmup retcode=$(sol_warmup.retcode)  elapsed=$(round(t_warmup; digits=1))s")
flush(stdout)
string(sol_warmup.retcode) == "Success" || @warn "Warmup run did not report Success."

GC.gc()
Profile.clear()
Profile.init(n=10_000_000, delay=0.001)

println("Profiled run starting...")
flush(stdout)
local sol
t_profiled = @elapsed begin
    global sol = Profile.@profile withenv(ENV_PAIRS...) do
        run_once(args)
    end
end
println("profiled run retcode=$(sol.retcode)  elapsed=$(round(t_profiled; digits=1))s")
flush(stdout)

outdir = joinpath(REPO_ROOT, "output")
mkpath(outdir)

overhead_path = joinpath(outdir, "leo_2048_gram_profile_overhead_density_$(density_mode).txt")
open(overhead_path, "w") do io
    Profile.print(io; format=:flat, C=true, sortedby=:overhead, mincount=20)
end
println("wrote $(overhead_path)")

count_path = joinpath(outdir, "leo_2048_gram_profile_count_density_$(density_mode).txt")
open(count_path, "w") do io
    Profile.print(io; format=:flat, C=true, sortedby=:count, mincount=40)
end
println("wrote $(count_path)")

thread_path = joinpath(outdir, "leo_2048_gram_profile_by_thread_density_$(density_mode).txt")
open(thread_path, "w") do io
    Profile.print(io; format=:flat, C=true, sortedby=:overhead, groupby=:thread, mincount=20)
end
println("wrote $(thread_path)")

# Also print the summary utilization line to stdout directly (Profile.print's final
# "Total snapshots: N. Utilization: X%" line), so it shows up in captured run output
# without needing to open the files.
println()
println("=== Overhead (self-time) summary, all threads ===")
Profile.print(format=:flat, C=true, sortedby=:overhead, mincount=200)
flush(stdout)
