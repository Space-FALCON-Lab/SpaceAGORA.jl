# Mars variant of leo_2048_constellation_gram_scaling.jl, at the same
# altitude above the surface (150 km) as the Earth case.
#
# This script only works because of the GRAM Mars ephemeris-state bypass
# (GRAMSuite.jl's _GRAM_EPHEMERIS_STATE_FN hook, populated by
# ext/SpaceAGORAGRAMSuiteExt.jl's _gram_spice_ephemeris_state). Without that
# fix, every Mars density query fails immediately with "GRAM update failed
# (code=1): Error: A Spice error occurred. Error in longitude of the sun
# calculation." -- libGRAM.dylib statically links its own private CSPICE
# instance, isolated from SpaceAGORA's own SPICE.jl, and that isolated
# instance's default kernels only resolve solar geometry for Earth/Venus (see
# [[project_gram_mars_isolated_cspice_bug]] in memory for the full story).
# This script is a working end-to-end demonstration of that fix, and a
# regression check: if it starts failing again with the error above, the
# ephemeris bypass has broken.
#
#   julia --project=. --threads=4 mars_150km_gram_scaling.jl parallel
#
# Same modes as the Earth original: serial | parallel | accuracy | ensemble.
# See leo_2048_constellation_gram_scaling.jl's own header for the full
# rationale behind gram_method/freeze-per-step/vacuum-cache knobs -- unchanged
# here, only the planet/harmonics-file/GRAM planet_name were swapped to Mars.

mode = length(ARGS) >= 1 ? ARGS[1] : ""
mode in ("serial", "parallel", "accuracy", "ensemble") || error(
    "Usage: julia --project=. --threads=<N> $(basename(@__FILE__)) <serial|parallel|accuracy|ensemble>\n" *
    "  serial   : standard/direct GRAM, run with --threads=1\n" *
    "  parallel : GRAM trajectory look-ahead cache, run with --threads=8\n" *
    "  accuracy : compare standard vs. look-ahead trajectories in one process\n" *
    "  ensemble : run_constellation_ensemble (one worker per satellite), run with --threads=8"
)

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(REPO_ROOT, "examples", "common.jl"))
ensure_gramsuite_loaded!()
const GRAMAtmosphereModel = SimulationModel.GRAMAtmosphereModel

using LinearAlgebra
using StaticArrays
using Statistics

include(joinpath(REPO_ROOT, "benchmarks", "studies", "parallelization_performance", "trajectory_parity.jl"))

const N_SATS = parse(Int, get(ENV, "LEO_SCALING_N_SATS", "1024"))
const ALT_M = 150e3  # same altitude above the surface as the Earth case
const MISSION_TIME_S = parse(Float64, get(ENV, "LEO_SCALING_MISSION_TIME_S", "600.0"))
const N_REPEATS = 1

function build_constellation_config()
    planet = Mars("", SPICE_PATH)
    harmonics = GravitationalHarmonicsModel(
        20, 20, joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "Mars50c.csv"), planet
    )

    spacecraft = SpacecraftModel[]
    for i in 1:N_SATS
        root = Link{0}(root=true, m=500.0, ref_area=12.0)
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
            density_model=GRAMAtmosphereModel(planet_name="mars"),
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

function env_pairs_for(; rhs_mode::String, gram_method::String)
    gram_method in ("standard", "lookahead") || error("gram_method must be \"standard\" or \"lookahead\"")
    vacuum_enabled = gram_method == "lookahead" ? "1" : "0"
    freeze_per_step = gram_method == "standard" ? "1" : "0"
    density_parallel = gram_method == "standard" ? "off" : (rhs_mode == "serial" ? "off" : "auto")
    return [
        "SPACEAGORA_RHS_EXECUTION_MODE" => rhs_mode,
        "SPACEAGORA_HARMONICS_BATCH_ENABLED" => (rhs_mode == "serial" ? "0" : "1"),
        "SPACEAGORA_EFFECTOR_PARALLEL" => (rhs_mode == "serial" ? "off" : "auto"),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => density_parallel,
        "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => freeze_per_step,
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_RHS_CALIBRATE" => "off",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
        "SPACEAGORA_VACUUM_GRAM_CACHE" => vacuum_enabled,
        "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "20",
        "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => string(MISSION_TIME_S + 500.0),
        "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1e8",
    ]
end

function run_once(args)
    result = SpaceAGORA.run_simulation(args; isolate_state=false, return_solution=true, return_solver_metadata=true)
    return result.solution
end

args = build_constellation_config()

function timed_run(label::String, pairs)
    sol_warmup = withenv(pairs...) do
        run_once(args)
    end
    println("$(label) warmup retcode=$(sol_warmup.retcode)")
    flush(stdout)
    string(sol_warmup.retcode) == "Success" ||
        @warn "$(label) warmup run did not report Success -- results below may not reflect a full propagation."

    GC.gc()
    times = Float64[]
    local last_sol
    for r in 1:N_REPEATS
        GC.gc()
        t = @elapsed begin
            last_sol = withenv(pairs...) do
                run_once(args)
            end
        end
        push!(times, t)
        println("$(label) repeat $(r): $(round(t; digits=3)) s  retcode=$(last_sol.retcode)")
        flush(stdout)
    end
    median_t = sort(times)[cld(length(times), 2)]
    return (solution=last_sol, median_s=median_t)
end

println("mode=$(mode)  julia_threads=$(Threads.nthreads())  planet=mars  n_sat=$(N_SATS)  alt_km=$(ALT_M/1e3)  mission_s=$(MISSION_TIME_S)")
flush(stdout)

function timed_ensemble_run(label::String, pairs)
    ensemble_once = () -> SpaceAGORA.run_constellation_ensemble(
        args; threads=Threads.nthreads(), return_solution=true, return_solver_metadata=true
    )

    warmup = withenv(pairs...) do
        ensemble_once()
    end
    n_ok = length(warmup.successful)
    println("$(label) warmup: $(n_ok)/$(N_SATS) members succeeded")
    flush(stdout)
    n_ok == N_SATS || @warn "$(label) warmup: $(length(warmup.failed)) member(s) failed -- results below may not reflect a full propagation."

    GC.gc()
    times = Float64[]
    local last_result
    for r in 1:N_REPEATS
        GC.gc()
        t = @elapsed begin
            last_result = withenv(pairs...) do
                ensemble_once()
            end
        end
        push!(times, t)
        println("$(label) repeat $(r): $(round(t; digits=3)) s  $(length(last_result.successful))/$(N_SATS) succeeded")
        flush(stdout)
    end
    median_t = sort(times)[cld(length(times), 2)]
    return (result=last_result, median_s=median_t)
end

if mode in ("serial", "parallel")
    rhs_mode = mode == "serial" ? "serial" : "flat"
    gram_method = mode == "serial" ? "standard" : "lookahead"
    pairs = env_pairs_for(rhs_mode=rhs_mode, gram_method=gram_method)
    println("gram_method=$(gram_method)")
    flush(stdout)
    result = timed_run(mode, pairs)
    println("median wall time ($(mode), $(gram_method) GRAM, $(Threads.nthreads()) threads): $(round(result.median_s; digits=3)) s")
elseif mode == "ensemble"
    pairs = env_pairs_for(rhs_mode="serial", gram_method="standard")
    result = timed_ensemble_run("ensemble", pairs)
    println("median wall time (ensemble, standard GRAM, $(Threads.nthreads()) outer workers): $(round(result.median_s; digits=3)) s")
else # mode == "accuracy"
    pairs_standard = env_pairs_for(rhs_mode="flat", gram_method="standard")
    pairs_lookahead = env_pairs_for(rhs_mode="flat", gram_method="lookahead")

    result_ref = timed_run("standard (reference)", pairs_standard)
    result_cmp = timed_run("lookahead", pairs_lookahead)

    println()
    println("median wall time (standard):  $(round(result_ref.median_s; digits=3)) s")
    println("median wall time (lookahead): $(round(result_cmp.median_s; digits=3)) s")

    comparison = ppc_compare_trajectories(result_ref.solution, result_cmp.solution, args; sample_count=31)
    println()
    println("Trajectory parity (standard GRAM = reference, look-ahead = candidate), $(N_SATS) satellites, $(comparison.samples) time samples:")
    println("  position error:  rms=$(comparison.pos_rel_rms)  p90=$(comparison.pos_rel_p90)  max=$(comparison.pos_rel_max)  (relative)")
    println("  velocity error:  rms=$(comparison.vel_rel_rms)  p90=$(comparison.vel_rel_p90)  max=$(comparison.vel_rel_max)  (relative)")
    println("  quaternion angle max: $(comparison.q_angle_max_rad) rad")
    println("  omega error max: $(comparison.omega_rel_max) (relative)")
    println("  mass error max:  $(comparison.mass_rel_max) (relative)")
    println("  periapsis crossings: ref=$(comparison.ref_periapsis_count) cmp=$(comparison.cmp_periapsis_count)")
    println("  EI crossings:        ref=$(comparison.ref_interface_count) cmp=$(comparison.cmp_interface_count)")
    println("  event time delta max: $(comparison.event_time_abs_max_s) s")
    println("  PASS (within parallelization_performance's parity thresholds): $(comparison.pass)")
end
