# Surrogate-GRAM variant of benchmarks/studies/leo_2048_constellation_gram_scaling.jl
# for a direct runtime comparison against real/native GRAM. Only the density
# model differs (GRAMAtmosphereModelSurrogate instead of GRAMAtmosphereModel);
# everything else (N_SATS, altitude, mission time, rhs_mode=flat threading) is
# held identical to the "parallel" mode real-GRAM case for an apples-to-apples
# comparison.
#
#   julia --project=. --threads=4 benchmarks/studies/gram_mars_fix_and_constellation_scaling/earth_surrogate_gram_scaling.jl

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO_ROOT, "examples", "common.jl"))
ensure_gramsuite_loaded!()
const GRAMAtmosphereModel = SimulationModel.GRAMAtmosphereModel
const GRAMAtmosphereModelSurrogate = SimulationModel.GRAMAtmosphereModelSurrogate

using LinearAlgebra
using StaticArrays
using Statistics

const N_SATS = parse(Int, get(ENV, "LEO_SCALING_N_SATS", "32"))
const ALT_M = 150e3
const MISSION_TIME_S = parse(Float64, get(ENV, "LEO_SCALING_MISSION_TIME_S", "600.0"))
const N_REPEATS = 1

function build_constellation_config()
    planet = Earth("", SPICE_PATH)
    harmonics = GravitationalHarmonicsModel(
        20, 20, joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv"), planet
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
            density_model=GRAMAtmosphereModelSurrogate(planet_name="earth"),
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

# Same routing knobs as the real-GRAM "parallel" case (rhs_mode=flat), minus
# the native-GRAM-specific vacuum-cache/freeze-per-step knobs, which don't
# apply to the surrogate path (surrogate_density_state is a pure array
# interpolation with no native calls to serialize or noise to freeze).
const ENV_PAIRS = [
    "SPACEAGORA_RHS_EXECUTION_MODE" => "flat",
    "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "1",
    "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
    "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
    "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
    "SPACEAGORA_RHS_CALIBRATE" => "off",
    "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
]

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

println("mode=surrogate  julia_threads=$(Threads.nthreads())  planet=earth  n_sat=$(N_SATS)  alt_km=$(ALT_M/1e3)  mission_s=$(MISSION_TIME_S)")
flush(stdout)

result = timed_run("surrogate", ENV_PAIRS)
println("median wall time (surrogate GRAM, $(Threads.nthreads()) threads): $(round(result.median_s; digits=3)) s")
