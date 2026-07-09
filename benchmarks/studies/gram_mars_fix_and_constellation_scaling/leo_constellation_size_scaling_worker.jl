# Worker for leo_constellation_size_scaling.jl -- runs ONE (satellite count,
# density-model mode) configuration and prints its median wall time. Each
# invocation is a separate process (spawned by the driver script) so that the
# per-(N_SATS, mode) JIT specialization this workload requires -- ODEParams
# is parameterized on both N_sats and the density-model type, so each
# distinct combination compiles its own RHS/solver specialization -- doesn't
# accumulate compiled code and memory across the whole sweep in one process.
#
#   julia --project=. --threads=<N> leo_constellation_size_scaling_worker.jl <n_sats> <mode>
#
# mode is one of:
#   standard  -- real/native GRAM, no vacuum-predicted lookahead cache
#                (SPACEAGORA_DENSITY_FREEZE_PER_STEP=1, density calls
#                serialized -- see leo_2048_constellation_gram_scaling.jl's
#                own header for why concurrent direct native GRAM calls are
#                unsafe). Measured faster than the lookahead cache for this
#                600s mission length (see project findings): the lookahead
#                cache's upfront drag-free-trajectory + spline-fit cost isn't
#                amortized by a mission this short.
#   surrogate -- GRAM offline surrogate (GRAMAtmosphereModelSurrogate):
#                precomputed interpolation grid, no native calls, no lock.
#   no_gram   -- NoAtmosphereModel(), and the aerodynamic effector dropped
#                entirely -- baseline gravity-only dynamics cost, isolating
#                constellation-size scaling with zero atmosphere-model
#                overhead.

mode = length(ARGS) >= 2 ? ARGS[2] : ""
mode in ("standard", "surrogate", "no_gram") || error(
    "Usage: julia --project=. --threads=<N> $(basename(@__FILE__)) <n_sats> <standard|surrogate|no_gram>"
)
const N_SATS = parse(Int, ARGS[1])

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO_ROOT, "examples", "common.jl"))
ensure_gramsuite_loaded!()
const GRAMAtmosphereModel = SimulationModel.GRAMAtmosphereModel
const GRAMAtmosphereModelSurrogate = SimulationModel.GRAMAtmosphereModelSurrogate

const ALT_M = 150e3
const MISSION_TIME_S = parse(Float64, get(ENV, "LEO_SCALING_MISSION_TIME_S", "600.0"))
const N_REPEATS = 1

function build_constellation_config()
    planet = Earth("", SPICE_PATH)
    harmonics = GravitationalHarmonicsModel(
        20, 20, joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv"), planet
    )

    density_model = if mode == "standard"
        GRAMAtmosphereModel(planet_name="earth")
    elseif mode == "surrogate"
        GRAMAtmosphereModelSurrogate(planet_name="earth")
    else
        NoAtmosphereModel()
    end

    # no_gram drops the aerodynamic effector entirely (matching this repo's
    # own no-GRAM baseline convention, e.g. examples/AGORA_Earth_NoGRAM.jl) so
    # the comparison isolates atmosphere-model overhead rather than measuring
    # a drag effector that always computes zero force against a zero-density
    # model.
    effectors = mode == "no_gram" ? (harmonics,) : (harmonics, AerodynamicCoefficientfM())

    spacecraft = SpacecraftModel[]
    for i in 1:N_SATS
        root = Link{0}(root=true, m=500.0, ref_area=12.0)
        jitter = 50.0 * (i - 1) / max(N_SATS, 1)
        ic = InitialCondition(
            ra=planet.Rp_e + ALT_M + jitter,
            rp=planet.Rp_e + ALT_M + jitter,
            i=53.0,
            ω=0.0,
            Ω=10.0,
            ν=360.0 * (i - 1) / max(N_SATS, 1)
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
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel(spacecraft, effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=2.0)
    )
end

function env_pairs_for(mode::String)
    common = [
        "SPACEAGORA_RHS_EXECUTION_MODE" => "flat",
        "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "1",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_RHS_CALIBRATE" => "off",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
    ]
    if mode == "standard"
        return vcat(common, [
            "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "1",
            "SPACEAGORA_VACUUM_GRAM_CACHE" => "0",
        ])
    elseif mode == "surrogate"
        return vcat(common, ["SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto"])
    else
        return common
    end
end

function run_once(args)
    result = SpaceAGORA.run_simulation(args; isolate_state=false, return_solution=true, return_solver_metadata=true)
    return result.solution
end

args = build_constellation_config()
pairs = env_pairs_for(mode)

sol_warmup = withenv(pairs...) do
    run_once(args)
end
println("warmup retcode=$(sol_warmup.retcode)")
flush(stdout)
string(sol_warmup.retcode) == "Success" ||
    @warn "warmup run did not report Success -- results below may not reflect a full propagation."

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
    println("repeat $(r): $(round(t; digits=4)) s  retcode=$(last_sol.retcode)")
    flush(stdout)
end
median_t = sort(times)[cld(length(times), 2)]
println("median wall time (mode=$(mode), n_sat=$(N_SATS), $(Threads.nthreads()) threads): $(round(median_t; digits=4)) s")
