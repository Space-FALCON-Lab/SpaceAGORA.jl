const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using CSV
using DataFrames
using SPICE
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SimulationModel = SpaceAGORA.SimulationModel
const quat_mult = SimulationModel.quat_mult
const run_simulation = SpaceAGORA.run_simulation

if Threads.nthreads() < 2
    error("Threaded smoke requires at least 2 Julia threads; got $(Threads.nthreads())")
end

spice_path = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
planet = Earth("", spice_path)

function make_sc(id::Int64, ν_deg::Float64)
    root = Link{0}(root=true, m=140.0, ref_area=1.2)
    ic = InitialCondition(
        ra=planet.Rp_e + 520e3,
        rp=planet.Rp_e + 500e3,
        i=28.0,
        ω=15.0,
        Ω=20.0,
        ν=ν_deg
    )
    return SpacecraftModel(
        joints=Joint[],
        links=Link[root],
        root=root,
        instant_actuation=true,
        prop_mass=15.0,
        inertia_tensor=root.inertia,
        n_reaction_wheels=0,
        n_thrusters=0,
        initial_condition=ic,
        id=id
    )
end

sc1 = make_sc(1, 170.0)
sc2 = make_sc(2, 160.0)

thruster = BaseThrusterModel(
    thrust=[0.9, 0.8],
    direction=[0.0, π],
    Δv=[0.0, 0.0],
    start_burn_time=[0.0, 0.0],
    stop_burn_time=[15.0, 15.0],
    Isp=[300.0, 280.0]
)

args = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        normalize=false
    ),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime,
        keplerian=true,
        number_of_orbits=1,
        mission_time=60.0,
        orientation_sim=false,
        num_steps_to_save=200
    ),
    environment_model=EnvironmentModel(
        planet=planet,
        EI=120.0,
        density_model=NoAtmosphereModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    ),
    dynamics_model=DynamicsModel([sc1, sc2], (InverseSquaredGravityModel(),)),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(thruster,), control_rates=[1.0]),
    initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
)

mktempdir() do tmp
    cd(tmp) do
        run_simulation(args)
        csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
        if !isfile(csv_path)
            error("Expected simulation_results.csv to be written by threaded smoke run")
        end
        df = CSV.read(csv_path, DataFrame)
        if nrow(df) < 6
            error("Threaded smoke produced too few rows: $(nrow(df))")
        end

        required_cols = [
            "sc1_pos_1", "sc1_vel_1", "sc1_mass",
            "sc2_pos_1", "sc2_vel_1", "sc2_mass"
        ]
        for col in required_cols
            if !(col in names(df))
                error("Missing required column $col in threaded smoke output")
            end
            if !all(isfinite, Float64.(df[!, col]))
                error("Non-finite values found in $col during threaded smoke")
            end
        end

        if !(minimum(Float64.(df.sc1_mass)) < Float64(df.sc1_mass[1]))
            error("Expected sc1 mass to decrease during threaded smoke run")
        end
        if !(minimum(Float64.(df.sc2_mass)) < Float64(df.sc2_mass[1]))
            error("Expected sc2 mass to decrease during threaded smoke run")
        end
    end
end

println("threaded_smoke_ok")

# ---------------------------------------------------------------------------
# Thread-count independence. The effector loop, the n-body loop and the
# multibody aero loop collect per-item results and sum them in a fixed order on
# one thread, so the same configuration must give bit-identical trajectories
# whether the effectors are evaluated serially or on 2 or 4 workers, and
# whether the n-body bodies are evaluated serially or on 2 workers. Three
# effectors with different costs make the partition non-trivial.
# ---------------------------------------------------------------------------
harmonics_file = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
det_effectors = (
    InverseSquaredGravityModel(),
    NBodyGravityModel(["Sun", "Moon"], "Earth", spice_path),
    GravitationalHarmonicsModel(4, 4, harmonics_file, planet),
)
det_args = SimulationConfiguration(
    simulation_settings=args.simulation_settings,
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime,
        keplerian=true,
        number_of_orbits=1,
        mission_time=300.0,
        orientation_sim=false,
        num_steps_to_save=50
    ),
    environment_model=args.environment_model,
    dynamics_model=DynamicsModel([sc1, sc2], det_effectors),
    guidance_model=args.guidance_model,
    navigation_model=args.navigation_model,
    control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
    initial_time=args.initial_time,
    integration_tolerances=args.integration_tolerances
)

function det_run(env_pairs::Vector{Pair{String, String}})::Matrix{Float64}
    return withenv(env_pairs...) do
        mktempdir() do tmp
            cd(tmp) do
                run_simulation(det_args)
                df = CSV.read(joinpath(det_args.simulation_settings.results_directory, "simulation_results.csv"), DataFrame)
                cols = [c for c in names(df) if eltype(df[!, c]) <: Real]
                return Matrix{Float64}(df[:, cols])
            end
        end
    end
end

det_force_on = [
    "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY" => "0",
    "SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "1",
    "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "1",
]
det_ref = det_run(["SPACEAGORA_EFFECTOR_PARALLEL" => "off", "SPACEAGORA_MULTIBODY_PARALLEL" => "off"])
size(det_ref, 1) >= 10 || error("Determinism check produced too few rows: $(size(det_ref, 1))")
det_variants = [
    "effectors on 2 workers" => vcat(det_force_on, ["SPACEAGORA_EFFECTOR_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_MAX_THREADS" => "2", "SPACEAGORA_MULTIBODY_PARALLEL" => "off"]),
    "effectors on 4 workers" => vcat(det_force_on, ["SPACEAGORA_EFFECTOR_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_MAX_THREADS" => "4", "SPACEAGORA_MULTIBODY_PARALLEL" => "off"]),
    "n-body on 2 workers"    => vcat(det_force_on, ["SPACEAGORA_EFFECTOR_PARALLEL" => "off", "SPACEAGORA_MULTIBODY_PARALLEL" => "on", "SPACEAGORA_MULTIBODY_MAX_THREADS" => "2"]),
    "everything on"          => vcat(det_force_on, ["SPACEAGORA_EFFECTOR_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_MAX_THREADS" => "4", "SPACEAGORA_MULTIBODY_PARALLEL" => "on", "SPACEAGORA_MULTIBODY_MAX_THREADS" => "2"]),
]
for (label, env_pairs) in det_variants
    out = det_run(env_pairs)
    if size(out) != size(det_ref)
        error("Threaded determinism: $label produced $(size(out)) rows/cols, serial produced $(size(det_ref))")
    end
    if !isequal(out, det_ref)
        worst = maximum(abs.(out .- det_ref))
        error("Threaded determinism: $label differs from the serial run (max abs difference $worst)")
    end
end
println("threaded_determinism_ok variants=$(length(det_variants)) rows=$(size(det_ref, 1))")
