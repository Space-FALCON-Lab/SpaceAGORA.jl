const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

using CSV
using DataFrames
using SPICE

include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
using .SimulationModel

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
include(joinpath(REPO_ROOT, "src", "simulation", "run_simulation.jl"))

if Threads.nthreads() < 2
    error("Threaded smoke requires at least 2 Julia threads; got $(Threads.nthreads())")
end

spice_path = joinpath(REPO_ROOT, "GRAM_Data", "SPICE")
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
        if !isfile("simulation_results.csv")
            error("Expected simulation_results.csv to be written by threaded smoke run")
        end
        df = CSV.read("simulation_results.csv", DataFrame)
        if nrow(df) < 20
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
