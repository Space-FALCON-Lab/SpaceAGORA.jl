const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

using CSV
using DataFrames
using LinearAlgebra
using Random
using SPICE
using StaticArrays

include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
using .SimulationModel

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
include(joinpath(REPO_ROOT, "src", "simulation", "run_simulation.jl"))

struct StressMetrics
    rows::Int
    t_final::Float64
    pos_norm::Float64
    vel_norm::Float64
    mass_final::Float64
    mass_min::Float64
end

function make_spacecraft(rng::AbstractRNG, planet::Earth)
    root = Link{0}(root=true, m=160.0, ref_area=1.5)

    ra_alt = 530e3 + randn(rng) * 20e3
    rp_alt = 490e3 + randn(rng) * 20e3
    if rp_alt >= ra_alt
        ra_alt = rp_alt + 8e3
    end

    ic = InitialCondition(
        ra=planet.Rp_e + ra_alt,
        rp=planet.Rp_e + rp_alt,
        i=28.0 + randn(rng) * 0.8,
        ω=15.0 + randn(rng) * 2.0,
        Ω=20.0 + randn(rng) * 2.0,
        ν=160.0 + randn(rng) * 8.0
    )

    return SpacecraftModel(
        joints=Joint[],
        links=Link[root],
        root=root,
        instant_actuation=true,
        prop_mass=20.0,
        inertia_tensor=root.inertia,
        n_reaction_wheels=0,
        n_thrusters=0,
        initial_condition=ic,
        id=1
    )
end

function build_case(seed::Int)
    rng = MersenneTwister(seed)
    spice_path = joinpath(REPO_ROOT, "GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
    planet = Earth("", spice_path)
    spacecraft = make_spacecraft(rng, planet)

    thrust = 0.6 + rand(rng) * 0.5
    burn_time = 180.0 + rand(rng) * 40.0
    isp = 280.0 + rand(rng) * 40.0
    thruster = BaseThrusterModel(
        thrust=[thrust],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[0.0],
        stop_burn_time=[burn_time],
        Isp=[isp]
    )

    return SimulationConfiguration(
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
            mission_time=900.0,
            orientation_sim=false,
            num_steps_to_save=400
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel([spacecraft], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(thruster,), control_rates=[1.0]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
end

function run_seed(seed::Int)
    args = build_case(seed)
    return mktempdir() do tmp
        cd(tmp) do
            run_simulation(args; isolate_state=true)
            isfile("simulation_results.csv") || error("Missing simulation_results.csv for seed=$seed")
            df = CSV.read("simulation_results.csv", DataFrame)

            nrow(df) >= 40 || error("Too few rows for seed=$seed: $(nrow(df))")
            required_cols = ("time", "sc1_pos_1", "sc1_pos_2", "sc1_pos_3", "sc1_vel_1", "sc1_vel_2", "sc1_vel_3", "sc1_mass")
            for col in required_cols
                col in names(df) || error("Missing column $col for seed=$seed")
                all(isfinite, Float64.(df[!, col])) || error("Non-finite values in column $col for seed=$seed")
            end

            mass = Float64.(df.sc1_mass)
            minimum(mass) > 0.0 || error("Non-positive mass encountered for seed=$seed")
            minimum(mass) < mass[1] || error("Mass did not decrease for seed=$seed")

            pos = SVector{3, Float64}(Float64(df.sc1_pos_1[end]), Float64(df.sc1_pos_2[end]), Float64(df.sc1_pos_3[end]))
            vel = SVector{3, Float64}(Float64(df.sc1_vel_1[end]), Float64(df.sc1_vel_2[end]), Float64(df.sc1_vel_3[end]))
            return StressMetrics(
                nrow(df),
                Float64(df.time[end]),
                norm(pos),
                norm(vel),
                mass[end],
                minimum(mass)
            )
        end
    end
end

function metrics_close(a::StressMetrics, b::StressMetrics)
    return a.rows == b.rows &&
           isapprox(a.t_final, b.t_final; atol=1e-9, rtol=0.0) &&
           isapprox(a.pos_norm, b.pos_norm; atol=0.0, rtol=1e-10) &&
           isapprox(a.vel_norm, b.vel_norm; atol=0.0, rtol=1e-10) &&
           isapprox(a.mass_final, b.mass_final; atol=0.0, rtol=1e-10) &&
           isapprox(a.mass_min, b.mass_min; atol=0.0, rtol=1e-10)
end

seeds = [11, 23, 37, 53, 71]
metrics = Dict{Int, StressMetrics}()

for seed in seeds
    m1 = run_seed(seed)
    m2 = run_seed(seed)
    metrics_close(m1, m2) || error("Deterministic replay mismatch for seed=$seed")
    metrics[seed] = m1
    println("nightly_mc_seed_ok seed=$(seed) rows=$(m1.rows) tf=$(m1.t_final) |r|=$(m1.pos_norm) |v|=$(m1.vel_norm) m=$(m1.mass_final)")
end

rounded_position_norms = Set(round(metrics[s].pos_norm; sigdigits=8) for s in seeds)
length(rounded_position_norms) > 1 || error("Monte Carlo stress has no cross-seed diversity in final position norms.")

println("nightly_montecarlo_stress_ok seeds=$(length(seeds))")
