const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using CSV
using DataFrames
using LinearAlgebra
using SPICE
using StaticArrays

if !isdefined(@__MODULE__, :SimulationModel)
    include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
end
using .SimulationModel

# SimulationEngine uses SimulationModel and provides the canonical run_simulation entrypoint.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :SimulationEngine)
    include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SimulationEngine.run_simulation
end

Base.@kwdef struct FlakeCase
    name::String
    ν_deg::Float64
    thrust::Float64
    burn_time::Float64
    isp::Float64
end

struct FlakeMetrics
    rows::Int
    t_final::Float64
    pos_norm::Float64
    vel_norm::Float64
    mass_final::Float64
end

const POS_RTOL = let raw = get(ENV, "SPACEAGORA_FLAKE_POS_RTOL", "1e-10")
    parsed = tryparse(Float64, raw)
    parsed === nothing && error("Invalid SPACEAGORA_FLAKE_POS_RTOL=$raw")
    parsed
end

const VEL_RTOL = let raw = get(ENV, "SPACEAGORA_FLAKE_VEL_RTOL", "1e-10")
    parsed = tryparse(Float64, raw)
    parsed === nothing && error("Invalid SPACEAGORA_FLAKE_VEL_RTOL=$raw")
    parsed
end

const MASS_RTOL = let raw = get(ENV, "SPACEAGORA_FLAKE_MASS_RTOL", "1e-12")
    parsed = tryparse(Float64, raw)
    parsed === nothing && error("Invalid SPACEAGORA_FLAKE_MASS_RTOL=$raw")
    parsed
end

const RUNS = let raw = get(ENV, "SPACEAGORA_FLAKE_GUARD_RUNS", "3")
    parsed = tryparse(Int, raw)
    parsed === nothing && error("Invalid SPACEAGORA_FLAKE_GUARD_RUNS=$raw")
    parsed < 2 && error("SPACEAGORA_FLAKE_GUARD_RUNS must be >= 2")
    parsed
end

const CASES = [
    FlakeCase(name="baseline", ν_deg=175.0, thrust=0.75, burn_time=110.0, isp=300.0),
    FlakeCase(name="retrograde_bias", ν_deg=160.0, thrust=0.90, burn_time=95.0, isp=290.0),
    FlakeCase(name="long_burn", ν_deg=150.0, thrust=0.60, burn_time=150.0, isp=320.0),
]

function build_case(case::FlakeCase)
    spice_path = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
    planet = Earth("", spice_path)

    root = Link{0}(root=true, m=150.0, ref_area=1.4)
    ic = InitialCondition(
        ra=planet.Rp_e + 520e3,
        rp=planet.Rp_e + 500e3,
        i=28.0,
        ω=15.0,
        Ω=20.0,
        ν=case.ν_deg
    )

    spacecraft = SpacecraftModel(
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

    thruster = BaseThrusterModel(
        thrust=[case.thrust],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[0.0],
        stop_burn_time=[case.burn_time],
        Isp=[case.isp]
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
            mission_time=240.0,
            orientation_sim=false,
            num_steps_to_save=250
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
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=8.0)
    )
end

function run_case(case::FlakeCase)
    args = build_case(case)
    return mktempdir() do tmp
        cd(tmp) do
            run_simulation(args; isolate_state=true)
            csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            isfile(csv_path) || error("Missing simulation_results.csv")
            df = CSV.read(csv_path, DataFrame)
            nrow(df) >= 25 || error("Too few output rows: $(nrow(df))")

            required = ("time", "sc1_pos_1", "sc1_pos_2", "sc1_pos_3", "sc1_vel_1", "sc1_vel_2", "sc1_vel_3", "sc1_mass")
            for col in required
                col in names(df) || error("Missing required column $col")
                all(isfinite, Float64.(df[!, col])) || error("Non-finite values in column $col")
            end

            mass = Float64.(df.sc1_mass)
            minimum(mass) > 0.0 || error("Encountered non-positive mass")
            minimum(mass) < mass[1] || error("Mass never decreased")

            pos = SVector{3, Float64}(Float64(df.sc1_pos_1[end]), Float64(df.sc1_pos_2[end]), Float64(df.sc1_pos_3[end]))
            vel = SVector{3, Float64}(Float64(df.sc1_vel_1[end]), Float64(df.sc1_vel_2[end]), Float64(df.sc1_vel_3[end]))

            return FlakeMetrics(
                nrow(df),
                Float64(df.time[end]),
                norm(pos),
                norm(vel),
                mass[end]
            )
        end
    end
end

function metrics_consistent(a::FlakeMetrics, b::FlakeMetrics)
    return a.rows == b.rows &&
           isapprox(a.t_final, b.t_final; atol=1e-9, rtol=0.0) &&
           isapprox(a.pos_norm, b.pos_norm; atol=0.0, rtol=POS_RTOL) &&
           isapprox(a.vel_norm, b.vel_norm; atol=0.0, rtol=VEL_RTOL) &&
           isapprox(a.mass_final, b.mass_final; atol=0.0, rtol=MASS_RTOL)
end

function metrics_delta(a::FlakeMetrics, b::FlakeMetrics)
    return (
        rows=(a.rows - b.rows),
        t_final=(a.t_final - b.t_final),
        pos_norm=(a.pos_norm - b.pos_norm),
        vel_norm=(a.vel_norm - b.vel_norm),
        mass_final=(a.mass_final - b.mass_final),
    )
end

failures = String[]

for case in CASES
    baseline = nothing
    for run_idx in 1:RUNS
        current = try
            run_case(case)
        catch err
            push!(failures, "case=$(case.name) run=$run_idx failed with $(typeof(err)): $(err)")
            break
        end

        if run_idx == 1
            baseline = current
            println("flake_case_baseline_ok case=$(case.name) rows=$(current.rows) tf=$(current.t_final) |r|=$(current.pos_norm) |v|=$(current.vel_norm) m=$(current.mass_final)")
            continue
        end

        metrics_consistent(baseline, current) || begin
            Δ = metrics_delta(current, baseline)
            push!(failures, "case=$(case.name) run=$run_idx inconsistent metrics: Δrows=$(Δ.rows), Δtf=$(Δ.t_final), Δ|r|=$(Δ.pos_norm), Δ|v|=$(Δ.vel_norm), Δm=$(Δ.mass_final)")
        end
    end
end

if !isempty(failures)
    println("flake_guard_failures=$(length(failures))")
    for msg in failures
        println("  - $msg")
    end
    error("Flake guard failed")
end

println("flake_guard_ok cases=$(length(CASES)) runs=$(RUNS)")
