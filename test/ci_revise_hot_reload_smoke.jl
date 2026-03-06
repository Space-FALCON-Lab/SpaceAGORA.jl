const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

using CSV
using DataFrames
using SPICE
using Revise

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel

# Revise-based hot-swap workflows require callback invokelatest dispatch.
ENV["SPACEAGORA_DEV_HOT_RELOAD"] = "1"

# SimulationEngine uses SimulationModel and provides the canonical run_simulation entrypoint.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :SimulationEngine)
    include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SimulationEngine.run_simulation
end

const spice_path = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const planet = Earth("", spice_path)

struct ReviseHotReloadEffector <: SimulationModel.AbstractControlEffectorModel end

function write_effector_methods(path::String; thrust_mag::Float64, mass_rate::Float64)
    src = """
    using LinearAlgebra
    using StaticArrays

    function SimulationModel.calcControlForceTorque(
        ::ReviseHotReloadEffector,
        u::AbstractVector{Float64},
        p::SimulationModel.ODEParams,
        i::Int64,
        t::Float64
    )
        v = SVector{3, Float64}(u.vel)
        vm = norm(v)
        if vm <= 1e-12 || !isfinite(vm)
            return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
        end
        force = $(thrust_mag) * (v / vm)
        return force, SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    function SimulationModel.calcControlMassFlowRate(
        ::ReviseHotReloadEffector,
        u::AbstractVector{Float64},
        p::SimulationModel.ODEParams,
        i::Int64,
        t::Float64
    )
        return $(mass_rate)
    end

    function SimulationModel.calcControlEffect!(
        ::ReviseHotReloadEffector,
        u,
        p::SimulationModel.ODEParams,
        t::Float64,
        i::Int64
    )
        return nothing
    end
    """
    write(path, src)
end

function run_once()::DataFrame
    root = Link{0}(root=true, m=120.0, ref_area=1.0)
    ic = InitialCondition(planet.Rp_e + 500e3, 0.0, 0.0, 0.0, 0.0, 0.0)
    spacecraft = SpacecraftModel(
        Joint[],
        [root],
        root,
        true,
        root.m,
        5.0,
        root.inertia,
        0,
        0,
        ic,
        1
    )

    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true,
            verbose=false,
            generate_plots=false,
            results_directory="output",
            save_csv=true,
            normalize=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=40.0,
            orientation_sim=false,
            num_steps_to_save=100
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
        control_model=ControlModel(control_effectors=(ReviseHotReloadEffector(),), control_rates=[1.0]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=2.0
        )
    )

    mktempdir() do tmp
        cd(tmp) do
            run_simulation(args; isolate_state=true)
            csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            if !isfile(csv_path)
                error("Expected simulation_results.csv was not written")
            end
            return CSV.read(csv_path, DataFrame)
        end
    end
end

mktempdir() do tmp
    methods_path = joinpath(tmp, "revise_hot_reload_methods.jl")

    write_effector_methods(methods_path; thrust_mag=0.2, mass_rate=-0.01)
    Revise.includet(methods_path)
    df1 = run_once()
    m1 = Float64(df1.sc1_mass[end])

    # Ensure file timestamp changes before revising.
    sleep(0.2)
    write_effector_methods(methods_path; thrust_mag=0.2, mass_rate=-0.05)
    Revise.revise()

    df2 = run_once()
    m2 = Float64(df2.sc1_mass[end])

    if !(m2 < m1 - 1e-3)
        error("Revise hot-reload smoke failed: expected lower final mass after method update (m1=$(m1), m2=$(m2)).")
    end
end

println("revise_hot_reload_smoke_ok")
