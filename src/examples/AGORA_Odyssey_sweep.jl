if !isdefined(@__MODULE__, :SimulationModel)
    include("../simulation_model/SimulationModel.jl")
end
using .SimulationModel
using Distributed
using SPICE
using StaticArrays
using LinearAlgebra
using CSV
using DataFrames

const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :run_simulation)
    include("../simulation/run_simulation.jl")
end
if !isdefined(@__MODULE__, :make_example_config)
    include("typed_example_utils.jl")
end

struct ConstantDensityModel <: AbstractDensityModel
    rho::Float64
    temp::Float64
end

function SimulationModel.EnvironmentModels.getDensity(
    model::ConstantDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)
    return model.rho, model.temp, SVector{3, Float64}(0.0, 0.0, 0.0)
end

struct TimedVelocityThrusterModel
    thrust::Float64
    start_time::Float64
    stop_time::Float64
end

function SimulationModel.calcControlForceTorque(
    model::TimedVelocityThrusterModel,
    u,
    p::ODEParams,
    i::Int64,
    t::Float64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    if t < model.start_time || t > model.stop_time
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    vel = SVector{3, Float64}(u.vel)
    speed = norm(vel)
    if speed <= 1e-9
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    force = model.thrust * vel / speed
    return force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcControlEffect!(
    model::TimedVelocityThrusterModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    return nothing
end

const ODYSSEY_SWEEP_ROOT = joinpath(REPO_ROOT, "output", "odyssey_sweep")
const ODYSSEY_SINGLE_ROOT = joinpath(REPO_ROOT, "output", "odyssey_single")
const ODYSSEY_WORKER_CONTEXT = Ref{Any}(nothing)
const ODYSSEY_SWEEP_CONFIG = (
    root_directory=ODYSSEY_SWEEP_ROOT,
    max_parallel=min(Sys.CPU_THREADS, 8),
    ra_m=[28_559.615e3],
    rp_altitude_m=[87_000.0],
    inclination_deg=[93.522, 95.0, 100.0, 110.0, 85.0, 80.0],
    aop_deg=[109.7454, 90.0, 45.0, 0.0, 315.0, 270.0],
    raan_deg=[28.1517, 0.0, 90.0, 180.0, 270.0],
    ta_deg=[175.0],
    thrust_n=[4.0],
    thrust_start_s=[900.0],
    thrust_stop_s=[1_200.0]
)

@inline _smoke_mode() = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
@inline function _active_project_arg()::String
    active = Base.active_project()
    active === nothing && return REPO_ROOT
    return dirname(String(active))
end

@inline function _case_slug(case_id::Int)::String
    return "case_" * lpad(string(case_id), 4, '0')
end

function _build_odyssey_worker_context()
    planet = Mars("", SPICE_PATH)
    smoke_mode = _smoke_mode()
    mars_harmonics_file = joinpath(REPO_ROOT, "Gravity_harmonics_data", "Mars50c.csv")
    dynamic_effectors = smoke_mode ? (InverseSquaredGravityModel(),) : (
        InverseSquaredGravityModel(),
        GravitationalHarmonicsModel(20, 20, mars_harmonics_file, planet),
        AerodynamicCoefficientfM()
    )
    density_model = GRAMAtmosphereModel(planet_name="mars")
    return (
        planet=planet,
        smoke_mode=smoke_mode,
        dynamic_effectors=dynamic_effectors,
        density_model=density_model
    )
end

function _odyssey_worker_context()
    if ODYSSEY_WORKER_CONTEXT[] === nothing
        ODYSSEY_WORKER_CONTEXT[] = _build_odyssey_worker_context()
    end
    return ODYSSEY_WORKER_CONTEXT[]
end

function build_odyssey_args(;
    ra_m::Float64=28_559.615e3,
    rp_altitude_m::Float64=87_000.0,
    inclination_deg::Float64=93.522,
    aop_deg::Float64=109.7454,
    raan_deg::Float64=28.1517,
    ta_deg::Float64=175.0,
    thrust_n::Float64=4.0,
    thrust_start_s::Float64=900.0,
    thrust_stop_s::Float64=1_200.0,
    results_directory::String=ODYSSEY_SINGLE_ROOT,
    verbose::Bool=false,
    context=_odyssey_worker_context()
)::SimulationConfiguration
    planet = context.planet
    smoke_mode = context.smoke_mode

    ic = InitialCondition(
        ra=ra_m,
        rp=planet.Rp_e + rp_altitude_m,
        i=inclination_deg,
        ω=aop_deg,
        Ω=raan_deg,
        ν=ta_deg
    )

    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.2, 2.6, 1.7),
        panel_dims=(0.01, 3.89 / 2.0, 1.7),
        bus_mass=391.0,
        panel_mass_each=10.0,
        panel_offset_y=2.6 / 2.0 + 3.89 / 4.0,
        ic=ic,
        prop_mass=50.0,
        id=100
    )

    base_args = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=3_600.0 * 100.0,
        initial_time=InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=32.0),
        dynamic_effectors=context.dynamic_effectors,
        density_model=context.density_model,
        orientation_sim=false,
        keplerian=smoke_mode,
        EI_km=125.0,
        verbose=verbose
    )

    thruster = TimedVelocityThrusterModel(thrust_n, thrust_start_s, thrust_stop_s)
    sim_settings = SimulationSettings(
        results=base_args.simulation_settings.results,
        verbose=verbose,
        results_directory=results_directory,
        generate_plots=base_args.simulation_settings.generate_plots,
        generate_filenames=base_args.simulation_settings.generate_filenames,
        normalize=base_args.simulation_settings.normalize,
        save_csv=base_args.simulation_settings.save_csv,
        checkpoint_enabled=base_args.simulation_settings.checkpoint_enabled,
        checkpoint_interval_s=base_args.simulation_settings.checkpoint_interval_s,
        checkpoint_directory=base_args.simulation_settings.checkpoint_directory,
        resume_from_checkpoint=base_args.simulation_settings.resume_from_checkpoint
    )

    return SimulationConfiguration(
        file_paths=base_args.file_paths,
        simulation_settings=sim_settings,
        mission_configuration=base_args.mission_configuration,
        environment_model=base_args.environment_model,
        dynamics_model=base_args.dynamics_model,
        guidance_model=base_args.guidance_model,
        navigation_model=base_args.navigation_model,
        control_model=ControlModel(control_effectors=(thruster,), control_rates=[30.0]),
        initial_time=base_args.initial_time,
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=20.0,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_atmosphere=0.2
        )
    )
end

function run_odyssey_case(;
    case_id::Int=1,
    ra_m::Float64=28_559.615e3,
    rp_altitude_m::Float64=87_000.0,
    inclination_deg::Float64=93.522,
    aop_deg::Float64=109.7454,
    raan_deg::Float64=28.1517,
    ta_deg::Float64=175.0,
    thrust_n::Float64=4.0,
    thrust_start_s::Float64=900.0,
    thrust_stop_s::Float64=1_200.0,
    root_directory::String=ODYSSEY_SINGLE_ROOT,
    verbose::Bool=false,
    context=_odyssey_worker_context()
)::NamedTuple
    case_dir = joinpath(root_directory, _case_slug(case_id))
    mkpath(case_dir)

    args = build_odyssey_args(
        ra_m=ra_m,
        rp_altitude_m=rp_altitude_m,
        inclination_deg=inclination_deg,
        aop_deg=aop_deg,
        raan_deg=raan_deg,
        ta_deg=ta_deg,
        thrust_n=thrust_n,
        thrust_start_s=thrust_start_s,
        thrust_stop_s=thrust_stop_s,
        results_directory=case_dir,
        verbose=verbose,
        context=context
    )

    csv_path = joinpath(case_dir, "simulation_results.csv")
    runtime_s = @elapsed begin
        cd(case_dir) do
            run_simulation(args)
        end
    end

    samples = 0
    if args.simulation_settings.results && isfile(csv_path)
        df = CSV.read(csv_path, DataFrame)
        samples = nrow(df)
        println("[$(_case_slug(case_id))] Saved $samples samples to $(abspath(csv_path))")
    end
    println("[$(_case_slug(case_id))] COMPUTATIONAL TIME = $(runtime_s) s")

    return (
        case_id=case_id,
        case_dir=case_dir,
        ra_m=ra_m,
        rp_altitude_m=rp_altitude_m,
        inclination_deg=inclination_deg,
        aop_deg=aop_deg,
        raan_deg=raan_deg,
        ta_deg=ta_deg,
        thrust_n=thrust_n,
        thrust_start_s=thrust_start_s,
        thrust_stop_s=thrust_stop_s,
        runtime_s=runtime_s,
        samples=samples,
        csv_path=csv_path
    )
end

function _sweep_cases_from_env()
    cfg = ODYSSEY_SWEEP_CONFIG
    cases = NamedTuple[]
    case_id = 1
    for ra_m in cfg.ra_m,
        rp_altitude_m in cfg.rp_altitude_m,
        inclination_deg in cfg.inclination_deg,
        aop_deg in cfg.aop_deg,
        raan_deg in cfg.raan_deg,
        ta_deg in cfg.ta_deg,
        thrust_n in cfg.thrust_n,
        thrust_start_s in cfg.thrust_start_s,
        thrust_stop_s in cfg.thrust_stop_s
        push!(cases, (
            case_id=case_id,
            ra_m=ra_m,
            rp_altitude_m=rp_altitude_m,
            inclination_deg=inclination_deg,
            aop_deg=aop_deg,
            raan_deg=raan_deg,
            ta_deg=ta_deg,
            thrust_n=thrust_n,
            thrust_start_s=thrust_start_s,
            thrust_stop_s=thrust_stop_s,
            root_directory=cfg.root_directory
        ))
        case_id += 1
    end
    return cases
end

function write_sweep_summary(results::Vector{<:NamedTuple}, root_directory::String)::String
    summary_path = joinpath(root_directory, "sweep_summary.csv")
    mkpath(root_directory)
    CSV.write(summary_path, DataFrame(results))
    return summary_path
end

function _ensure_odyssey_workers(count::Int)
    count <= 0 && return Int[]
    current = workers()
    missing = count - length(current)
    if missing > 0
        addprocs(missing; exeflags="--project=$(_active_project_arg())")
    end
    current = workers()
    isempty(current) && error("No distributed workers are available for the Odyssey sweep.")
    selected_count = min(count, length(current))
    selected = current[1:selected_count]
    selected_count == count || @warn "Requested $count Odyssey workers, but only $selected_count are available. Continuing with the available workers."
    script_path = abspath(@__FILE__)
    for pid in selected
        remotecall_wait(pid) do
            include(script_path)
            nothing
        end
    end
    return selected
end

function run_odyssey_sweep(; max_parallel::Int=ODYSSEY_SWEEP_CONFIG.max_parallel)
    cases = _sweep_cases_from_env()
    isempty(cases) && error("Odyssey sweep produced zero cases.")

    root_directory = first(cases).root_directory
    mkpath(root_directory)
    max_parallel = clamp(max_parallel, 1, length(cases))

    println("Running $(length(cases)) Odyssey cases with max_parallel=$max_parallel")
    results = NamedTuple[]
    walltime_s = @elapsed begin
        initial_workers = Set(workers())
        pool_workers = Int[]
        try
            pool_workers = _ensure_odyssey_workers(max_parallel)
            wp = WorkerPool(pool_workers)
            results_or_errors = pmap(wp, cases) do case
                try
                    run_odyssey_case(
                        case_id=case.case_id,
                        ra_m=case.ra_m,
                        rp_altitude_m=case.rp_altitude_m,
                        inclination_deg=case.inclination_deg,
                        aop_deg=case.aop_deg,
                        raan_deg=case.raan_deg,
                        ta_deg=case.ta_deg,
                        thrust_n=case.thrust_n,
                        thrust_start_s=case.thrust_start_s,
                        thrust_stop_s=case.thrust_stop_s,
                        root_directory=case.root_directory
                    )
                catch err
                    case.case_id => sprint(showerror, err)
                end
            end

            errors = Pair{Int, Any}[]
            empty!(results)
            for item in results_or_errors
                if item isa Pair
                    push!(errors, item)
                else
                    push!(results, item)
                end
            end

            isempty(errors) || error("Odyssey sweep failed for case ids $(join(first.(errors), ", ")).")
            sort!(results; by=r -> r.case_id)

            summary_path = write_sweep_summary(results, root_directory)
            println("Wrote sweep summary to $(abspath(summary_path))")
        finally
            added_workers = [pid for pid in pool_workers if !(pid in initial_workers)]
            isempty(added_workers) || rmprocs(added_workers)
        end
    end

    println("ODYSSEY SWEEP WALLTIME = $(walltime_s) s")
    return results
end
if abspath(PROGRAM_FILE) == @__FILE__
    cases = _sweep_cases_from_env()
    if length(cases) > 1
        run_odyssey_sweep()
    else
        only_case = first(cases)
        run_odyssey_case(
            case_id=only_case.case_id,
            ra_m=only_case.ra_m,
            rp_altitude_m=only_case.rp_altitude_m,
            inclination_deg=only_case.inclination_deg,
            aop_deg=only_case.aop_deg,
            raan_deg=only_case.raan_deg,
            ta_deg=only_case.ta_deg,
            thrust_n=only_case.thrust_n,
            thrust_start_s=only_case.thrust_start_s,
            thrust_stop_s=only_case.thrust_stop_s
        )
    end
end
