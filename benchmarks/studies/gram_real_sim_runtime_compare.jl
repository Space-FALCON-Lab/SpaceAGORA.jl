const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using Dates
using DataFrames
using CSV
using SPICE
using Profile

include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
using .SimulationModel

const QUAT_MULT = SimulationModel.quat_mult
include(joinpath(REPO_ROOT, "src", "simulation", "execution", "run_simulation.jl"))

function _with_env(vars::Dict{String, String}, f::Function)
    old = Dict{String, Union{Nothing, String}}()
    for (k, v) in vars
        old[k] = haskey(ENV, k) ? ENV[k] : nothing
        ENV[k] = v
    end
    try
        return f()
    finally
        for (k, ov) in old
            if ov === nothing
                delete!(ENV, k)
            else
                ENV[k] = ov
            end
        end
    end
end

_with_env(f::Function, vars::Dict{String, String}) = _with_env(vars, f)

function _make_args(; mission_time_s::Float64)
    spice_path = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
    planet = Mars("", spice_path)

    ic = InitialCondition(
        ra=planet.Rp_e + 500_000.0,
        rp=planet.Rp_e + 100_000.0,
        i=25.0,
        ω=30.0,
        Ω=20.0,
        ν=-20.0
    )

    root = Link{0}(root=true, m=391.0, ref_area=2.2 * 1.7)
    spacecraft = SpacecraftModel(
        joints=Joint[],
        links=Link[root],
        root=root,
        instant_actuation=true,
        prop_mass=0.0,
        inertia_tensor=root.inertia,
        n_reaction_wheels=0,
        n_thrusters=0,
        initial_condition=ic,
        id=100
    )

    density_model = Base.invokelatest(
        GRAMAtmosphereModel;
        planet_name="mars",
        initial_time=InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=32.0)
    )

    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=mission_time_s,
            orientation_sim=false,
            num_steps_to_save=1000
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=125.0,
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=true
        ),
        dynamics_model=DynamicsModel([spacecraft], (InverseSquaredGravityModel(), AerodynamicCoefficientfM())),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=32.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=1.0,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_atmosphere=0.2
        )
    )
end

function _run_case(label::String, envvars::Dict{String, String}; mission_time_s::Float64, profile_run::Bool=false)
    args = _make_args(mission_time_s=mission_time_s)
    elapsed_s = 0.0
    n_steps = 0
    retcode = ""
    profile_file = ""
    stats = (
        density_calls=0,
        cache_enabled_calls=0,
        cache_hits=0,
        cache_misses=0,
        miss_time_window=0,
        miss_state_tolerance=0,
        direct_calls=0,
        refresh_calls=0,
        refresh_points_total=0,
        refresh_points_max=0,
        refresh_failures=0,
        refresh_elapsed_s=0.0,
        state_error_samples=0,
        alt_err_abs_max_m=0.0,
        lat_err_abs_max_deg=0.0,
        lon_err_abs_max_deg=0.0,
        alt_err_abs_sum_m=0.0,
        lat_err_abs_sum_deg=0.0,
        lon_err_abs_sum_deg=0.0
    )
    _with_env(envvars) do
        SimulationModel.SimulationCallbacks._gram_runtime_stats_reset!()
        if profile_run
            Profile.clear()
            elapsed_s = @elapsed begin
                sol = Profile.@profile run_simulation(args; isolate_state=false, return_solution=true)
                n_steps = length(sol.t)
                retcode = string(sol.retcode)
            end
            profile_file = joinpath(REPO_ROOT, "output", "gram_profile_$(label).txt")
            mkpath(dirname(profile_file))
            open(profile_file, "w") do io
                Profile.print(io; format=:flat, sortedby=:count, maxdepth=14, mincount=40)
            end
        else
            elapsed_s = @elapsed begin
                sol = run_simulation(args; isolate_state=false, return_solution=true)
                n_steps = length(sol.t)
                retcode = string(sol.retcode)
            end
        end
        stats = SimulationModel.SimulationCallbacks._gram_runtime_stats_snapshot()
    end
    return merge(
        (
            mode=label,
            elapsed_s=elapsed_s,
            n_steps=n_steps,
            retcode=retcode,
            profile_run=profile_run,
            profile_file=profile_file
        ),
        stats
    )
end

function run_benchmark()
    mission_time_s = 1_800.0

    off_env = Dict(
        "SPACEAGORA_GRAM_TRACK_CACHE" => "off",
        "SPACEAGORA_GRAM_PROFILE" => "1"
    )
    on_env = Dict(
        "SPACEAGORA_GRAM_TRACK_CACHE" => "on",
        "SPACEAGORA_GRAM_PROFILE" => "1"
    )
    on_unconstrained_env = Dict(
        "SPACEAGORA_GRAM_TRACK_CACHE" => "on",
        "SPACEAGORA_GRAM_PROFILE" => "1",
        "SPACEAGORA_GRAM_SEGMENT_CACHE_ENTRY_ALT_TOL_M" => "1e9",
        "SPACEAGORA_GRAM_SEGMENT_CACHE_ORBIT_ALT_TOL_M" => "1e9",
        "SPACEAGORA_GRAM_SEGMENT_CACHE_ENTRY_ANG_TOL_DEG" => "180.0",
        "SPACEAGORA_GRAM_SEGMENT_CACHE_ORBIT_ANG_TOL_DEG" => "180.0"
    )
    on_no_time_env = Dict(
        "SPACEAGORA_GRAM_TRACK_CACHE" => "on",
        "SPACEAGORA_GRAM_PROFILE" => "1",
        "SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW" => "1"
    )

    println("GRAM real-simulation runtime compare")
    println("mission_time = $(mission_time_s) s")
    println("Warmup...")
    _run_case("warmup_off", off_env; mission_time_s=60.0)
    _run_case("warmup_on", on_env; mission_time_s=60.0)

    rows = NamedTuple[]
    println("Measured run: point_to_point (cache=off)")
    push!(rows, _run_case("point_to_point", off_env; mission_time_s=mission_time_s))
    println("Measured run: cache_interpolation (cache=on)")
    push!(rows, _run_case("cache_interpolation", on_env; mission_time_s=mission_time_s))
    println("Measured run: cache_no_time_window (cache=on, ignore time window)")
    push!(rows, _run_case("cache_no_time_window", on_no_time_env; mission_time_s=mission_time_s))
    println("Measured run: cache_unconstrained (cache=on, huge alt/ang tol)")
    push!(rows, _run_case("cache_unconstrained", on_unconstrained_env; mission_time_s=mission_time_s))
    println("Profiled run: point_to_point (cache=off)")
    push!(rows, _run_case("point_to_point_profile", off_env; mission_time_s=mission_time_s, profile_run=true))
    println("Profiled run: cache_interpolation (cache=on)")
    push!(rows, _run_case("cache_interpolation_profile", on_env; mission_time_s=mission_time_s, profile_run=true))
    println("Profiled run: cache_no_time_window (cache=on, ignore time window)")
    push!(rows, _run_case("cache_no_time_window_profile", on_no_time_env; mission_time_s=mission_time_s, profile_run=true))
    println("Profiled run: cache_unconstrained (cache=on, huge alt/ang tol)")
    push!(rows, _run_case("cache_unconstrained_profile", on_unconstrained_env; mission_time_s=mission_time_s, profile_run=true))

    df = DataFrame(rows)
    mkpath(joinpath(REPO_ROOT, "output"))
    out = joinpath(REPO_ROOT, "output", "gram_real_sim_runtime_compare.csv")
    CSV.write(out, df)

    println("\nResults:")
    show(df, allrows=true, allcols=true)
    println("\nSaved CSV: $out")
    return df
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_benchmark()
end
