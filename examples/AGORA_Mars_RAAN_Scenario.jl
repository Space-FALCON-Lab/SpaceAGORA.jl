include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "aerobraking_mission_plot_utils.jl"))

using Printf
using DiffEqBase: ContinuousCallback, terminate!

setup_gram_example!()

const MARS_RAAN_TERMINATE_ON_PERIAPSIS_RAISE = true
const MARS_RAAN_INITIAL_PERIAPSIS_ALTITUDE_M = 117.3e3
const MARS_RAAN_PERIAPSIS_RAISE_LIMIT_M = 1.0e3
const MARS_RAAN_TARGET_APOAPSIS_RADIUS_M = 9_000.0e3

const EMPTY_MARS_RAAN_MANEUVER_SCHEDULE = (
    maneuver_orbit_number=Int[],
    maneuver_Δv=Float64[]
)

function _mars_raan_values_deg()
    raw = strip(get(ENV, "SPACEAGORA_MARS_RAAN_VALUES", "90"))
    isempty(raw) && return [90.0]

    values = Float64[]
    for token in split(raw, ',')
        value = tryparse(Float64, strip(token))
        value === nothing && throw(ArgumentError(
            "SPACEAGORA_MARS_RAAN_VALUES must be a comma-separated list of degrees; got $(repr(raw))."
        ))
        push!(values, value)
    end
    return values
end

function _results_directory_for_raan(base_dir::String, raan_deg::Float64, n_raan::Int)
    n_raan == 1 && return base_dir
    return joinpath(base_dir, "mars_raan_" * replace(@sprintf("%.3f", raan_deg), "." => "p"))
end

function _make_mars_raan_scenario(raan_deg::Float64, results_directory::String)
    planet = Mars("", SPICE_PATH)
    smoke_mode = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
    initial_time = InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0)
    ephemerides_model = SpiceEphemeridesModel()

    ENV["SPACEAGORA_SOLVER_MODE"] = get(ENV, "SPACEAGORA_SOLVER_MODE", "split_imex")
    ENV["SPACEAGORA_SPLIT_IMEX_SOLVER"] = get(ENV, "SPACEAGORA_SPLIT_IMEX_SOLVER", "kencarp4")
    ENV["SPACEAGORA_VACUUM_GRAM_CACHE"] = get(ENV, "SPACEAGORA_VACUUM_GRAM_CACHE", "1")

    ic = InitialCondition(
        planet;
        ra=18_477.0e3,  # Mars radius + apoapsis altitude
        hp=MARS_RAAN_INITIAL_PERIAPSIS_ALTITUDE_M,
        i=63.43,
        ω=220.0,
        Ω=raan_deg,
        initial_time=initial_time,
        ephemerides_model=ephemerides_model
    )

    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.2, 2.6, 1.7),
        panel_dims=(0.01, 5.5 / 1.35, 2.6),
        bus_mass=391.0,
        panel_mass_each=10.0,
        panel_offset_y=2.6 / 2.0 + 5.5 / 4.0,
        ic=ic,
        reflection_coefficient=0.9,
        prop_mass=50.0,
        id=110
    )

    mars_harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "Mars50c.csv")
    sun_gravity = NBodyGravityModel(body_names=("Sun",), primary_body_name="Mars", planet=planet)
    solar_radiation_pressure = SolarRadiationPressureModel(
        spacecraft.root.reflection_coefficient,
        spacecraft.root.ref_area
    )
    dynamic_effectors = smoke_mode ? (
        InverseSquaredGravityModel(),
        sun_gravity
    ) : (
        sun_gravity,
        GravitationalHarmonicsModel(2, 0, mars_harmonics_file, planet),
        # solar_radiation_pressure,
        AerodynamicCoefficientfM()
    )

    thruster = BaseThrusterModel(
        thrust=[4.0],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[-1.0],
        stop_burn_time=[-1.0],
        Isp=[300.0]
    )
    guidance_effector = ApoapsisTargetPeriapsisRaiseGuidanceModel(
        target_apoapsis_radius_m=MARS_RAAN_TARGET_APOAPSIS_RADIUS_M,
        target_periapsis_altitude_m=200.0e3
    )

    base_args = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=3_600.0 * 24.0 * 142.0,
        initial_time=initial_time,
        dynamic_effectors=dynamic_effectors,
        density_model=Base.invokelatest(GRAMAtmosphereModel; planet_name="mars"),
        ephemerides_model=ephemerides_model,
        orientation_sim=false,
        keplerian=smoke_mode,
        EI_km=160.0,
        results_directory=results_directory
    )

    return SimulationConfiguration(
        file_paths=base_args.file_paths,
        simulation_settings=base_args.simulation_settings,
        mission_configuration=base_args.mission_configuration,
        environment_model=base_args.environment_model,
        dynamics_model=base_args.dynamics_model,
        guidance_model=GuidanceModel(guidance_effectors=(guidance_effector,), guidance_rates=[30.0]),
        navigation_model=base_args.navigation_model,
        control_model=ControlModel(control_effectors=(thruster,), control_rates=[10.0]),
        initial_time=base_args.initial_time,
        integration_tolerances=base_args.integration_tolerances,
        solver_config=base_args.solver_config
    )
end

function _mars_raan_oblate_altitude_m(integrator, sat_idx::Int=1)::Float64
    pos = SVector{3, Float64}(integrator.u.sc[sat_idx].pos)
    vel = SVector{3, Float64}(integrator.u.sc[sat_idx].vel)
    args = integrator.p.args
    planet = args.environment_model.planet
    ephemerides_model = args.environment_model.ephemerides_model
    et = integrator.p.shared_buffers.et_start[] + Float64(integrator.t)
    pos_p, _ = SimulationModel.SimulationCallbacks.r_intor_p!(pos, vel, planet, et, ephemerides_model)
    return SimulationModel.SimulationCallbacks.rtolatlong(pos_p, planet, ephemerides_model)[1]
end

function _mars_raan_periapsis_raise_termination_callback()
    limit_m = MARS_RAAN_INITIAL_PERIAPSIS_ALTITUDE_M + MARS_RAAN_PERIAPSIS_RAISE_LIMIT_M

    condition(u, t, integrator) = begin
        pos = SVector{3, Float64}(u.sc[1].pos)
        vel = SVector{3, Float64}(u.sc[1].vel)
        return dot(pos, vel)
    end

    function affect!(integrator)
        periapsis_altitude_m = _mars_raan_oblate_altitude_m(integrator)
        if periapsis_altitude_m > limit_m
            println(@sprintf(
                "Stopping Mars RAAN scenario at t = %.3f s: periapsis altitude %.3f km is more than %.3f km above initial %.3f km.",
                Float64(integrator.t),
                periapsis_altitude_m * 1e-3,
                MARS_RAAN_PERIAPSIS_RAISE_LIMIT_M * 1e-3,
                MARS_RAAN_INITIAL_PERIAPSIS_ALTITUDE_M * 1e-3
            ))
            terminate!(integrator)
        end
        return nothing
    end

    return ContinuousCallback(condition, affect!, nothing)
end

function _mars_raan_extra_callbacks()
    return MARS_RAAN_TERMINATE_ON_PERIAPSIS_RAISE ?
        (_mars_raan_periapsis_raise_termination_callback(),) :
        ()
end

function _save_mars_raan_ground_track_plot(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    isfile(csv_path) || throw(ArgumentError("Simulation results CSV not found at $(abspath(csv_path))."))

    df = CSV.read(csv_path, DataFrame)
    longitude_deg = _require_float_column(df, :sc1_longitude_deg)
    latitude_deg = _require_float_column(df, :sc1_latitude_deg)
    altitude_km = _require_float_column(df, :sc1_altitude) ./ 1e3

    p = scatter(
        longitude_deg,
        latitude_deg;
        zcolor=altitude_km,
        colorbar_title="Altitude (km)",
        xlabel="Longitude (deg)",
        ylabel="Latitude (deg)",
        label=false,
        markerstrokewidth=0,
        markersize=3,
        color=:viridis,
        xlims=(-180, 180),
        ylims=(-90, 90),
        grid=true,
        size=(1100, 600),
        title="Mars Ground Track"
    )
    plot!(p, longitude_deg, latitude_deg; label=false, color=:gray35, linewidth=1.2, alpha=0.55)

    plot_path = joinpath(results_dir, "mars_ground_track.png")
    savefig(p, plot_path)
    println("Saved Mars ground track plot to $(abspath(plot_path))")
    return plot_path
end

function _save_mars_raan_apoapsis_periapsis_plot(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    extrema = _derive_orbit_extrema_from_results(csv_path, args.environment_model.planet)
    peri_events = _load_periapsis_events(results_dir)
    peri_orbit = nrow(peri_events) > 0 ? peri_events.orbit : extrema.peri.orbit
    peri_altitude_km = nrow(peri_events) > 0 ? peri_events.altitude_km : extrema.peri.altitude_km

    p = plot(
        peri_orbit,
        peri_altitude_km;
        xlabel="Orbit Number",
        ylabel="Periapsis Altitude (km)",
        label="Periapsis",
        linewidth=2.5,
        marker=:circle,
        color=:dodgerblue,
        grid=true,
        legend=:topright
    )
    apo_axis = twinx()
    plot!(
        apo_axis,
        extrema.apo.orbit,
        extrema.apo.altitude_km;
        ylabel="Apoapsis Altitude (km)",
        label=false,
        linewidth=2.5,
        marker=:diamond,
        color=:crimson
    )
    plot!(
        p,
        [NaN],
        [NaN];
        label="Apoapsis",
        linewidth=2.5,
        marker=:diamond,
        color=:crimson
    )

    plot_path = joinpath(results_dir, "apoapsis_periapsis_vs_orbit.png")
    savefig(p, plot_path)
    println("Saved apoapsis/periapsis plot to $(abspath(plot_path))")
    return plot_path
end

function _save_mars_raan_periapsis_latlon_plot(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    peri_events = _load_periapsis_events(results_dir)
    if nrow(peri_events) == 0
        @warn "Skipping periapsis latitude/longitude plot because no event-located periapsis crossings were found." results_dir
        return nothing
    end

    p = plot(
        peri_events.orbit,
        peri_events.latitude_deg;
        xlabel="Orbit Number",
        ylabel="Periapsis Latitude (deg)",
        label="Latitude",
        linewidth=2.5,
        marker=:circle,
        color=:royalblue,
        grid=true,
        legend=:topright
    )
    lon_axis = twinx()
    plot!(
        lon_axis,
        peri_events.orbit,
        peri_events.longitude_deg;
        ylabel="Periapsis Longitude (deg)",
        label=false,
        linewidth=2.5,
        marker=:diamond,
        color=:firebrick
    )
    plot!(
        p,
        [NaN],
        [NaN];
        label="Longitude",
        linewidth=2.5,
        marker=:diamond,
        color=:firebrick
    )

    plot_path = joinpath(results_dir, "periapsis_latitude_longitude_vs_orbit.png")
    savefig(p, plot_path)
    println("Saved periapsis latitude/longitude plot to $(abspath(plot_path))")
    return plot_path
end

function _save_mars_raan_orbital_elements_plot(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    time_s, x_km, y_km, z_km = _simulation_position_samples(args)
    time_v_s, vx_kms, vy_kms, vz_kms = _simulation_velocity_samples(args)
    length(time_s) == length(time_v_s) && all(time_s .== time_v_s) ||
        throw(ArgumentError("Position and velocity samples do not share the same time grid."))

    oe = _orbital_elements_from_state_samples(
        x_km,
        y_km,
        z_km,
        vx_kms,
        vy_kms,
        vz_kms,
        args.environment_model.planet
    )
    time_hr = time_s ./ 3600.0
    specs = (
        (:semi_major_axis_km, "Semi-major Axis (km)", :royalblue),
        (:eccentricity, "Eccentricity", :firebrick),
        (:inclination_deg, "Inclination (deg)", :seagreen),
        (:raan_deg, "RAAN (deg)", :darkorange2),
        (:arg_periapsis_deg, "Argument of Periapsis (deg)", :purple),
        (:true_anomaly_deg, "True Anomaly (deg)", :teal)
    )

    subplots = Plots.Plot[]
    for (field, ylabel, color) in specs
        push!(
            subplots,
            plot(
                time_hr,
                getproperty(oe, field);
                ylabel=ylabel,
                label=false,
                linewidth=2.0,
                color=color,
                grid=true
            )
        )
    end
    xlabel!(subplots[5], "Mission Time (hr)")
    xlabel!(subplots[6], "Mission Time (hr)")

    p = plot(
        subplots...;
        layout=(3, 2),
        size=(1300, 1000),
        plot_title="Mars RAAN Scenario Orbital Elements"
    )

    plot_path = joinpath(results_dir, "orbital_elements_3x2.png")
    savefig(p, plot_path)
    println("Saved orbital-elements 3x2 plot to $(abspath(plot_path))")
    return plot_path
end

function _try_save_mars_raan_artifact(f::Function, description::String)
    try
        return f()
    catch err
        @warn "Skipping $(description)." reason=sprint(showerror, err)
        return nothing
    end
end

function _save_fallback_periapsis_event_table_from_extrema(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    extrema = _derive_orbit_extrema_from_results(csv_path, args.environment_model.planet)
    events = DataFrame(
        orbit=extrema.peri.orbit,
        time_s=fill(NaN, length(extrema.peri.orbit)),
        altitude_km=extrema.peri.altitude_km,
        latitude_deg=extrema.peri.latitude_deg,
        longitude_deg=extrema.peri.longitude_deg,
        radius_km=fill(NaN, length(extrema.peri.orbit)),
        radial_velocity_kms=fill(NaN, length(extrema.peri.orbit))
    )
    event_path = joinpath(results_dir, "periapsis_events.csv")
    CSV.write(event_path, events)
    println("Saved fallback periapsis table from sampled extrema to $(abspath(event_path))")
    return event_path
end

function _save_mars_raan_plot_suite(sol, args::SimulationConfiguration, smoke_mode::Bool)
    if !args.simulation_settings.results
        @warn "Skipping Mars RAAN plots because simulation results are disabled." results_dir=args.simulation_settings.results_directory
        return nothing
    end

    csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
    if !isfile(csv_path)
        @warn "Skipping Mars RAAN plots because no simulation results CSV was saved." csv_path=abspath(csv_path)
        return nothing
    end

    if sol === nothing
        @warn "No solution object is available after solve failure; deriving fallback periapsis table from saved samples."
        _try_save_mars_raan_artifact("fallback periapsis table") do
            _save_fallback_periapsis_event_table_from_extrema(args)
        end
    else
        _try_save_mars_raan_artifact("event-located periapsis table") do
            _save_periapsis_event_table(
                sol,
                args,
                args.environment_model.planet;
                allow_empty=smoke_mode
            )
        end
    end

    if _has_orbit_extrema_for_plots(args, args.environment_model.planet)
        _try_save_mars_raan_artifact("apoapsis/periapsis altitude plot") do
            _save_mars_raan_apoapsis_periapsis_plot(args)
        end
        _try_save_mars_raan_artifact("apoapsis radius plot") do
            _save_apoapsis_radius_plot(
                args,
                args.environment_model.planet,
                EMPTY_MARS_RAAN_MANEUVER_SCHEDULE
            )
        end
        _try_save_mars_raan_artifact("periapsis latitude/longitude plot") do
            _save_mars_raan_periapsis_latlon_plot(args)
        end
    end

    _try_save_mars_raan_artifact("ground track plot") do
        _save_mars_raan_ground_track_plot(args)
    end
    _try_save_mars_raan_artifact("orbital elements plot") do
        _save_mars_raan_orbital_elements_plot(args)
    end
    return nothing
end

raan_values_deg = _mars_raan_values_deg()
base_results_directory = joinpath(REPO_ROOT, "output", "mars_raan_scenario")

for raan_deg in raan_values_deg
    println("Running Mars RAAN scenario with RAAN = $(raan_deg) deg")
    args = _make_mars_raan_scenario(
        raan_deg,
        _results_directory_for_raan(base_results_directory, raan_deg, length(raan_values_deg))
    )
    args_eff = SpaceAGORA.TelemetryVerification._example_smoke_args(args)
    print_thread_diagnostics(args_eff; label="AGORA_Mars_RAAN_Scenario RAAN=$(raan_deg) deg")
    sol = nothing
    sim_elapsed = @elapsed begin
        try
            sol = run_simulation(args_eff; return_solution=true, extra_callbacks=_mars_raan_extra_callbacks())
        catch err
            @error "Mars RAAN scenario solve failed; saving available results and plots before exiting." raan_deg exception=(err, catch_backtrace())
            csv_path = joinpath(args_eff.simulation_settings.results_directory, "simulation_results.csv")
            if args_eff.simulation_settings.results && isfile(csv_path)
                df = CSV.read(csv_path, DataFrame)
                println("Saved $(nrow(df)) samples to $(abspath(csv_path))")
            end
            _save_mars_raan_plot_suite(nothing, args_eff, get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1")
            rethrow()
        end
    end
    csv_path = joinpath(args_eff.simulation_settings.results_directory, "simulation_results.csv")
    if args_eff.simulation_settings.results && isfile(csv_path)
        df = CSV.read(csv_path, DataFrame)
        println("Saved $(nrow(df)) samples to $(abspath(csv_path))")
    end
    println("COMPUTATIONAL TIME = $(sim_elapsed) s")
    _save_mars_raan_plot_suite(sol, args_eff, get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1")
end
