include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "aerobraking_mission_plot_utils.jl"))

setup_gram_example!()

const EARTH_AEROBRAKING_PERIAPSIS_ALTITUDE_M = 125.0e3
const EARTH_AEROBRAKING_APOAPSIS_ALTITUDE_M = 60_500.0e3#20_000.0e3
const EARTH_AEROBRAKING_TARGET_ORBITS = 50
const EARTH_AEROBRAKING_TICK_FONT_SIZE = 12
const EARTH_AEROBRAKING_GUIDE_FONT_SIZE = 12
const EARTH_AEROBRAKING_LEGEND_FONT_SIZE = 12
const EARTH_AEROBRAKING_TITLE_FONT_SIZE = 15

planet = Earth("", SPICE_PATH)
initial_time = InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
ephemerides_model = SpiceEphemeridesModel()

ENV["SPACEAGORA_SOLVER_MODE"] = get(ENV, "SPACEAGORA_SOLVER_MODE", "split_imex")
ENV["SPACEAGORA_SPLIT_IMEX_SOLVER"] = get(ENV, "SPACEAGORA_SPLIT_IMEX_SOLVER", "kencarp4")
ENV["SPACEAGORA_VACUUM_GRAM_CACHE"] = get(ENV, "SPACEAGORA_VACUUM_GRAM_CACHE", "1")

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.05, 2.05, 2.8),
    panel_dims=(0.01, 5.7 / 2.0, 1.0),
    bus_mass=620.0,
    panel_mass_each=10.0,
    panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=InitialCondition(
        planet;
        ra=EARTH_AEROBRAKING_APOAPSIS_ALTITUDE_M,
        hp=EARTH_AEROBRAKING_PERIAPSIS_ALTITUDE_M,
        # i=90.0,
        # ω=90.0,
        i=45.0,
        ω=0.0,
        Ω=0.0,
        ν=180.0,
        initial_time=initial_time,
        ephemerides_model=ephemerides_model,
    ),
    prop_mass=200.0,
    id=1
)

earth_harmonics_file = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
sun_moon_gravity = NBodyGravityModel(body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet)
solar_radiation_pressure = SolarRadiationPressureModel(
    spacecraft.root.reflection_coefficient,
    spacecraft.root.ref_area
)

dynamic_effectors = (
    sun_moon_gravity,
    GravitationalHarmonicsModel(50, 50, earth_harmonics_file, planet),
    solar_radiation_pressure,
    AerodynamicCoefficientfM(),
)

function _save_earth_aerobraking_apoapsis_periapsis_plot(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    isfile(csv_path) || begin
        @warn "Skipping apoapsis/periapsis plot because simulation_results.csv was not found." csv_path
        return nothing
    end

    results = CSV.read(csv_path, DataFrame)
    if nrow(results) < 3
        @warn "Skipping apoapsis/periapsis plot because at least three saved samples are required." csv_path samples=nrow(results)
        return nothing
    end

    extrema = _derive_orbit_extrema_from_results(csv_path, args.environment_model.planet)
    peri_events = try
        _load_periapsis_events(results_dir)
    catch err
        if err isa ArgumentError
            @warn "Using sampled periapsis extrema for apoapsis/periapsis plot because no event-located table is available." reason=sprint(showerror, err)
            DataFrame()
        else
            rethrow()
        end
    end
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
        legend=:topright,
        tickfontsize=EARTH_AEROBRAKING_TICK_FONT_SIZE,
        guidefontsize=EARTH_AEROBRAKING_GUIDE_FONT_SIZE,
        legendfontsize=EARTH_AEROBRAKING_LEGEND_FONT_SIZE
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
        color=:crimson,
        tickfontsize=EARTH_AEROBRAKING_TICK_FONT_SIZE,
        guidefontsize=EARTH_AEROBRAKING_GUIDE_FONT_SIZE
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

    plot_path = joinpath(results_dir, "apoapsis_periapsis_vs_orbit.pdf")
    savefig(p, plot_path)
    println("Saved apoapsis/periapsis plot to $(abspath(plot_path))")
    return plot_path
end

function _save_earth_aerobraking_ground_track_plot(args::SimulationConfiguration)
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
        title="Earth Aerobraking Ground Track",
        tickfontsize=EARTH_AEROBRAKING_TICK_FONT_SIZE,
        guidefontsize=EARTH_AEROBRAKING_GUIDE_FONT_SIZE,
        titlefontsize=EARTH_AEROBRAKING_TITLE_FONT_SIZE,
        colorbar_tickfontsize=EARTH_AEROBRAKING_TICK_FONT_SIZE,
        colorbar_titlefontsize=EARTH_AEROBRAKING_GUIDE_FONT_SIZE
    )
    plot!(p, longitude_deg, latitude_deg; label=false, color=:gray35, linewidth=1.2, alpha=0.55)

    plot_path = joinpath(results_dir, "earth_ground_track.png")
    savefig(p, plot_path)
    println("Saved Earth ground track plot to $(abspath(plot_path))")
    return plot_path
end

function _save_earth_aerobraking_orbital_elements_plot(args::SimulationConfiguration)
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
                grid=true,
                tickfontsize=EARTH_AEROBRAKING_TICK_FONT_SIZE,
                guidefontsize=EARTH_AEROBRAKING_GUIDE_FONT_SIZE
            )
        )
    end
    xlabel!(subplots[5], "Mission Time (hr)")
    xlabel!(subplots[6], "Mission Time (hr)")

    p = plot(
        subplots...;
        layout=(3, 2),
        size=(1300, 1000),
        plot_title="Earth Aerobraking Orbital Elements",
        plot_titlefontsize=EARTH_AEROBRAKING_TITLE_FONT_SIZE
    )

    plot_path = joinpath(results_dir, "orbital_elements_3x2.png")
    savefig(p, plot_path)
    println("Saved orbital-elements 3x2 plot to $(abspath(plot_path))")
    return plot_path
end

function _try_save_earth_aerobraking_artifact(f::Function, description::String)
    try
        return f()
    catch err
        @warn "Skipping $(description)." reason=sprint(showerror, err)
        return nothing
    end
end

function _save_earth_aerobraking_fallback_periapsis_event_table_from_extrema(args::SimulationConfiguration)
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
        radial_velocity_mps=fill(NaN, length(extrema.peri.orbit))
    )
    event_path = joinpath(results_dir, "periapsis_events.csv")
    CSV.write(event_path, events)
    println("Saved fallback periapsis table from sampled extrema to $(abspath(event_path))")
    return event_path
end

function _save_earth_aerobraking_periapsis_event_table(sol, args::SimulationConfiguration, smoke_mode::Bool)
    if sol === nothing
        @warn "No solution object is available; deriving fallback periapsis table from saved samples."
        return _save_earth_aerobraking_fallback_periapsis_event_table_from_extrema(args)
    end
    return _save_periapsis_event_table(
        sol,
        args,
        args.environment_model.planet;
        allow_empty=smoke_mode
    )
end

function _save_earth_aerobraking_plot_suite(sol, args::SimulationConfiguration, smoke_mode::Bool)
    if !args.simulation_settings.results
        @warn "Skipping Earth aerobraking plots because simulation results are disabled." results_dir=args.simulation_settings.results_directory
        return nothing
    end

    csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
    if !isfile(csv_path)
        @warn "Skipping Earth aerobraking plots because no simulation results CSV was saved." csv_path=abspath(csv_path)
        return nothing
    end

    if _has_orbit_extrema_for_plots(args, args.environment_model.planet)
        _try_save_earth_aerobraking_artifact("periapsis event table") do
            _save_earth_aerobraking_periapsis_event_table(sol, args, smoke_mode)
        end
        _try_save_earth_aerobraking_artifact("apoapsis/periapsis altitude plot") do
            _save_earth_aerobraking_apoapsis_periapsis_plot(args)
        end
    end
    _try_save_earth_aerobraking_artifact("ground track plot") do
        _save_earth_aerobraking_ground_track_plot(args)
    end
    _try_save_earth_aerobraking_artifact("orbital elements plot") do
        _save_earth_aerobraking_orbital_elements_plot(args)
    end
    return nothing
end

orbital_period_s = 2π * sqrt(spacecraft.initial_condition.a^3 / planet.μ)
mission_time_cap_s = EARTH_AEROBRAKING_TARGET_ORBITS * orbital_period_s

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=mission_time_cap_s,
    initial_time=initial_time,
    dynamic_effectors=dynamic_effectors,
    density_model=GRAMAtmosphereModel(planet_name="earth"),
    ephemerides_model=ephemerides_model,
    orientation_sim=false,
    keplerian=false,
    EI_km=225.0,
    verbose=true,
    results_directory=joinpath(REPO_ROOT, "output", "earth_aerobraking")
)

args = SimulationConfiguration(
    file_paths=args.file_paths,
    simulation_settings=args.simulation_settings,
    mission_configuration=MissionConfiguration(
        mission_type=MissionOrbits,
        keplerian=args.mission_configuration.keplerian,
        number_of_orbits=EARTH_AEROBRAKING_TARGET_ORBITS,
        mission_time=mission_time_cap_s,
        orientation_sim=args.mission_configuration.orientation_sim,
        num_steps_to_save=args.mission_configuration.num_steps_to_save,
        data_rate=args.mission_configuration.data_rate
    ),
    environment_model=args.environment_model,
    dynamics_model=args.dynamics_model,
    guidance_model=args.guidance_model,
    navigation_model=args.navigation_model,
    control_model=args.control_model,
    initial_time=args.initial_time,
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-8,
        abstol_orbit=1e-8,
        dt_max_orbit=30.0,
        reltol_atmosphere=1e-8,
        abstol_atmosphere=1e-8,
        dt_max_atmosphere=5.0
    ),
    solver_config=SolverConfig(
        solver_mode=Symbol(get(ENV, "SPACEAGORA_SOLVER_MODE", "split_imex"))
    )
)

args_eff = SpaceAGORA.TelemetryVerification._example_smoke_args(args)
smoke_mode = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
print_thread_diagnostics(args_eff; label="AGORA_Earth_Aerobraking")
sol_ref = Ref{Any}(nothing)
sim_elapsed = @elapsed begin
    try
        sol_ref[] = run_simulation(args_eff; return_solution=true)
    catch err
        @error "Earth aerobraking solve failed; saving available results and plots before exiting." exception=(err, catch_backtrace())
        csv_path = joinpath(args_eff.simulation_settings.results_directory, "simulation_results.csv")
        if args_eff.simulation_settings.results && isfile(csv_path)
            df = CSV.read(csv_path, DataFrame)
            println("Saved $(nrow(df)) samples to $(abspath(csv_path))")
        end
        _save_earth_aerobraking_plot_suite(sol_ref[], args_eff, smoke_mode)
        rethrow()
    end
end
csv_path = joinpath(args_eff.simulation_settings.results_directory, "simulation_results.csv")
if args_eff.simulation_settings.results && isfile(csv_path)
    df = CSV.read(csv_path, DataFrame)
    println("Saved $(nrow(df)) samples to $(abspath(csv_path))")
end
println("COMPUTATIONAL TIME = $(sim_elapsed) s")
_save_earth_aerobraking_plot_suite(sol_ref[], args_eff, smoke_mode)
