include(joinpath(@__DIR__, "common.jl"))

using CSV
using DataFrames
using Plots

# This example is the smallest repo-only onboarding path:
# - no GRAM
# - no SPICE
# - no licensed external assets
# - simple J2 orbital dynamics

function _has_columns(df::DataFrame, cols::Tuple{Vararg{Symbol}})
    names = propertynames(df)
    return all(col -> col in names, cols)
end

function _time_hours(df::DataFrame)
    t = Float64.(df[!, :time])
    return (t .- first(t)) ./ 3600.0
end

function _save_quickstart_plot!(saved_paths::Vector{String}, p, plots_dir::String, name::String)
    path = joinpath(plots_dir, name)
    savefig(p, path)
    push!(saved_paths, path)
    return nothing
end

function save_quickstart_plots(csv_path)
    csv_path === nothing && return String[]
    isfile(csv_path) || return String[]

    df = CSV.read(csv_path, DataFrame)
    nrow(df) == 0 && return String[]

    plots_dir = joinpath(dirname(csv_path), "plots")
    mkpath(plots_dir)

    t_hr = _time_hours(df)
    saved_paths = String[]

    if _has_columns(df, (:sc1_altitude, :sc1_periapsis_altitude, :sc1_vel_1, :sc1_vel_2, :sc1_vel_3))
        altitude_km = df[!, :sc1_altitude] ./ 1e3
        periapsis_altitude_km = df[!, :sc1_periapsis_altitude] ./ 1e3
        speed_km_s = sqrt.(df[!, :sc1_vel_1].^2 .+ df[!, :sc1_vel_2].^2 .+ df[!, :sc1_vel_3].^2) ./ 1e3

        p_alt = plot(
            t_hr,
            altitude_km;
            label="Altitude",
            xlabel="Time since start (hr)",
            ylabel="Altitude (km)",
            title="Altitude",
            lw=2
        )
        plot!(p_alt, t_hr, periapsis_altitude_km; label="Periapsis altitude", lw=2, ls=:dash)

        p_speed = plot(
            t_hr,
            speed_km_s;
            label=false,
            xlabel="Time since start (hr)",
            ylabel="Speed (km/s)",
            title="Inertial Speed",
            lw=2
        )

        _save_quickstart_plot!(
            saved_paths,
            plot(p_alt, p_speed; layout=(2, 1), size=(900, 650)),
            plots_dir,
            "quickstart_altitude_speed.png"
        )
    end

    if _has_columns(df, (:sc1_pos_1, :sc1_pos_2))
        p = plot(
            df[!, :sc1_pos_1] ./ 1e3,
            df[!, :sc1_pos_2] ./ 1e3;
            label=false,
            xlabel="Inertial x (km)",
            ylabel="Inertial y (km)",
            title="Inertial Trajectory",
            aspect_ratio=:equal,
            lw=2
        )
        scatter!(p, [df[1, :sc1_pos_1] / 1e3], [df[1, :sc1_pos_2] / 1e3]; label="Start", ms=4)
        scatter!(p, [df[end, :sc1_pos_1] / 1e3], [df[end, :sc1_pos_2] / 1e3]; label="End", ms=4)
        _save_quickstart_plot!(saved_paths, p, plots_dir, "quickstart_inertial_trajectory.png")
    end

    if _has_columns(df, (:sc1_pos_1, :sc1_pos_2, :sc1_pos_3))
        x_km = df[!, :sc1_pos_1] ./ 1e3
        y_km = df[!, :sc1_pos_2] ./ 1e3
        z_km = df[!, :sc1_pos_3] ./ 1e3
        p = plot(
            x_km,
            y_km,
            z_km;
            label="Orbit",
            xlabel="Inertial x (km)",
            ylabel="Inertial y (km)",
            zlabel="Inertial z (km)",
            title="3D Inertial Orbit",
            aspect_ratio=:equal,
            camera=(35, 25),
            lw=2
        )
        scatter!(p, [first(x_km)], [first(y_km)], [first(z_km)]; label="Start", ms=4)
        scatter!(p, [last(x_km)], [last(y_km)], [last(z_km)]; label="End", ms=4)
        _save_quickstart_plot!(saved_paths, p, plots_dir, "quickstart_3d_orbit.png")
    end

    if _has_columns(df, (:sc1_longitude_deg, :sc1_latitude_deg))
        p = plot(
            df[!, :sc1_longitude_deg],
            df[!, :sc1_latitude_deg];
            label=false,
            xlabel="Longitude (deg)",
            ylabel="Latitude (deg)",
            title="Ground Track",
            xlim=(-180, 180),
            ylim=(-90, 90),
            lw=2
        )
        _save_quickstart_plot!(saved_paths, p, plots_dir, "quickstart_ground_track.png")
    end

    if !isempty(saved_paths)
        println("Saved quickstart plots:")
        foreach(path -> println("  ", abspath(path)), saved_paths)
    end

    return saved_paths
end

planet = make_no_gram_planet(:earth)

# Build a simple bus + two-panel spacecraft from the shared example helper.
spacecraft = make_three_body_spacecraft(
    bus_dims=(2.05, 2.05, 2.8),
    panel_dims=(0.01, 5.7 / 2.0, 1.0),
    bus_mass=620.0,
    panel_mass_each=10.0,
    panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=InitialCondition(
        ra=56_378.7978559e3,
        rp=planet.Rp_e + 200_590.0,
        i=89.876,
        ω=75.505,
        Ω=104.115,
        ν=175.0
    ),
    prop_mass=200.0,
    id=1
)

# `make_example_config` supplies the standard example boilerplate.
# We keep the environment intentionally simple here so installation is easy to
# validate before moving to GRAM-backed workflows.
args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=12.0 * 3600.0,
    initial_time=InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredJ2GravityModel(),),
    density_model=NoAtmosphereModel(),
    ephemerides_model=SimpleEphemeridesModel(),
    orientation_sim=false,
    keplerian=true,
    EI_km=120.0,
    verbose=true
)

csv_path = run_and_report(args)
save_quickstart_plots(csv_path)
