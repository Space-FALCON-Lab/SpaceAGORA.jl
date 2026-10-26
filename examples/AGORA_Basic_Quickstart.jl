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

function save_quickstart_plots(csv_path; planet_radius_km=nothing)
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
            label=false,
            xlabel="Time since start (hr)",
            ylabel="Altitude (km)",
            title="Instantaneous geodetic altitude",
            lw=2
        )
        # The periapsis altitude is a few hundred km while the altitude spans tens
        # of thousands, so it gets its own panel and scale.
        p_peri = plot(
            t_hr,
            periapsis_altitude_km;
            label=false,
            xlabel="Time since start (hr)",
            ylabel="Altitude (km)",
            title="Osculating periapsis altitude (spherical)",
            lw=2
        )

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
            plot(p_alt, p_peri, p_speed; layout=(3, 1), size=(900, 900), left_margin=5Plots.mm),
            plots_dir,
            "quickstart_altitude_speed.png"
        )
    end

    if _has_columns(df, (:sc1_pos_1, :sc1_pos_2, :sc1_pos_3))
        x_km = df[!, :sc1_pos_1] ./ 1e3
        y_km = df[!, :sc1_pos_2] ./ 1e3
        z_km = df[!, :sc1_pos_3] ./ 1e3
        # Two projections share one planet-centred scale: a polar orbit seen in a
        # single x-y projection collapses to a nearly straight line.
        halfspan = 1.05 * maximum(maximum(abs, v) for v in (x_km, y_km, z_km))
        θ = range(0.0, 2π; length=181)
        function projection(a, b, a_label, b_label, title; legend)
            p = plot(
                a,
                b;
                label="Orbit",
                xlabel=a_label,
                ylabel=b_label,
                title=title,
                aspect_ratio=:equal,
                xlims=(-halfspan, halfspan),
                ylims=(-halfspan, halfspan),
                legend=legend,
                lw=2
            )
            if planet_radius_km !== nothing
                plot!(p, planet_radius_km .* cos.(θ), planet_radius_km .* sin.(θ);
                      label="Planet", color=:gray, lw=1.5)
            end
            scatter!(p, [first(a)], [first(b)]; label="Start", ms=4)
            scatter!(p, [last(a)], [last(b)]; label="End", ms=4)
            return p
        end
        p = plot(
            projection(x_km, y_km, "Inertial x (km)", "Inertial y (km)", "Inertial trajectory, x-y"; legend=:bottomright),
            projection(x_km, z_km, "Inertial x (km)", "Inertial z (km)", "Inertial trajectory, x-z"; legend=false);
            layout=(1, 2),
            size=(1100, 540),
            margin=5Plots.mm
        )
        _save_quickstart_plot!(saved_paths, p, plots_dir, "quickstart_inertial_trajectory.png")
    end

    if _has_columns(df, (:sc1_pos_1, :sc1_pos_2, :sc1_pos_3))
        x_km = df[!, :sc1_pos_1] ./ 1e3
        y_km = df[!, :sc1_pos_2] ./ 1e3
        z_km = df[!, :sc1_pos_3] ./ 1e3
        # Use the same span on every axis so a thin orbit remains readable.
        bounds = extrema.((x_km, y_km, z_km))
        halfspan = max(maximum(hi - lo for (lo, hi) in bounds) * 0.55, 1.0)
        limits = [((lo + hi) / 2 - halfspan, (lo + hi) / 2 + halfspan)
                  for (lo, hi) in bounds]
        tick_values = [round.([(lo + hi) / 2 - 0.8 * halfspan,
                               (lo + hi) / 2,
                               (lo + hi) / 2 + 0.8 * halfspan]; sigdigits=2)
                       for (lo, hi) in limits]
        axis_ticks = [(values, string.(values)) for values in tick_values]
        p = plot(
            x_km,
            y_km,
            z_km;
            label="Orbit",
            xlabel="x (km)",
            ylabel="y (km)",
            zlabel="z (km)",
            title="3D Inertial Orbit",
            aspect_ratio=:equal,
            camera=(35, 25),
            lw=2,
            size=(900, 700),
            xlims=limits[1],
            ylims=limits[2],
            zlims=limits[3],
            xticks=axis_ticks[1],
            yticks=axis_ticks[2],
            zticks=axis_ticks[3],
            legend=:outertopright,
            margin=5Plots.mm,
            tickfontsize=9,
            guidefontsize=11
        )
        scatter!(p, [first(x_km)], [first(y_km)], [first(z_km)]; label="Start", ms=4)
        scatter!(p, [last(x_km)], [last(y_km)], [last(z_km)]; label="End", ms=4)
        _save_quickstart_plot!(saved_paths, p, plots_dir, "quickstart_3d_orbit.png")
    end

    if _has_columns(df, (:sc1_longitude_deg, :sc1_latitude_deg))
        longitude = Float64.(df[!, :sc1_longitude_deg])
        latitude = Float64.(df[!, :sc1_latitude_deg])
        track_lon = Float64[]
        track_lat = Float64[]
        for k in eachindex(longitude)
            if k > 1 && abs(longitude[k] - longitude[k - 1]) > 180.0
                push!(track_lon, NaN)
                push!(track_lat, NaN)
            end
            push!(track_lon, longitude[k])
            push!(track_lat, latitude[k])
        end
        p = plot(
            track_lon,
            track_lat;
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

# Parallel execution is one flag in the solver setup: with
# `SolverConfig(parallel=true)` SpaceAGORA decides how to use the threads Julia
# was started with (`julia --threads=auto ...`). It is off here by default so
# the quickstart stays a plain serial run; set SPACEAGORA_EXAMPLE_PARALLEL=1 to
# try it. See docs/src/user/parallel_execution.md.
const QUICKSTART_PARALLEL = get(ENV, "SPACEAGORA_EXAMPLE_PARALLEL", "0") == "1"

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
    verbose=true,
    solver_config=QUICKSTART_PARALLEL ? SolverConfig(parallel=true) : nothing
)

csv_path = run_and_report(args)
save_quickstart_plots(csv_path; planet_radius_km=planet.Rp_e / 1e3)
