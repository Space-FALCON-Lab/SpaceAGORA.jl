include(joinpath(@__DIR__, "common.jl"))
using CSV
using DataFrames
using LinearAlgebra
using PlotlyJS: Plot, Layout, attr, scatter, scatter3d, surface
using SPICE
using StaticArrays

function inbound_hyperbolic_ic(
    planet;
    v_inf_kms::Float64,
    periapsis_altitude_km::Float64,
    entry_interface_altitude_km::Float64,
    inclination_deg::Float64,
    arg_periapsis_deg::Float64,
    raan_deg::Float64
)
    v_inf = v_inf_kms * 1e3
    rp = planet.Rp_e + periapsis_altitude_km * 1e3
    r_entry = planet.Rp_e + entry_interface_altitude_km * 1e3

    a = -planet.μ / v_inf^2
    e = 1.0 + rp * v_inf^2 / planet.μ
    p = a * (1.0 - e^2)
    ν_entry = -acos((p / r_entry - 1.0) / e)

    return InitialCondition(
        a,
        e,
        inclination_deg,
        arg_periapsis_deg,
        raan_deg,
        rad2deg(ν_entry)
    )
end

function summarize_capture(csv_path::String, planet)
    df = CSV.read(csv_path, DataFrame)
    isempty(df) && error("No saved trajectory samples found at $(abspath(csv_path)).")

    first_row = df[1, :]
    last_row = df[end, :]

    r0 = SVector{3, Float64}(first_row.sc1_pos_1, first_row.sc1_pos_2, first_row.sc1_pos_3)
    v0 = SVector{3, Float64}(first_row.sc1_vel_1, first_row.sc1_vel_2, first_row.sc1_vel_3)
    rf = SVector{3, Float64}(last_row.sc1_pos_1, last_row.sc1_pos_2, last_row.sc1_pos_3)
    vf = SVector{3, Float64}(last_row.sc1_vel_1, last_row.sc1_vel_2, last_row.sc1_vel_3)

    ε0 = dot(v0, v0) / 2.0 - planet.μ / norm(r0)
    εf = dot(vf, vf) / 2.0 - planet.μ / norm(rf)
    oe_f = SpaceAGORA.TelemetryVerification.rvtoorbitalelement(rf, vf, planet)

    min_altitude_km = minimum(df.sc1_altitude) * 1e-3
    max_heat_rate = maximum(df.sc1_heat_rate)
    captured = εf < 0.0

    println("Initial specific orbital energy = $(round(ε0 / 1e6, digits=6)) MJ/kg")
    println("Final specific orbital energy   = $(round(εf / 1e6, digits=6)) MJ/kg")
    println("Minimum altitude                = $(round(min_altitude_km, digits=3)) km")
    println("Peak heat rate                  = $(round(max_heat_rate, digits=3)) W/m^2")
    println("Captured                        = $(captured)")

    if captured
        apoapsis_altitude_km = oe_f[1] * (1.0 + oe_f[2]) * 1e-3 - planet.Rp_e * 1e-3
        periapsis_altitude_km = oe_f[1] * (1.0 - oe_f[2]) * 1e-3 - planet.Rp_e * 1e-3
        println("Final apoapsis altitude         = $(round(apoapsis_altitude_km, digits=3)) km")
        println("Final periapsis altitude        = $(round(periapsis_altitude_km, digits=3)) km")
        println("Final eccentricity              = $(round(oe_f[2], digits=6))")
    end
end

function save_trajectory_plot(csv_path::String, planet; output_path::Union{Nothing, String}=nothing)
    df = CSV.read(csv_path, DataFrame)
    isempty(df) && error("No saved trajectory samples found at $(abspath(csv_path)).")

    pos_x_km = Float64.(df.sc1_pos_1) ./ 1e3
    pos_y_km = Float64.(df.sc1_pos_2) ./ 1e3
    pos_z_km = Float64.(df.sc1_pos_3) ./ 1e3
    altitude_km = Float64.(df.sc1_altitude) ./ 1e3

    θ = range(0.0, 2pi; length=72)
    ϕ = range(0.0, pi; length=36)
    earth_radius_km = planet.Rp_e / 1e3
    earth_x = [earth_radius_km * cos(θj) * sin(ϕi) for ϕi in ϕ, θj in θ]
    earth_y = [earth_radius_km * sin(θj) * sin(ϕi) for ϕi in ϕ, θj in θ]
    earth_z = [earth_radius_km * cos(ϕi) for ϕi in ϕ, θj in θ]

    traces = [
        surface(
            x=earth_x,
            y=earth_y,
            z=earth_z,
            name="Earth",
            opacity=0.45,
            showscale=false,
            colorscale=[[0.0, "#2c7fb8"], [1.0, "#7fcdbb"]],
            hoverinfo="skip"
        ),
        scatter3d(
            x=pos_x_km,
            y=pos_y_km,
            z=pos_z_km,
            mode="lines",
            name="Trajectory",
            line=attr(color="#d95f02", width=6),
            text=["Altitude: $(round(h, digits=2)) km" for h in altitude_km],
            hovertemplate="x: %{x:.2f} km<br>y: %{y:.2f} km<br>z: %{z:.2f} km<br>%{text}<extra></extra>"
        ),
        scatter3d(
            x=[pos_x_km[1]],
            y=[pos_y_km[1]],
            z=[pos_z_km[1]],
            mode="markers",
            name="Start",
            marker=attr(color="#1b9e77", size=5),
            hovertemplate="Start<br>x: %{x:.2f} km<br>y: %{y:.2f} km<br>z: %{z:.2f} km<extra></extra>"
        ),
        scatter3d(
            x=[pos_x_km[end]],
            y=[pos_y_km[end]],
            z=[pos_z_km[end]],
            mode="markers",
            name="End",
            marker=attr(color="#7570b3", size=5),
            hovertemplate="End<br>x: %{x:.2f} km<br>y: %{y:.2f} km<br>z: %{z:.2f} km<extra></extra>"
        )
    ]

    layout = Layout(
        title="Earth Aerocapture Trajectory",
        scene=attr(
            xaxis=attr(title="x (km)", showgrid=true),
            yaxis=attr(title="y (km)", showgrid=true),
            zaxis=attr(title="z (km)", showgrid=true),
            aspectmode="data"
        ),
        legend=attr(x=0.02, y=0.98),
        margin=attr(l=0, r=0, b=0, t=40)
    )
    plot = Plot(traces, layout)

    plot_path = isnothing(output_path) ? joinpath(dirname(csv_path), "trajectory_3d.html") : output_path
    open(plot_path, "w") do io
        show(io, MIME("text/html"), plot; include_plotlyjs="cdn", full_html=true)
    end
    println("Saved interactive trajectory plot to $(abspath(plot_path))")
    return plot_path
end

function save_specific_energy_plot(csv_path::String, planet; output_path::Union{Nothing, String}=nothing)
    df = CSV.read(csv_path, DataFrame)
    isempty(df) && error("No saved trajectory samples found at $(abspath(csv_path)).")

    times_hr = Float64.(df.time) ./ 3600.0
    pos_x = Float64.(df.sc1_pos_1)
    pos_y = Float64.(df.sc1_pos_2)
    pos_z = Float64.(df.sc1_pos_3)
    vel_x = Float64.(df.sc1_vel_1)
    vel_y = Float64.(df.sc1_vel_2)
    vel_z = Float64.(df.sc1_vel_3)

    specific_energy_mj_kg = similar(times_hr)
    for idx in eachindex(times_hr)
        r = SVector{3, Float64}(pos_x[idx], pos_y[idx], pos_z[idx])
        v = SVector{3, Float64}(vel_x[idx], vel_y[idx], vel_z[idx])
        specific_energy_mj_kg[idx] = (dot(v, v) / 2.0 - planet.μ / norm(r)) / 1e6
    end

    capture_idx = findfirst(<(0.0), specific_energy_mj_kg)
    traces = [
        scatter(
            x=times_hr,
            y=specific_energy_mj_kg,
            mode="lines",
            name="Specific orbital energy",
            line=attr(color="#d95f02", width=3),
            hovertemplate="Time: %{x:.3f} hr<br>Energy: %{y:.6f} MJ/kg<extra></extra>"
        )
    ]
    if !isnothing(capture_idx)
        push!(traces, scatter(
            x=[times_hr[capture_idx]],
            y=[specific_energy_mj_kg[capture_idx]],
            mode="markers",
            name="Capture threshold crossing",
            marker=attr(color="#1b9e77", size=9),
            hovertemplate="Capture crossing<br>Time: %{x:.3f} hr<br>Energy: %{y:.6f} MJ/kg<extra></extra>"
        ))
    end

    layout = Layout(
        title="Earth Aerocapture Specific Orbital Energy",
        xaxis=attr(title="Elapsed time (hr)"),
        yaxis=attr(title="Specific orbital energy (MJ/kg)", zeroline=true, zerolinecolor="#444"),
        legend=attr(x=0.02, y=0.98),
        margin=attr(l=70, r=30, b=60, t=40)
    )
    plot = Plot(traces, layout)

    plot_path = isnothing(output_path) ? joinpath(dirname(csv_path), "specific_energy.html") : output_path
    open(plot_path, "w") do io
        show(io, MIME("text/html"), plot; include_plotlyjs="cdn", full_html=true)
    end
    println("Saved interactive specific energy plot to $(abspath(plot_path))")
    return plot_path
end

function open_plot_in_browser(plot_path::String)
    try
        run(`xdg-open $plot_path`; wait=false)
    catch err
        @warn "Failed to open interactive plot automatically." plot_path exception=(err, catch_backtrace())
    end
    return nothing
end

function run_with_progress(args_eff::SimulationConfiguration, smoke_mode::Bool)
    if smoke_mode || haskey(ENV, "SPACEAGORA_PROGRESS_INTERVAL_S")
        return @elapsed run_simulation(args_eff)
    end

    default_progress_interval_s = 600.0
    println("Progress logging enabled every $(Int(default_progress_interval_s)) simulated seconds.")
    return withenv("SPACEAGORA_PROGRESS_INTERVAL_S" => string(default_progress_interval_s)) do
        @elapsed run_simulation(args_eff)
    end
end

planet = Earth("", SPICE_PATH)
smoke_mode = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
initial_time = InitialTime(year=2025, month=1, day=15, hour=12, minute=0, second=0.0)

ic = inbound_hyperbolic_ic(
    planet;
    v_inf_kms=smoke_mode ? 0.85 : 1.05,
    periapsis_altitude_km=smoke_mode ? 72.0 : 68.0,
    entry_interface_altitude_km=160.0,
    inclination_deg=28.5,
    arg_periapsis_deg=25.0,
    raan_deg=40.0
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(3.2, 2.7, 2.1),
    panel_dims=(0.02, 2.4, 1.4),
    bus_mass=2_000.0,
    panel_mass_each=45.0,
    panel_offset_y=2.15,
    ic=ic,
    reflection_coefficient=0.9,
    prop_mass=0.0
)

earth_harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
dynamic_effectors = smoke_mode ? (
    InverseSquaredGravityModel(),
    AerodynamicCoefficientfM()
) : (
    InverseSquaredGravityModel(),
    GravitationalHarmonicsModel(8, 8, earth_harmonics_file, planet),
    AerodynamicCoefficientfM()
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=smoke_mode ? 2.0 * 3600.0 : 10.0 * 3600.0,
    initial_time=initial_time,
    dynamic_effectors=dynamic_effectors,
    density_model=GRAMAtmosphereModel(planet_name="earth", initial_time=initial_time),
    orientation_sim=false,
    keplerian=false,
    EI_km=300.0,
    verbose=!smoke_mode
)

args_eff = SpaceAGORA.TelemetryVerification._example_smoke_args(args)
if args_eff.mission_configuration.mission_time != args.mission_configuration.mission_time
    println(
        "Smoke mode active: mission time reduced from $(args.mission_configuration.mission_time / 3600.0) hr " *
        "to $(args_eff.mission_configuration.mission_time / 3600.0) hr."
    )
end
sim_elapsed = run_with_progress(args_eff, smoke_mode)
csv_path = joinpath(args_eff.simulation_settings.results_directory, "simulation_results.csv")

if args_eff.simulation_settings.results && isfile(csv_path)
    df = CSV.read(csv_path, DataFrame)
    println("Saved $(nrow(df)) samples to $(abspath(csv_path))")
    summarize_capture(csv_path, planet)
    trajectory_plot_path = save_trajectory_plot(csv_path, planet)
    specific_energy_plot_path = save_specific_energy_plot(csv_path, planet)
    open_plot_in_browser(trajectory_plot_path)
    open_plot_in_browser(specific_energy_plot_path)
end

println("COMPUTATIONAL TIME = $(sim_elapsed) s")
