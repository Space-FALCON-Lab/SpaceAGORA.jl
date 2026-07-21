include(joinpath(@__DIR__, "common.jl"))
using CSV
using DataFrames
using LinearAlgebra
using Plots
import PlotlyJS
using StaticArrays

planet = Mars()

# ABTS Odyssey control-test limits and orbit seed.
limit_qdot = 0.15
limit_q = 30.0
limit_dyn_pressure = 0.5
ra = 28_559.615e3
rp = planet.Rp_e + 77e3
orbit_period_s = 2pi * sqrt(((ra + rp) / 2)^3 / planet.μ)

ic = InitialCondition(
    ra=ra,
    rp=rp,
    i=93.6,
    ω=109.7454,
    Ω=28.1517,
    ν=180.0,
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.2, 2.6, 1.7),
    panel_dims=(0.01, 3.76 / 2.0, 1.93),
    bus_mass=391.0,
    panel_mass_each=10.0,
    panel_offset_y=2.6 / 2.0 + 3.76 / 4.0,
    ic=ic,
    reflection_coefficient=0.9,
    prop_mass=50.0,
    id=101,
)

energy_state = AerobrakingEnergyDepletionState(num_sats=1)
energy_config = AerobrakingEnergyDepletionConfig(
    guidance_modes=(:targeting, :max_energy_depletion),
    max_energy_submodes=(:heat_rate, :heat_load),
    heat_load_switch_solver=:tpbvp_integration,
    controlled_panel_links=(2, 3),
    target_apoapsis_radius_m=26_750e3,
    max_alpha_rad=pi / 2,
    min_alpha_rad=1e-4,
    heat_rate_limit_w_cm2=limit_qdot,
    heat_load_limit_j_cm2=limit_q,
    structural_load_limit_pa=limit_dyn_pressure,
)

for link in spacecraft.links
    link.α = energy_config.min_alpha_rad
end

guidance = AerobrakingEnergyDepletionGuidanceModel(energy_config, energy_state)
control = AerobrakingEnergyDepletionControlModel(energy_config, energy_state)

base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=orbit_period_s,
    initial_time=InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
    density_model=ExponentialAtmosphereModel(planet),
    ephemerides_model=SimpleEphemeridesModel(),
    orientation_sim=false,
    keplerian=true,
    EI_km=160.0,
)

args = SimulationConfiguration(
    file_paths=base_args.file_paths,
    simulation_settings=base_args.simulation_settings,
    mission_configuration=MissionConfiguration(
        mission_type=base_args.mission_configuration.mission_type,
        keplerian=base_args.mission_configuration.keplerian,
        number_of_orbits=base_args.mission_configuration.number_of_orbits,
        mission_time=base_args.mission_configuration.mission_time,
        orientation_sim=base_args.mission_configuration.orientation_sim,
        num_steps_to_save=base_args.mission_configuration.num_steps_to_save,
        data_rate=1.0,
    ),
    environment_model=base_args.environment_model,
    dynamics_model=base_args.dynamics_model,
    guidance_model=GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[3.0]),
    navigation_model=base_args.navigation_model,
    control_model=ControlModel(control_effectors=(control,), control_rates=[0.1]),
    initial_time=base_args.initial_time,
    integration_tolerances=base_args.integration_tolerances,
)

panel_alpha_save_field = SaveField(
    :panel_alpha_rad,
    (u, t, integrator) -> [
        SVector{2, Float64}(
            integrator.p.args.dynamics_model.spacecraft[i].links[2].α,
            integrator.p.args.dynamics_model.spacecraft[i].links[3].α,
        )
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="panel_alpha_rad",
)

function _edg_runtime_control(integrator)
    return integrator.p.args.control_model.control_effectors[1]
end

function _edg_mode_code(mode::Symbol)::Float64
    mode == :targeting && return 2.0
    mode == :max_energy_depletion && return 1.0
    mode == :safe_low_drag && return -1.0
    return 0.0
end

function _edg_mode_name(code::Real)
    value = round(Int, Float64(code))
    value == 2 && return :targeting
    value == 1 && return :max_energy_depletion
    value == -1 && return :safe_low_drag
    return :inactive
end

mode_save_field = SaveField(
    :edg_mode_code,
    (u, t, integrator) -> [
        _edg_mode_code(_edg_runtime_control(integrator).state.selected_mode[i])
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_mode_code",
)

targeting_active_save_field = SaveField(
    :edg_targeting_active,
    (u, t, integrator) -> [
        _edg_runtime_control(integrator).state.targeting_active[i] ? 1.0 : 0.0
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_targeting_active",
)

targeting_switch_save_field = SaveField(
    :edg_targeting_switch_s,
    (u, t, integrator) -> [
        _edg_runtime_control(integrator).state.targeting_switch_s[i]
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_targeting_switch_s",
)

target_energy_save_field = SaveField(
    :edg_target_energy_jkg,
    (u, t, integrator) -> [
        _edg_runtime_control(integrator).state.target_energy_jkg[i]
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_target_energy_jkg",
)

bracket_min_save_field = SaveField(
    :edg_bracket_min_energy_jkg,
    (u, t, integrator) -> [
        _edg_runtime_control(integrator).state.bracket_min_energy_jkg[i]
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_bracket_min_energy_jkg",
)

bracket_max_save_field = SaveField(
    :edg_bracket_max_energy_jkg,
    (u, t, integrator) -> [
        _edg_runtime_control(integrator).state.bracket_max_energy_jkg[i]
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_bracket_max_energy_jkg",
)

save_fields = vcat(
    default_save_fields(args),
    [
        panel_alpha_save_field,
        mode_save_field,
        targeting_active_save_field,
        targeting_switch_save_field,
        target_energy_save_field,
        bracket_min_save_field,
        bracket_max_save_field,
    ],
)

function _series_magnitude(x, y, z)
    return sqrt.(x .^ 2 .+ y .^ 2 .+ z .^ 2)
end

function _orbit_metrics_from_row(row, planet)
    pos = SVector{3, Float64}(row.sc1_pos_1, row.sc1_pos_2, row.sc1_pos_3)
    vel = SVector{3, Float64}(row.sc1_vel_1, row.sc1_vel_2, row.sc1_vel_3)
    energy = dot(vel, vel) / 2 - planet.μ / norm(pos)
    oe = SimulationModel.ControlHooks.rvtoorbitalelement(pos, vel, planet)
    a, e = oe[1], oe[2]
    return (
        energy_jkg=energy,
        apoapsis_radius_m=a * (1.0 + e),
        periapsis_radius_m=a * (1.0 - e),
    )
end

function _drag_passage_indices(df::DataFrame)
    indices = findall(df.sc1_altitude .<= 160_000.0)
    isempty(indices) && error("No drag-passage samples found at or below 160 km.")
    return indices
end

function _remove_stale_energy_depletion_plots!(plot_dir::AbstractString)
    keep = Set([
        "drag_passage_summary.png",
        "drag_passage_trajectory_summary.png",
        "spacecraft_orbit_3d_mars_apo_to_apo.png",
    ])
    isdir(plot_dir) || return nothing
    for name in readdir(plot_dir)
        if endswith(name, ".png") && !(name in keep)
            rm(joinpath(plot_dir, name); force=true)
        end
    end
    return nothing
end

function _save_energy_depletion_plots(csv_path::AbstractString, plot_dir::AbstractString, planet)
    df = CSV.read(csv_path, DataFrame)
    mkpath(plot_dir)
    _remove_stale_energy_depletion_plots!(plot_dir)
    gr()

    plot_theme = (
        linewidth=2.6,
        grid=true,
        gridalpha=0.25,
        framestyle=:box,
        foreground_color_legend=nothing,
        background_color_legend=:white,
        guidefontsize=12,
        tickfontsize=10,
        legendfontsize=9,
        titlefontsize=13,
        left_margin=12Plots.mm,
        right_margin=8Plots.mm,
        top_margin=7Plots.mm,
        bottom_margin=11Plots.mm,
        dpi=180,
    )

    drag_indices = _drag_passage_indices(df)
    t_drag = df.time[drag_indices] .- first(df.time[drag_indices])

    drag = _series_magnitude(df.sc1_drag_1, df.sc1_drag_2, df.sc1_drag_3)
    velocity = _series_magnitude(df.sc1_vel_1, df.sc1_vel_2, df.sc1_vel_3)
    position = _series_magnitude(df.sc1_pos_1, df.sc1_pos_2, df.sc1_pos_3)
    radial_velocity = (
        df.sc1_pos_1 .* df.sc1_vel_1 .+
        df.sc1_pos_2 .* df.sc1_vel_2 .+
        df.sc1_pos_3 .* df.sc1_vel_3
    ) ./ max.(position .* velocity, eps(Float64))
    flight_path_angle_deg = rad2deg.(asin.(clamp.(radial_velocity, -1.0, 1.0)))

    heat_rate_plot = plot(
        t_drag,
        df.sc1_heat_rate[drag_indices];
        xlabel="Time from entry [s]",
        ylabel="Heat rate [W/cm^2]",
        label="Heat rate",
        title="Heat Rate",
        color=:royalblue,
        plot_theme...,
    )
    hline!(heat_rate_plot, [limit_qdot]; label="Limit", linestyle=:dash, color=:crimson, linewidth=1.8)

    heat_load_plot = plot(
        t_drag,
        df.sc1_heat_load[drag_indices];
        xlabel="Time from entry [s]",
        ylabel="Heat load [J/cm^2]",
        label="Heat load",
        title="Heat Load",
        color=:darkorange3,
        plot_theme...,
    )
    hline!(heat_load_plot, [limit_q]; label="Limit", linestyle=:dash, color=:crimson, linewidth=1.8)

    drag_plot = plot(
        t_drag,
        drag[drag_indices];
        xlabel="Time from entry [s]",
        ylabel="Drag [N]",
        label="Drag",
        title="Drag",
        color=:seagreen4,
        plot_theme...,
    )

    panel_alpha_plot = plot(
        t_drag,
        rad2deg.(df.sc1_panel_alpha_rad_1[drag_indices]);
        xlabel="Time from entry [s]",
        ylabel="Panel AoA [deg]",
        label="Panel link 2",
        title="Solar Panel AoA",
        color=:purple4,
        plot_theme...,
    )
    plot!(
        panel_alpha_plot,
        t_drag,
        rad2deg.(df.sc1_panel_alpha_rad_2[drag_indices]);
        label="Panel link 3",
        linewidth=2.3,
        linestyle=:dash,
        color=:deeppink3,
    )

    control_summary = plot(
        heat_rate_plot,
        heat_load_plot,
        drag_plot,
        panel_alpha_plot;
        layout=(2, 2),
        size=(1500, 1000),
        margin=8Plots.mm,
    )
    savefig(control_summary, joinpath(plot_dir, "drag_passage_summary.png"))

    altitude_plot = plot(
        t_drag,
        df.sc1_altitude[drag_indices] ./ 1_000.0;
        xlabel="Time from entry [s]",
        ylabel="Altitude [km]",
        label="Altitude",
        title="Altitude",
        color=:royalblue,
        plot_theme...,
    )
    velocity_plot = plot(
        t_drag,
        velocity[drag_indices] ./ 1_000.0;
        xlabel="Time from entry [s]",
        ylabel="Velocity [km/s]",
        label="Velocity",
        title="Velocity",
        color=:seagreen4,
        plot_theme...,
    )
    flight_path_angle_plot = plot(
        t_drag,
        flight_path_angle_deg[drag_indices];
        xlabel="Time from entry [s]",
        ylabel="Flight path angle [deg]",
        label="FPA",
        title="Flight Path Angle",
        color=:darkorange3,
        plot_theme...,
    )
    trajectory_summary = plot(
        altitude_plot,
        velocity_plot,
        flight_path_angle_plot;
        layout=(3, 1),
        size=(1300, 1100),
        margin=8Plots.mm,
    )
    savefig(trajectory_summary, joinpath(plot_dir, "drag_passage_trajectory_summary.png"))

    theta = range(0, 2pi; length=72)
    phi = range(0, pi; length=36)
    mars_x = [planet.Rp_e * cos(th) * sin(ph) / 1_000.0 for ph in phi, th in theta]
    mars_y = [planet.Rp_e * sin(th) * sin(ph) / 1_000.0 for ph in phi, th in theta]
    mars_z = [planet.Rp_e * cos(ph) / 1_000.0 for ph in phi, th in theta]

    axis_points = (
        df.sc1_pos_1 ./ 1_000.0,
        df.sc1_pos_2 ./ 1_000.0,
        df.sc1_pos_3 ./ 1_000.0,
        vec(mars_x),
        vec(mars_y),
        vec(mars_z),
    )
    axis_min = minimum(minimum, axis_points)
    axis_max = maximum(maximum, axis_points)
    axis_mid = 0.5 * (axis_min + axis_max)
    axis_half_width = 0.5 * (axis_max - axis_min)
    axis_limits = (axis_mid - axis_half_width, axis_mid + axis_half_width)
    axis_range = [axis_limits[1], axis_limits[2]]

    mars_trace = PlotlyJS.surface(
        x=mars_x,
        y=mars_y,
        z=mars_z,
        showscale=false,
        colorscale=[[0.0, "orangered"], [1.0, "orangered"]],
        opacity=1.0,
        name="Mars",
        hoverinfo="skip",
    )
    orbit_trace = PlotlyJS.scatter3d(
        x=df.sc1_pos_1 ./ 1_000.0,
        y=df.sc1_pos_2 ./ 1_000.0,
        z=df.sc1_pos_3 ./ 1_000.0,
        mode="lines",
        name="Orbit",
        line=PlotlyJS.attr(color="blue", width=4),
    )
    drag_trace = PlotlyJS.scatter3d(
        x=df.sc1_pos_1[drag_indices] ./ 1_000.0,
        y=df.sc1_pos_2[drag_indices] ./ 1_000.0,
        z=df.sc1_pos_3[drag_indices] ./ 1_000.0,
        mode="markers",
        name="Drag passage",
        marker=PlotlyJS.attr(color="gold", size=3),
    )
    orbit_layout = PlotlyJS.Layout(
        title="Odyssey Apoapsis-to-Apoapsis Orbit",
        width=1100,
        height=950,
        margin=PlotlyJS.attr(l=0, r=0, b=0, t=60),
        paper_bgcolor="white",
        plot_bgcolor="white",
        legend=PlotlyJS.attr(x=0.82, y=0.95),
        scene=PlotlyJS.attr(
            aspectmode="cube",
            bgcolor="white",
            xaxis=PlotlyJS.attr(title="X [km]", range=axis_range, backgroundcolor="white", gridcolor="lightgrey", zerolinecolor="grey"),
            yaxis=PlotlyJS.attr(title="Y [km]", range=axis_range, backgroundcolor="white", gridcolor="lightgrey", zerolinecolor="grey"),
            zaxis=PlotlyJS.attr(title="Z [km]", range=axis_range, backgroundcolor="white", gridcolor="lightgrey", zerolinecolor="grey"),
        ),
    )
    orbit_plot = PlotlyJS.plot([mars_trace, orbit_trace, drag_trace], orbit_layout)
    PlotlyJS.savefig(orbit_plot, joinpath(plot_dir, "spacecraft_orbit_3d_mars_apo_to_apo.png"))

    return (
        control_summary=joinpath(plot_dir, "drag_passage_summary.png"),
        trajectory_summary=joinpath(plot_dir, "drag_passage_trajectory_summary.png"),
        orbit=joinpath(plot_dir, "spacecraft_orbit_3d_mars_apo_to_apo.png"),
    )
end

t = @elapsed run_simulation(args; save_fields=save_fields)
csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
if args.simulation_settings.results && isfile(csv_path)
    df = CSV.read(csv_path, DataFrame)
    println("Saved $(nrow(df)) samples to $(abspath(csv_path))")
    final_orbit = _orbit_metrics_from_row(df[end, :], planet)
    println("Selected guidance mode = $(_edg_mode_name(df.sc1_edg_mode_code[end]))")
    println("Targeting active = $(df.sc1_edg_targeting_active[end] > 0.5)")
    println("Targeting switch time = $(df.sc1_edg_targeting_switch_s[end]) s")
    println("Target energy = $(df.sc1_edg_target_energy_jkg[end]) J/kg")
    println("Targeting energy bracket = ($(df.sc1_edg_bracket_min_energy_jkg[end]), $(df.sc1_edg_bracket_max_energy_jkg[end])) J/kg")
    println("Final specific orbital energy = $(final_orbit.energy_jkg) J/kg")
    println("Final apoapsis radius = $(final_orbit.apoapsis_radius_m) m")
    println("Target apoapsis radius = $(energy_config.target_apoapsis_radius_m) m")
    plot_paths = _save_energy_depletion_plots(
        csv_path,
        joinpath(args.simulation_settings.results_directory, "energy_depletion_plots"),
        planet,
    )
    println("Saved plots:")
    for path in values(plot_paths)
        println("  $(abspath(path))")
    end
end
println("COMPUTATIONAL TIME = $(t) s")
