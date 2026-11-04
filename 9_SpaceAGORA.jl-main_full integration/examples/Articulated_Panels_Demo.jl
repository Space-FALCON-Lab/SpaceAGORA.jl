include(joinpath(@__DIR__, "common.jl"))

using DataFrames
using StaticArrays
using Plots

# Articulated solar panels: a bus with two hinged panels, no GRAM and no SPICE.
#
# The panels attach to the bus through `:hinge` joints with a spring and a damper. They start
# folded away from their rest angle (+0.20 rad and -0.10 rad about the bus x axis) while the
# bus is at rest, so the panels flap and settle and their reaction rocks the bus. The loads are
# J2 gravity only. The panels' motion is integrated in joint coordinates (`joint_q`); the bus
# keeps its own position and attitude. The script saves the hinge angles and the bus roll
# angle to a PNG.
#
# Declaring joints: build the links with their configured geometry (a panel's `r` is its COM
# in the bus frame), then connect them with `Joint(parent, p_parent, child, p_child; joint_type=...)`,
# where `p_parent`/`p_child` are the same joint point in the parent and child frames.

planet = make_no_gram_planet(:earth)

bus_dims = (2.05, 2.05, 2.8)
panel_dims = (0.01, 5.7 / 2.0, 1.0)
panel_offset_y = 2.05 / 2.0 + 5.7 / 4.0
spacecraft = make_three_body_spacecraft(
    bus_dims=bus_dims,
    panel_dims=panel_dims,
    bus_mass=620.0,
    panel_mass_each=10.0,
    panel_offset_y=panel_offset_y,
    ic=InitialCondition(ra=planet.Rp_e + 700e3, rp=planet.Rp_e + 700e3, i=51.6, ω=0.0, Ω=0.0, ν=0.0),
    prop_mass=0.0,
    id=1
)
bus, left_panel, right_panel = spacecraft.links
hinge_y = bus_dims[2] / 2.0                      # joint point on the bus edge (bus frame)
panel_half_span = panel_dims[2]                  # joint point to panel COM along y (panel frame)
for (panel, side, angle0) in ((left_panel, -1.0, 0.20), (right_panel, 1.0, -0.10))
    add_joint!(spacecraft, Joint(
        bus, SVector{3, Float64}(0.0, side * hinge_y, 0.0),
        panel, SVector{3, Float64}(0.0, -side * (panel_offset_y - hinge_y), 0.0);
        joint_type=:hinge,
        axis=[1.0, 0.0, 0.0],        # hinge axis in the bus frame
        stiffness=60.0,              # N m/rad
        damping=4.0,                 # N m s/rad
        initial_q=angle0,            # rad
    ))
end

mission_time = 60.0
args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=mission_time,
    initial_time=InitialTime(year=2024, month=3, day=1, hour=12, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredJ2GravityModel(),),
    density_model=NoAtmosphereModel(),
    ephemerides_model=SimpleEphemeridesModel(),
    orientation_sim=true,                # required for articulated spacecraft
    keplerian=true,
    EI_km=120.0,
    verbose=false,
    results=false
)
args = SimulationModel.SimConfig._with_configuration(args;
    mission_configuration=SimulationModel.MissionConfiguration(
        mission_type=args.mission_configuration.mission_type, keplerian=true, number_of_orbits=1,
        mission_time=mission_time, orientation_sim=true, num_steps_to_save=1000, data_rate=0.1),
    integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9,
        reltol_quaternion=1e-10, abstol_quaternion=1e-12, reltol_angular_rate=1e-9, abstol_angular_rate=1e-11,
        dt_max_orbit=0.2, dt_max_atmosphere=0.2))

args_eff = SpaceAGORA.TelemetryVerification._example_smoke_args(args)
t = @elapsed result = run_simulation(args_eff; return_results=true)
df = result.table
println("Saved $(nrow(df)) samples in $(round(t; digits=2)) s; columns include ",
    join(filter(n -> startswith(n, "sc1_joint") || startswith(n, "sc1_system_com"), names(df)), ", "))

# Joint coordinates follow the order of `spacecraft.joints` (left panel, right panel).
roll_deg = rad2deg.(2 .* atan.(df.sc1_q_1, df.sc1_q_4))
p1 = plot(df.time, rad2deg.(df.sc1_joint_q_1); label="left panel hinge", lw=2,
    xlabel="Time (s)", ylabel="Hinge angle (deg)", title="Panel hinge angles")
plot!(p1, df.time, rad2deg.(df.sc1_joint_q_2); label="right panel hinge", lw=2)
p2 = plot(df.time, roll_deg .- first(roll_deg); label="bus roll", lw=2, color=:black,
    xlabel="Time (s)", ylabel="Roll change (deg)", title="Bus attitude reaction")
plots_dir = joinpath(args_eff.simulation_settings.results_directory, "plots")
mkpath(plots_dir)
plot_path = joinpath(plots_dir, "articulated_panels.png")
savefig(plot(p1, p2; layout=(2, 1), size=(900, 650)), plot_path)
println("Saved plot to $(abspath(plot_path))")
