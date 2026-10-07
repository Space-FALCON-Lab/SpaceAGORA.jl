#!/usr/bin/env julia

# Cloth solar panel deployment inside `run_simulation`: a CompliantAttachment on a bus.
#
# The four-panel tile mesh of Solar_Panel_Cloth_Deployment_Demo.jl (48 compliant bodies, folded
# in a zigzag, deployed by hinge-line motors on a smoothstep rest schedule) is mounted on the
# +x face of a rigid bus with `CompliantAttachment`. The mesh and the bus exchange forces and
# torques in both directions: the mesh bodies are integrated together with the spacecraft, and
# the joint reactions at the mount act on the bus.
#
# Two runs, no GRAM and no SPICE:
#
# 1. Comparison. A bus of 1e8 kg in free space (no gravity): effectively the fixed base of the
#    standalone demo. The saved `attachment_pose` is compared with the standalone RK4 stepping
#    (`step_compliant_multibody_rk4`) of the same model, rest schedule and actuators.
# 2. Orbit. A 100 kg, 0.6 m bus in a circular orbit (point-mass gravity, inertial attitude hold). The
#    deployment reaction torques now rock the bus, and gravity acts on every tile at its own
#    position.
#
# Outputs (output/cloth_panel_attachment_demo/): a 3D deployment HTML built from the saved
# `attachment_pose` columns, a comparison/reaction HTML, and a CSV of the plotted curves. Set
# SPACEAGORA_EXAMPLE_SMOKE=1 for a 3 s run.
#
# Solver advice: cloth meshes are stiff. This mesh is mild (fastest mode ~30 rad/s, damped), so
# the explicit :tsit5 works; with joint stiffness much larger than the orbit and attitude rates
# use :auto_stiff or :rodas5p. The attachment state is relative to the bus, so stiff meshes carry no
# roundoff from orbital-radius positions.

using DataFrames
using LinearAlgebra
using Printf
using StaticArrays
using PlotlyJS

# The mesh builder, rest schedule and plot helpers of the standalone demo (it only runs when executed directly).
include(joinpath(@__DIR__, "Solar_Panel_Cloth_Deployment_Demo.jl"))

const SM = SpaceAGORA.SimulationModel
const ATT_OUT_DIR = joinpath(REPO_ROOT, "output", "cloth_panel_attachment_demo")
const ATT_SMOKE = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
const ATT_DURATION_S = ATT_SMOKE ? 3.0 : 12.0
const ATT_SAVE_DT_S = 0.05
const ATT_STANDALONE_DT_S = 0.005
const ATT_BUS_DIMS_M = 0.6
const ATT_MOUNT_POINT = SVector{3, Float64}(0.5 * ATT_BUS_DIMS_M, 0.0, 0.0)   # +x face of the bus

function _attachment_spacecraft(; bus_mass_kg, ic, build, actuators, rest_schedule)
    bus = SM.Link(root=true, m=bus_mass_kg, dims=MVector{3, Float64}(ATT_BUS_DIMS_M, ATT_BUS_DIMS_M, ATT_BUS_DIMS_M))
    attachment = SpaceAGORA.CompliantAttachment(;
        model=build,                      # the build's state is the initial (folded) mesh state
        link=bus,
        mount_point=ATT_MOUNT_POINT,      # mount frame at the bus +x face, axes aligned with the bus
        joint_actuators=actuators,
        rest_schedule=rest_schedule,      # (out, t) -> writes every joint's rest quaternion
    )
    return SM.SpacecraftModel(;
        links=[bus], root=bus, initial_condition=ic, inertia_tensor=bus.inertia, attachments=[attachment],
    )
end

function _attachment_run(spacecraft, effectors; duration_s)
    planet = SM.make_no_gram_planet(:earth)
    args = SpaceAGORA.TelemetryVerification.make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=duration_s,
        initial_time=SM.InitialTime(year=2024, month=3, day=1, hour=12, minute=0, second=0.0),
        dynamic_effectors=effectors,
        density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=true,            # required: the mount frame follows the bus attitude
        keplerian=true,
        verbose=false,
        results=false,
        solver_config=SM.SolverConfig(solver_mode=:tsit5),
    )
    args = SM.SimConfig._with_configuration(args;
        mission_configuration=SM.MissionConfiguration(
            mission_type=args.mission_configuration.mission_type, keplerian=true, number_of_orbits=1,
            mission_time=duration_s, orientation_sim=true, num_steps_to_save=1000, data_rate=ATT_SAVE_DT_S),
        # Attachment positions are bus-relative, so no force roundoff floor; the bus itself is still at orbital radius.
        integration_tolerances=SM.IntegrationTolerances(
            reltol_orbit=1e-10, abstol_orbit=1e-9, reltol_atmosphere=1e-10, abstol_atmosphere=1e-9,
            reltol_quaternion=1e-10, abstol_quaternion=1e-9, reltol_angular_rate=1e-10, abstol_angular_rate=1e-9,
            reltol_mass=1e-10, abstol_mass=1e-9, dt_max_orbit=0.05, dt_max_atmosphere=0.05),
    )
    return SpaceAGORA.run_simulation(args; return_results=true).table
end

# Mount-frame mesh states (13 numbers per tile, velocities left zero) from the saved inertial
# `attachment_pose` columns and the bus attitude, in the layout the standalone plot helpers read.
function _mount_frame_states(df, nbodies)
    states = Vector{Vector{Float64}}(undef, nrow(df))
    for k in 1:nrow(df)
        qb = SVector{4, Float64}(df[k, "sc1_q_1"], df[k, "sc1_q_2"], df[k, "sc1_q_3"], df[k, "sc1_q_4"])
        Rb = _rot(qb)
        pos = SVector{3, Float64}(df[k, "sc1_pos_1"], df[k, "sc1_pos_2"], df[k, "sc1_pos_3"])
        mount = pos + Rb * ATT_MOUNT_POINT
        x = zeros(13 * nbodies)
        for i in 1:nbodies
            r = SVector{3, Float64}(df[k, "sc1_attachment_pose_$(7i - 6)"], df[k, "sc1_attachment_pose_$(7i - 5)"], df[k, "sc1_attachment_pose_$(7i - 4)"])
            q = SVector{4, Float64}(df[k, "sc1_attachment_pose_$(7i - 3)"], df[k, "sc1_attachment_pose_$(7i - 2)"], df[k, "sc1_attachment_pose_$(7i - 1)"], df[k, "sc1_attachment_pose_$(7i)"])
            b = 13 * (i - 1)
            x[(b + 1):(b + 3)] .= Rb' * (r - mount)
            x[(b + 4):(b + 7)] .= _quat_mul(_quat_conj(qb), q)
        end
        states[k] = x
    end
    return states
end

# The standalone demo's stepping: RK4 on a fixed base with one set of rest quaternions per step,
# evaluated `rest_offset` of a step after the step start (0 = as the standalone demo does, 0.5 = mid-step).
function _standalone_states(build, groups, init_angles, actuators, times; rest_offset=0.0)
    states = Vector{Vector{Float64}}(undef, length(times))
    x = copy(build.initial_state)
    states[1] = copy(x)
    t = 0.0
    for k in 2:length(times)
        while t < times[k] - 1e-12
            h = min(ATT_STANDALONE_DT_S, times[k] - t)
            rests = _edge_rest_quats(build.model, groups, init_angles, t + rest_offset * h)
            x = SpaceAGORA.step_compliant_multibody_rk4(build.model, x, t, h; joint_rest_quaternions=rests, joint_actuators=actuators)
            t += h
        end
        states[k] = copy(x)
    end
    return states
end

function _max_tile_difference(a, b, nbodies)
    return maximum(norm(SpaceAGORA.compliant_state_parts(a, i).r - SpaceAGORA.compliant_state_parts(b, i).r) for i in 1:nbodies)
end

function run_cloth_panel_attachment_demo()
    mkpath(ATT_OUT_DIR)
    build, init_angles, groups, _, actuators = _build_folded_panel_model()
    nbodies = length(build.model.bodies)
    rest_schedule = (out, t) -> (copyto!(out, _edge_rest_quats(build.model, groups, init_angles, t)); nothing)
    planet = SM.make_no_gram_planet(:earth)
    q0 = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)

    # 1. Comparison: huge bus, free space, against the standalone stepping.
    still = SM.CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(0.0, 0.0, 0.0); q=q0, ang_vel=SVector(0.0, 0.0, 0.0))
    t_cmp = @elapsed cmp = _attachment_run(
        _attachment_spacecraft(bus_mass_kg=1.0e8, ic=still, build=build, actuators=actuators, rest_schedule=rest_schedule), (); duration_s=ATT_DURATION_S)
    engine_states = _mount_frame_states(cmp, nbodies)
    t_ref = @elapsed ref_states = _standalone_states(build, groups, init_angles, actuators, cmp.time)
    ref_mid_states = _standalone_states(build, groups, init_angles, actuators, cmp.time; rest_offset=0.5)
    tile_diff = [_max_tile_difference(engine_states[k], ref_states[k], nbodies) for k in eachindex(cmp.time)]
    tile_diff_mid = [_max_tile_difference(engine_states[k], ref_mid_states[k], nbodies) for k in eachindex(cmp.time)]
    angles_engine = hcat((rad2deg.(_relative_hinge_angles(x, groups)) for x in engine_states)...)
    angles_ref = hcat((rad2deg.(_relative_hinge_angles(x, groups)) for x in ref_states)...)

    # 2. Orbit: free 100 kg bus in a circular orbit, inertial attitude hold.
    r0 = 7.0e6
    vc = sqrt(planet.μ / r0)
    orbit_ic = SM.CartesianInitialCondition(SVector(r0, 0.0, 0.0), SVector(0.0, vc, 0.0); q=q0, ang_vel=SVector(0.0, 0.0, 0.0))
    t_orb = @elapsed orb = _attachment_run(
        _attachment_spacecraft(bus_mass_kg=100.0, ic=orbit_ic, build=build, actuators=actuators, rest_schedule=rest_schedule),
        (SM.InverseSquaredGravityModel(),); duration_s=ATT_DURATION_S)
    orbit_states = _mount_frame_states(orb, nbodies)
    bus_q = [SVector{4, Float64}(orb[k, "sc1_q_1"], orb[k, "sc1_q_2"], orb[k, "sc1_q_3"], orb[k, "sc1_q_4"]) for k in 1:nrow(orb)]
    bus_angle_deg = [rad2deg(2 * acos(clamp(abs(q[4]), 0.0, 1.0))) for q in bus_q]
    # angular rate from the rotation between consecutive saved attitudes
    bus_rate = [2 * acos(clamp(abs(_quat_mul(_quat_conj(bus_q[k]), bus_q[k + 1])[4]), 0.0, 1.0)) / (orb.time[k + 1] - orb.time[k]) for k in 1:(nrow(orb) - 1)]
    push!(bus_rate, bus_rate[end])
    tip = hcat((_panel_tip(x) for x in orbit_states)...)
    sim = (times=collect(orb.time), states=orbit_states, tip_m=tip)

    # Outputs.
    deployment_html = _save_html(joinpath(ATT_OUT_DIR, "cloth_panel_attachment_3d.html"), _deployment_3d_plot(sim))
    traces = PlotlyJS.GenericTrace[]
    for h in 1:(PANEL_COUNT - 1)
        push!(traces, PlotlyJS.scatter(x=cmp.time, y=angles_engine[h, :], mode="lines", name="hinge $h, in the engine", xaxis="x", yaxis="y"))
        push!(traces, PlotlyJS.scatter(x=cmp.time, y=angles_ref[h, :], mode="lines", line=PlotlyJS.attr(dash="dash"), name="hinge $h, standalone", xaxis="x", yaxis="y"))
    end
    push!(traces, PlotlyJS.scatter(x=cmp.time, y=tile_diff, mode="lines", name="max tile difference, standalone as in the old demo", xaxis="x2", yaxis="y2"))
    push!(traces, PlotlyJS.scatter(x=cmp.time, y=tile_diff_mid, mode="lines", name="max tile difference, standalone with mid-step rest", xaxis="x2", yaxis="y2"))
    push!(traces, PlotlyJS.scatter(x=orb.time, y=bus_angle_deg, mode="lines", name="bus attitude change", xaxis="x3", yaxis="y3"))
    push!(traces, PlotlyJS.scatter(x=orb.time, y=rad2deg.(bus_rate), mode="lines", name="bus angular rate", xaxis="x4", yaxis="y4"))
    comparison = PlotlyJS.Plot(traces, PlotlyJS.Layout(
        title="Cloth panel deployment as a CompliantAttachment",
        grid=PlotlyJS.attr(rows=2, columns=2, pattern="independent"),
        xaxis=PlotlyJS.attr(title="time (s)"), yaxis=PlotlyJS.attr(title="relative hinge angle (deg)"),
        xaxis2=PlotlyJS.attr(title="time (s)"), yaxis2=PlotlyJS.attr(title="engine - standalone (m)", type="log"),
        xaxis3=PlotlyJS.attr(title="time (s)"), yaxis3=PlotlyJS.attr(title="bus rotation angle (deg)"),
        xaxis4=PlotlyJS.attr(title="time (s)"), yaxis4=PlotlyJS.attr(title="bus angular rate (deg/s)"),
        legend=PlotlyJS.attr(orientation="h", y=-0.2),
    ))
    comparison_html = _save_html(joinpath(ATT_OUT_DIR, "cloth_panel_attachment_comparison.html"), comparison)
    csv = DataFrame(time_s=cmp.time, tile_difference_m=tile_diff, tile_difference_midstep_rest_m=tile_diff_mid)
    for h in 1:(PANEL_COUNT - 1)
        csv[!, "hinge_$(h)_engine_deg"] = angles_engine[h, :]
        csv[!, "hinge_$(h)_standalone_deg"] = angles_ref[h, :]
    end
    csv.bus_angle_deg = bus_angle_deg
    csv.bus_rate_deg_s = rad2deg.(bus_rate)
    csv_path = joinpath(ATT_OUT_DIR, "cloth_panel_attachment_comparison.csv")
    open(csv_path, "w") do io
        println(io, join(names(csv), ","))
        for row in eachrow(csv)
            println(io, join((@sprintf("%.10g", v) for v in row), ","))
        end
    end

    println("Cloth panel attachment demo complete.")
    println("  tile bodies in the attachment    = ", nbodies)
    println("  engine run (huge bus)            = ", @sprintf("%.2f s", t_cmp), ", standalone RK4 = ", @sprintf("%.2f s", t_ref))
    println("  max tile difference vs standalone = ", @sprintf("%.3e m", maximum(tile_diff)), " over ", ATT_DURATION_S, " s (rest quaternions at the step start, as the old demo)")
    println("  ... with the rest evaluated mid-step = ", @sprintf("%.3e m", maximum(tile_diff_mid)), " (the engine evaluates the schedule at each stage time)")
    println("  final hinge angles, engine deg    = ", join((@sprintf("%.2f", a) for a in angles_engine[:, end]), ", "))
    println("  final hinge angles, standalone deg = ", join((@sprintf("%.2f", a) for a in angles_ref[:, end]), ", "))
    println("  orbit run (100 kg bus)           = ", @sprintf("%.2f s", t_orb), ", bus rotates ", @sprintf("%.3f deg", maximum(bus_angle_deg)), ", peak rate ", @sprintf("%.4f deg/s", rad2deg(maximum(bus_rate))))
    println("  3D deployment HTML               = ", deployment_html)
    println("  comparison HTML                  = ", comparison_html)
    println("  curves CSV                       = ", csv_path)
    return (comparison=cmp, orbit=orb, tile_difference_m=tile_diff, deployment_html=deployment_html, comparison_html=comparison_html, csv=csv_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_cloth_panel_attachment_demo()
end
