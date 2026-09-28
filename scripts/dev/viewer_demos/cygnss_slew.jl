# CYGNSS FM01: measured-wheel diagnostic and recorded-command reconstruction.
# Run with scripts/dev/run.jl; see docs/spaceagora_cygnss_reconstruction_record.md.
# Raw telemetry, fitted/memo constants, generated pages and scores remain private.
# The wheels-only adequacy check is deliberately independent of execution success.
# Demands are not automatically executed torques. The rwCmd and tqDmdCtrl comparisons
# retain the same initial state and constants; neither includes a torque-rod model.
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "cygnss_slew_telemetry.jl"))
using .CygnssSlewTelemetry
include(joinpath(@__DIR__, "cygnss_command_replay.jl"))
using .CygnssCommandReplay
using Dates
using Printf
using Statistics
using TOML
using OrdinaryDiffEq

import SpaceAGORA.TelemetryVerification as TV

const OUTDIR = demo_outdir("cygnss_slew")
const TELEMETRY_DIR = get(ENV, "SPACEAGORA_CYGNSS_DATA", joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS"))
const ADCS_PATH = get(ENV, "SPACEAGORA_DEMO_CYGNSS_SLEW_ADCS", joinpath(TELEMETRY_DIR, "cyg01_slew_adcs.feather"))
const COMMAND_PATH = get(ENV, "SPACEAGORA_DEMO_CYGNSS_COMMANDS", joinpath(TELEMETRY_DIR, "cyg01_slew_commands.feather"))
const COMMAND_REPLAY = get(ENV, "SPACEAGORA_DEMO_CYGNSS_COMMAND_REPLAY", "1") != "0"
const PV_PATH = get(ENV, "SPACEAGORA_DEMO_CYGNSS_SLEW_PV", joinpath(TELEMETRY_DIR, "cyg01_slew_pv_eci.feather"))
const CONSTANTS_PATH = get(ENV, "SPACEAGORA_DEMO_CYGNSS_SLEW_CONSTANTS", joinpath(TELEMETRY_DIR, "cyg01_adcs_constants.toml"))

# Fixed comparison window; it is not selected again to improve a replay score.
const T_CAL_S = parse(Float64, get(ENV, "SPACEAGORA_DEMO_CYGNSS_SLEW_T_CAL", "890.0"))
const T_END_S = parse(Float64, get(ENV, "SPACEAGORA_DEMO_CYGNSS_SLEW_T_END", "1100.0"))
const COMMAND_STEP_S = SLEW_COMMAND_STEP_S
const FULL_HOUR = get(ENV, "SPACEAGORA_DEMO_CYGNSS_SLEW_FULL_HOUR", "1") != "0"

# --- the drawn vehicle ------------------------------------------------------
# NASA's own public 3D model of the observatory, as the constellation page uses
# it (see `cygnss_constellation.jl` for the scale and rotation derivation and
# `data/models/README.md` for the provenance). The nadir-face assumption is
# inherited from there: a rotation of +90 degrees about model Z puts the antenna
# deck on body +z and the solar cells on body -z.
const CYGNSS_MODEL = joinpath(MODELS_DIR, "cygnss_nasa_3d_resources.glb")
const CYGNSS_MODEL_SCALE = 1.0500
const CYGNSS_MODEL_ROTATION_DEG = (0.0, 0.0, 90.0)

# The one published dimension on this page: NASA's CYGNSS mission page gives the
# observatory a 1.67 m deployed span. The link box is a cube of that span, which
# bounds the drawn model and carries no mission geometry. It sets the drag and
# SRP reference area; over the 360 s window that area moves the spacecraft by
# far less than a metre, which the script prints from the run's own drag column.
const CYGNSS_SPAN_M = 1.67
const CYGNSS_REF_AREA_M2 = CYGNSS_SPAN_M^2
# 29.00 kg is the vehicle mass stated in the SpaceAGORA CYGNSS reconstruction
# record, section 2; it is the figure `cygnss_constellation.jl` already carries.
const CYGNSS_MASS_KG = 29.0
const CYGNSS_SRP_CR = 1.3          # ASSUMPTION: a generic small-satellite reflectivity; the record states only "fixed coefficients"
# The reconstruction record (section 3) calibrated an effective drag scale of
# about 0.3 against its own geometry. It is carried here as a named field
# because the record's configuration is the one being reproduced, not because
# this window can measure it: see the displacement the script prints.
const CYGNSS_DRAG_SCALE = parse(Float64, get(ENV, "SPACEAGORA_DEMO_CYGNSS_SLEW_DRAG_SCALE", "0.3"))

const RPM_TO_RAD_S = SLEW_RPM_TO_RAD_S

# --- interpolants over the telemetry ---------------------------------------
# The same natural cubic spline the wheel effector uses, so the reference
# solution and the simulation read the telemetry through one interpolant.
# The sample times are NOT uniform (the export's 4 Hz cadence jitters by up to
# 0.08 s, and it opens with a 1 s gap), so the spline is fitted on the recorded
# times; treating them as uniform, which an earlier standalone reconstruction
# did, misplaces the middle of the hour by up to 0.8 s.
spline_of(t, y) = SM.wheel_speed_spline(t, y)

"Telemetered attitude (SpaceAGORA convention) at a telemetry time."
function telemetry_quaternion(splines, t::Float64)::SVector{4, Float64}
    q = SVector{4, Float64}(SM.wheel_spline_value(splines.q[1], t), SM.wheel_spline_value(splines.q[2], t),
        SM.wheel_spline_value(splines.q[3], t), SM.wheel_spline_value(splines.q[4], t))
    return q / norm(q)
end

"Telemetered body rate at a telemetry time."
telemetry_omega(splines, t::Float64)::SVector{3, Float64} = SVector{3, Float64}(
    SM.wheel_spline_value(splines.omega[1], t), SM.wheel_spline_value(splines.omega[2], t),
    SM.wheel_spline_value(splines.omega[3], t))

"Telemetered inertial position at a telemetry time."
telemetry_position(splines, t::Float64)::SVector{3, Float64} = SVector{3, Float64}(
    SM.wheel_spline_value(splines.pos[1], t), SM.wheel_spline_value(splines.pos[2], t),
    SM.wheel_spline_value(splines.pos[3], t))

"Telemetered inertial velocity at a telemetry time."
telemetry_velocity(splines, t::Float64)::SVector{3, Float64} = SVector{3, Float64}(
    SM.wheel_spline_value(splines.vel[1], t), SM.wheel_spline_value(splines.vel[2], t),
    SM.wheel_spline_value(splines.vel[3], t))

# --- the reference solution ------------------------------------------------

"""
    reference_attitude(tel, splines, constants, t_cal, t_end) -> function

The algebraic conservation solution this scenario must agree with. For a
torque-free body exchanging momentum only with its wheels the total angular
momentum is constant IN INERTIAL SPACE, so the body rate at every instant
follows from the attitude and the measured wheel speed with no dynamics at all:

    omega(t) = inv(I) ( C_bi(q(t)) H_inertial - H_w(t) ),

and the attitude is that rate integrated. No torque and no derivative of the
wheel signal appears: it is a different formulation of the same physics as the
simulated run, which is exactly why it is the check. The two must agree far more
closely with each other than either agrees with telemetry; if they do not, the
simulated scenario has a bug.

The inertial rotation is not a detail. Holding the conserved momentum constant
in the BODY frame instead, which the standalone reconstruction in
`extra_examples/` does, is only right for a body that is not turning; this one
turns once per revolution, so over the 360 s window the body frame carries the
conserved vector through 23 degrees and the shortcut costs about 0.8 deg of
attitude.

Returns a callable giving the SpaceAGORA-convention quaternion at a telemetry
time in `[t_cal, t_end]`.
"""
function reference_attitude(splines, constants, t_cal::Float64, t_end::Float64)
    inertia = constants.inertia
    axes = constants.wheel_axes
    Jw = constants.wheel_inertia
    wheel_momentum(t) = SVector{3, Float64}(axes * (Jw .* SVector{3, Float64}(
        SM.wheel_spline_value(splines.speeds[1], t), SM.wheel_spline_value(splines.speeds[2], t),
        SM.wheel_spline_value(splines.speeds[3], t))))
    q0 = telemetry_quaternion(splines, t_cal)
    omega0 = telemetry_omega(splines, t_cal)
    h_inertial = SM.rot(q0)' * (inertia * omega0 + wheel_momentum(t_cal))
    function rhs!(du, u, _p, t)
        q = SVector{4, Float64}(u[1], u[2], u[3], u[4])
        qn = q / norm(q)
        omega = SVector{3, Float64}(inertia \ (SM.rot(qn) * h_inertial - wheel_momentum(t)))
        du .= SM.DynamicsRotational.quaternion_derivative(omega, q)
        return nothing
    end
    prob = ODEProblem(rhs!, collect(q0), (t_cal, t_end))
    sol = solve(prob, Tsit5(); reltol=1e-12, abstol=1e-14, dtmax=0.25)
    return t -> (u = sol(t); q = SVector{4, Float64}(u[1], u[2], u[3], u[4]); q / norm(q))
end

# --- the run ---------------------------------------------------------------

println("CYGNSS FM01 commanded slew")
tel = load_slew_telemetry(ADCS_PATH, PV_PATH)
constants = load_slew_constants(CONSTANTS_PATH)
splines = (
    q=[spline_of(tel.t_rel, tel.q[k, :]) for k in 1:4],
    omega=[spline_of(tel.t_rel, tel.omega[k, :]) for k in 1:3],
    pos=[spline_of(tel.pv_t_rel, tel.pos_m[k, :]) for k in 1:3],
    vel=[spline_of(tel.pv_t_rel, tel.vel_mps[k, :]) for k in 1:3],
    speeds=[spline_of(tel.t_rel, tel.speeds_rad_s[:, k]) for k in 1:3],
)
# The kernels first: the export's clock is SPICE ephemeris time (see
# `CygnssSlewTelemetry`), so every calendar label below goes through the
# leap-second kernel the planet constructor furnishes.
planet = Earth("", SPICE_PATH)
const EPOCH_ET = slew_epoch_et(tel.t_abs[1] + T_CAL_S)
epoch_dt = slew_epoch_utc(tel.t_abs[1] + T_CAL_S)
const EPOCH_UTC = Dates.format(epoch_dt, "yyyy-mm-ddTHH:MM:SS.sss")
const DATE_LABEL = Dates.format(epoch_dt, "yyyy-mm-dd")
@printf("  export: %d samples over %.1f s, first sample %s UTC (counter %.6f s ET)\n", length(tel.t_rel), tel.t_rel[end],
    Dates.format(slew_epoch_utc(tel.t_abs[1]), "yyyy-mm-ddTHH:MM:SS.s"), tel.t_abs[1])
@printf("  body-frame nadir spread (quaternion convention check): %.4f\n", tel.nadir_spread)
@printf("  window: t_rel %.1f to %.1f s (%.1f s), commanded step at %.2f s, run epoch %s UTC\n",
    T_CAL_S, T_END_S, T_END_S - T_CAL_S, COMMAND_STEP_S, EPOCH_UTC)

# The engine's clock is resolved from the ET directly, not from the formatted
# UTC label, so nothing is rounded twice.
initial_time = initial_time_of(EPOCH_ET)

r0 = telemetry_position(splines, T_CAL_S)
v0 = telemetry_velocity(splines, T_CAL_S)
q0 = telemetry_quaternion(splines, T_CAL_S)
omega0 = telemetry_omega(splines, T_CAL_S)
oe = TV.rvtoorbitalelement(r0, v0, planet)
@printf("  initial state: alt %.1f km, i %.3f deg, period %.2f min; |omega| %.3e rad/s\n",
    (oe[1] - planet.Rp_e) / 1e3, rad2deg(oe[3]), 2pi * sqrt(oe[1]^3 / planet.μ) / 60, norm(omega0))

"""
    build_run(mission_time, outdir; with_gravity_gradient) -> (args, wheels)

Translation uses degree/order 50 Earth gravity, Sun/Moon gravity, SRP and
NRLMSISE-00 with the retained effective drag assumption and DP8. The torque
model is either measured-wheel momentum or one recorded-command replay, with
an optional gravity-gradient comparison. Both omit torque rods.
"""
function build_run(mission_time::Float64, outdir::AbstractString; with_gravity_gradient::Bool=true, command_channel::Union{Nothing,Symbol}=nothing)
    root = SM.Link(root=true, m=CYGNSS_MASS_KG,
        dims=MVector{3, Float64}(CYGNSS_SPAN_M, CYGNSS_SPAN_M, CYGNSS_SPAN_M),
        ref_area=CYGNSS_REF_AREA_M2, q=MVector{4, Float64}(q0...), ω=MVector{3, Float64}(omega0...),
        reflection_coefficient=CYGNSS_SRP_CR)
    ic = SM.CartesianInitialCondition(r0, v0; q=q0, ang_vel=omega0)
    sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=0.0,
        inertia_tensor=constants.inertia, initial_condition=ic, id=1)

    wheels = SM.ReactionWheelMomentumModel(tel.t_rel, tel.speeds_rad_s,
        constants.wheel_axes, constants.wheel_inertia; spacecraft_index=1, time_offset_s=T_CAL_S)

    # The telemetry-wheel object remains a diagnostic when commands are applied.
    # Exactly one wheel-reaction effector enters the dynamics tuple.
    replay = if command_channel === nothing
        wheels
    else
        cmd = load_wheel_commands(COMMAND_PATH; channel=command_channel)
        # Include the exact start/end by interpolation on the recorded grid.
        # Initialize H once from tachometry at the run start, not at packet zero.
        A = constants.wheel_axes ./ reshape(constants.effective_inertia_scale, 1, :)
        h0 = constants.wheel_inertia .* constants.effective_inertia_scale .* SM.wheel_speeds_rad_s(wheels, 0.0)
        full = CommandedWheelReplay(cmd.times, cmd.torque, Matrix{Float64}(I,3,3), zeros(3))
        times = [T_CAL_S; filter(t -> T_CAL_S < t < T_CAL_S+mission_time, cmd.times); T_CAL_S+mission_time]
        torque = Matrix{Float64}(undef, length(times), 3)
        for (i,t) in enumerate(times)
            torque[i,:] .= command_values(full,t).rate
        end
        CommandedWheelReplay(times, torque, A, h0; time_offset_s=T_CAL_S)
    end
    effectors = Any[
        GravitationalHarmonicsModel(50, 50, joinpath(HARMONICS_DIR, "EarthGGM05C.csv"), planet),
        NBodyGravityModel(body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet),
        SolarRadiationPressureModel(CYGNSS_SRP_CR, CYGNSS_REF_AREA_M2),
        TV.ScaledAerodynamicCoefficientfM(AerodynamicCoefficientfM(), CYGNSS_DRAG_SCALE),
        replay,
    ]
    with_gravity_gradient && push!(effectors, SM.GravityGradientTorqueModel())
    effector_tuple = Tuple(effectors)

    base = make_example_config(planet=planet, spacecraft=sc, mission_time=mission_time,
        initial_time=initial_time, dynamic_effectors=effector_tuple,
        density_model=NRLMSISE00AtmosphereModel(use_space_indices=true),
        orientation_sim=true, keplerian=false, EI_km=120.0, verbose=false, results=true,
        results_directory=String(outdir))
    args = SM.SimConfig._with_configuration(base;
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=false,
            number_of_orbits=1, mission_time=mission_time, orientation_sim=true,
            num_steps_to_save=8000, data_rate=0.25),
        dynamics_model=SM.DynamicsModel([sc], effector_tuple),
        integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-11, abstol_orbit=1e-9,
            dt_max_orbit=2.0, reltol_atmosphere=1e-11, abstol_atmosphere=1e-9, dt_max_atmosphere=2.0),
        solver_config=SM.SolverConfig(solver_mode=:dp8))
    return args, wheels
end

"""
    slew_save_fields(args, wheels) -> Vector{SaveField}

The channels the page's plot panel tells the story with, per spacecraft:
the attitude error against telemetry, the LVLH pointing angle that goes from
zero to ten degrees, the three body rates, the three wheel speeds, the three
components of the wheel momentum in the body frame, and the magnitude of the
reaction torque the wheels apply. Every one of them is a flight observable or a
quantity derived from one; none of them is a mass property.
"""
function slew_save_fields(args, wheels)
    SC = SM.SimulationCallbacks
    q_of(u) = SVector{4, Float64}(u.sc[1].q)
    ω_of(u) = SVector{3, Float64}(u.sc[1].ω)
    r_of(u) = SVector{3, Float64}(u.sc[1][1], u.sc[1][2], u.sc[1][3])
    v_of(u) = SVector{3, Float64}(u.sc[1][4], u.sc[1][5], u.sc[1][6])
    fields = [
        SC.default_save_fields(args)...,
        SC.SaveField(:attitude_error_deg,
            (u, t, integ) -> [attitude_angle_deg(q_of(u), telemetry_quaternion(splines, t + T_CAL_S))];
            per_satellite=true),
        SC.SaveField(:lvlh_pointing_deg,
            (u, t, integ) -> [lvlh_pointing_angle_deg(q_of(u), r_of(u), v_of(u))];
            per_satellite=true),
        # The same angle from the flight record, so the page's plot carries the
        # commanded maneuver and the simulated one on one axis.
        SC.SaveField(:lvlh_pointing_flight_deg,
            (u, t, integ) -> [lvlh_pointing_angle_deg(telemetry_quaternion(splines, t + T_CAL_S),
                telemetry_position(splines, t + T_CAL_S), telemetry_velocity(splines, t + T_CAL_S))];
            per_satellite=true),
        SC.SaveField(:body_rate_rad_s, (u, t, integ) -> [ω_of(u)]; per_satellite=true, column_prefix="body_rate"),
        SC.SaveField(:wheel_speed_rpm,
            (u, t, integ) -> [SM.wheel_speeds_rad_s(wheels, t) ./ RPM_TO_RAD_S];
            per_satellite=true, column_prefix="wheel_speed_rpm"),
        SC.SaveField(:wheel_momentum_nms,
            (u, t, integ) -> [SM.wheel_momentum_body(wheels, t)];
            per_satellite=true, column_prefix="wheel_momentum_nms"),
        SC.SaveField(:wheel_torque_nm,
            (u, t, integ) -> [norm(SM.wheel_reaction_torque(SM.wheel_momentum_body(wheels, t),
                SM.wheel_momentum_rate_body(wheels, t), ω_of(u)))];
            per_satellite=true),
    ]
    idx = findfirst(e -> e isa CommandedWheelReplay, args.dynamics_model.dynamic_effectors)
    if idx !== nothing
        replay = args.dynamics_model.dynamic_effectors[idx]
        push!(fields, SC.SaveField(:command_wheel_momentum_nms,
            (u,t,integ)->[command_values(replay,t).momentum];per_satellite=true))
        push!(fields, SC.SaveField(:command_reaction_torque_nm,
            (u,t,integ)->[let x=command_values(replay,t)
                SM.wheel_reaction_torque(x.momentum,x.rate,ω_of(u))
            end];per_satellite=true))
    end
    return fields
end

"""
    run_case(name, mission_time; with_gravity_gradient) -> prefix

One propagation into its own results directory, reused when it is already there
(`SPACEAGORA_DEMO_FORCE=1` reruns).
"""
function run_case(name::AbstractString, mission_time::Float64; with_gravity_gradient::Bool=true, command_channel::Union{Nothing,Symbol}=nothing)
    dir = joinpath(OUTDIR, String(name))
    mkpath(dir)
    args, wheels = build_run(mission_time, dir; with_gravity_gradient=with_gravity_gradient, command_channel=command_channel)
    return run_or_reuse!(args, dir; save_fields=slew_save_fields(args, wheels), isolate_state=false)
end

"""
    score(prefix) -> NamedTuple

Attitude and body-rate error of a run against the telemetry it is reproducing,
plus the separation from the telemetered orbit, over the run's own saved times.
"""
function score(prefix::AbstractString)
    df = DataFrame(Arrow.Table(prefix * ".feather"))
    t = Float64.(df.time)
    att = Float64.(df.sc1_attitude_error_deg)
    rate = Float64[]
    sep = Float64[]
    across = Float64[]
    for i in eachindex(t)
        om = SVector{3, Float64}(df.sc1_body_rate_1[i], df.sc1_body_rate_2[i], df.sc1_body_rate_3[i])
        push!(rate, norm(om - telemetry_omega(splines, t[i] + T_CAL_S)))
        r = SVector{3, Float64}(df.sc1_pos_1[i], df.sc1_pos_2[i], df.sc1_pos_3[i])
        v = SVector{3, Float64}(df.sc1_vel_1[i], df.sc1_vel_2[i], df.sc1_vel_3[i])
        d = r - telemetry_position(splines, t[i] + T_CAL_S)
        push!(sep, norm(d))
        # The along-track component is the one the navigation fixes' quarter
        # second of time-stamp resolution cannot pin (see `CygnssSlewTelemetry`),
        # so the radial and cross-track part is reported separately.
        push!(across, norm(d - dot(d, v / norm(v)) * (v / norm(v))))
    end
    return (df=df, t=t, att=att, rate=rate, sep=sep, across=across)
end

init_nrlmsise_space_indices!()
window_s = T_END_S - T_CAL_S

# The measured external torque over the window, from the flight record alone,
# so the page states what the wheels-only replay leaves out rather than
# assuming it small. See the header and `cygnss_slew_checks.jl`.
ledger = momentum_ledger(tel, constants; window=(T_CAL_S, T_END_S))
@printf("  momentum ledger over the window: net external torque %.3e N m (gravity gradient alone %.3e N m); omitted %.3e against exchanged %.3e N m s, ratio %.2f\n",
    norm(ledger.drift_nm), norm(ledger.gravity_gradient_nm), ledger.omitted_nms, ledger.exchanged_nms, ledger.ratio)

println("\nwheels-only diagnostic (physical adequacy may FAIL): ", round(window_s; digits=1), " s")
prefix = run_case("window", window_s; with_gravity_gradient=false)
main = score(prefix)
@printf("  attitude error vs telemetry: mean %.3f  rms %.3f  max %.3f deg\n",
    mean(main.att), sqrt(mean(main.att .^ 2)), maximum(main.att))
@printf("  body-rate error vs telemetry: rms %.3e  max %.3e rad/s\n", sqrt(mean(main.rate .^ 2)), maximum(main.rate))
@printf("  position vs the telemetered orbit: start %.1f m, end %.1f m, max %.1f m\n",
    main.sep[1], main.sep[end], maximum(main.sep))
@printf("  the same, off the track (radial and cross-track only, which the fix time stamps do resolve): end %.1f m, max %.1f m\n",
    main.across[end], maximum(main.across))
@printf("  LVLH pointing angle, simulated: %.3f deg at the start, %.3f at the end (peak %.3f)\n",
    main.df.sc1_lvlh_pointing_deg[1], main.df.sc1_lvlh_pointing_deg[end], maximum(main.df.sc1_lvlh_pointing_deg))
@printf("  LVLH pointing angle, flight:    %.3f deg at the start, %.3f at the end (peak %.3f)\n",
    main.df.sc1_lvlh_pointing_flight_deg[1], main.df.sc1_lvlh_pointing_flight_deg[end],
    maximum(main.df.sc1_lvlh_pointing_flight_deg))
let d = hypot.(main.df.sc1_drag_1, main.df.sc1_drag_2, main.df.sc1_drag_3)
    @printf("  drag over the window: mean %.3e N on %.2f kg, so %.3f m of displacement at most\n",
        mean(d), CYGNSS_MASS_KG, 0.5 * (mean(d) / CYGNSS_MASS_KG) * window_s^2)
end

# Where the window could have ended, and what each choice costs. The run is one
# trajectory, so an earlier end is a truncation of it and the table needs no
# extra propagation.
println("\nwindow end against error and pointing angle (the wheels-only diagnostic, truncated):")
for te in 1000.0:25.0:T_END_S
    k = findall(x -> x <= te - T_CAL_S, main.t)
    isempty(k) && continue
    i = k[end]
    @printf("    end t_rel %6.0f s: error at the end %6.3f deg (mean %6.3f, max %6.3f) | LVLH flight %6.3f, simulated %6.3f deg\n",
        te, main.att[i], mean(main.att[k]), maximum(main.att[k]),
        main.df.sc1_lvlh_pointing_flight_deg[i], main.df.sc1_lvlh_pointing_deg[i])
end

println("\ncomparison run: the same, with gravity-gradient torque added")
prefix_gg = run_case("window_gravity_gradient", window_s; with_gravity_gradient=true)
with_gg = score(prefix_gg)
@printf("  attitude error vs telemetry: mean %.3f  rms %.3f  max %.3f deg\n",
    mean(with_gg.att), sqrt(mean(with_gg.att .^ 2)), maximum(with_gg.att))

println("\nagreement with the algebraic conservation solution (the reference)")
ref_q = reference_attitude(splines, constants, T_CAL_S, T_CAL_S + window_s)
ref_vs_sim = Float64[]
ref_vs_tel = Float64[]
gg_shift = Float64[]
for i in eachindex(main.t)
    tt = main.t[i] + T_CAL_S
    qr = ref_q(tt)
    qw = SVector{4, Float64}(main.df.sc1_q_1[i], main.df.sc1_q_2[i], main.df.sc1_q_3[i], main.df.sc1_q_4[i])
    qg = SVector{4, Float64}(with_gg.df.sc1_q_1[i], with_gg.df.sc1_q_2[i], with_gg.df.sc1_q_3[i], with_gg.df.sc1_q_4[i])
    push!(ref_vs_sim, attitude_angle_deg(qw, qr))
    push!(ref_vs_tel, attitude_angle_deg(qr, telemetry_quaternion(splines, tt)))
    push!(gg_shift, attitude_angle_deg(qg, qw))
end
@printf("  the wheels-only diagnostic vs the algebraic reference: mean %.2e  max %.2e deg\n",
    mean(ref_vs_sim), maximum(ref_vs_sim))
@printf("  the algebraic reference vs telemetry:         mean %.3f  max %.3f deg\n",
    mean(ref_vs_tel), maximum(ref_vs_tel))
@printf("  gravity gradient moves the simulated attitude by: mean %.3f  max %.3f deg\n",
    mean(gg_shift), maximum(gg_shift))

full_hour = nothing
full_hour_gg = nothing
if FULL_HOUR
    hour_s = tel.t_rel[end] - T_CAL_S
    println("\nrange of validity: the same setup carried to the end of the hour (", round(hour_s; digits=0), " s)")
    full_hour = score(run_case("full_hour", hour_s; with_gravity_gradient=false))
    full_hour_gg = score(run_case("full_hour_gravity_gradient", hour_s; with_gravity_gradient=true))
    @printf("  wheel momentum alone:   mean %.2f  rms %.2f  max %.2f deg\n",
        mean(full_hour.att), sqrt(mean(full_hour.att .^ 2)), maximum(full_hour.att))
    @printf("  with gravity gradient:  mean %.2f  rms %.2f  max %.2f deg\n",
        mean(full_hour_gg.att), sqrt(mean(full_hour_gg.att .^ 2)), maximum(full_hour_gg.att))
    for mark in (window_s, 600.0, 1200.0, 1800.0, 2400.0, hour_s)
        i = findlast(t -> t <= mark, full_hour.t)
        i === nothing && continue
        j = findlast(t -> t <= mark, full_hour_gg.t)
        @printf("    t_rel %6.0f s: %7.2f deg wheels only, %7.2f deg with gravity gradient\n",
            full_hour.t[i] + T_CAL_S, full_hour.att[i], full_hour_gg.att[j])
    end
end

# Keep completion, model adequacy and observed fit as separate outcomes.
replay_results = Dict{String,Any}()
if COMMAND_REPLAY
    isfile(COMMAND_PATH) || error("Missing command extract $COMMAND_PATH; see the reconstruction record")
    for channel in (:rwCmd, :tqDmdCtrl)
        label = "command_$(channel)"
        scored = score(run_case(label, window_s; with_gravity_gradient=true, command_channel=channel))
        replay_results[label] = Dict("mean_attitude_error_deg"=>mean(scored.att),
            "max_attitude_error_deg"=>maximum(scored.att), "channel"=>string(channel),
            "meaning"=>channel == :rwCmd ? "recorded executed wheel-command hypothesis" : "control-demand diagnostic, not executed torque")
        @printf("  %s: mean %.3f, max %.3f deg (reconstruction, not actuator validation)\n", label, mean(scored.att), maximum(scored.att))
    end
end
open(joinpath(OUTDIR,"reconstruction_summary.toml"), "w") do io
    TOML.print(io, Dict("wheels_only"=>Dict("adequacy_pass"=>isfinite(ledger.ratio) && ledger.ratio<1,
        "omitted_to_exchanged_momentum_ratio"=>ledger.ratio,
        "mean_attitude_error_deg"=>mean(main.att), "max_attitude_error_deg"=>maximum(main.att)),
        "command_replays"=>replay_results,
        "momentum_scale_qualification"=>"Initial wheel momentum uses the tachometer with memo-scale inertia and per-wheel multipliers; command momentum changes use gain-one recorded commands. These scales differ, so attitude scores are not a matched test of the executed-channel hypothesis. No inertia rebaseline is performed.",
        "qualification"=>"Private reconstruction with supplied constants; no independent actuator or flight-accuracy acceptance."))
end

# --- the page --------------------------------------------------------------
# The page is exported from a trimmed copy of the results, as the constellation
# page is: the per-spacecraft mass column is dropped and the link's mass is
# zeroed, so no mass property reaches the page in any field. The link box stays,
# because the viewer falls back to it when a model cannot be parsed, and it is
# the published deployed span rather than any mission geometry.
#
# The body-frame wheel-momentum columns go too. They are a diagnostic worth
# keeping in the run's own results, but on a page that also carries the wheel
# speeds their ratio is the wheel inertia times the spin-axis matrix, and this
# page carries no inertia of any kind. The three wheel SPEEDS stay: they are the
# flight observable, and they tell a reader the same story.
const PAGE_DIR = joinpath(OUTDIR, "page")
mkpath(PAGE_DIR)
page_prefix = joinpath(PAGE_DIR, "simulation_results")
let doc = JSON.parsefile(prefix * "_scene.json")
    doc["spacecraft"][1]["name"] = "CYGNSS FM01"
    for link in doc["spacecraft"][1]["links"]
        link["mass_kg"] = 0.0
    end
    open(page_prefix * "_scene.json", "w") do io
        JSON.print(io, doc)
    end
    full = DataFrame(Arrow.Table(prefix * ".feather"))
    dropped = [c for c in names(full) if occursin(r"^sc\d+_mass$", c) || occursin(r"^sc\d+_wheel_momentum_nms_\d$", c)]
    Arrow.write(page_prefix * ".feather", select(full, Not(dropped)))
    println("\npage copy: dropped ", length(dropped), " column(s): ", join(dropped, ", "))
end

# The ghost carries BOTH the flown orbit and the flown attitude, so the
# translucent copy beside the simulated spacecraft is what the vehicle actually
# did -- where it was and which way it was facing. At 4 Hz over six minutes the
# whole flown arc is 1441 samples, small enough to embed entire.
# The ghost's time base is the 4 Hz attitude, which is what this page is about;
# the position and velocity are the 1 Hz navigation fixes read onto it through
# the same splines the run scores against, because the position column of the
# export repeats each fix four times (see `CygnssSlewTelemetry`).
ghost_idx = [i for i in eachindex(tel.t_rel) if T_CAL_S - 1e-6 <= tel.t_rel[i] <= T_CAL_S + window_s + 1e-6]
references = [(
    name="CYGNSS FM01 flight telemetry",
    t_s=tel.t_rel[ghost_idx] .- T_CAL_S,
    pos_m=reduce(hcat, [telemetry_position(splines, tel.t_rel[i]) for i in ghost_idx]),
    vel_mps=reduce(hcat, [telemetry_velocity(splines, tel.t_rel[i]) for i in ghost_idx]),
    q=tel.q[:, ghost_idx],
    target=1, color="#ffb347", opacity=0.45,
)]
println("ghost: ", length(ghost_idx), " flown samples with position and attitude")

# Frame and trail. The commanded step is a rotation in INERTIAL space and the
# window is 0.064 of a revolution, so an inertial frame is the one in which the
# thing the page is about -- the vehicle turning about a fixed axis -- is the
# only motion that reads; a planet-fixed frame would add half a degree per
# second of Earth rotation to it for no gain. The trail is the whole window, so
# the arc behind the spacecraft shows where the six minutes started rather than
# a rolling stub. Ground tracks stay off: six minutes of sub-satellite track is
# a short arc that says nothing about an attitude maneuver.
# The channels the page carries beside the trajectory, in the order they tell
# the story: what the vehicle was commanded to do, how well the simulation did
# it, and the wheel speeds that drove it. Every one is a flight observable or an
# angle derived from one. The body-frame wheel momentum is deliberately not
# among them (see the page copy above).
channels = [
    (column="lvlh_pointing_flight_deg", label="LVLH pointing, flight", unit="deg", digits=3),
    (column="lvlh_pointing_deg", label="LVLH pointing, simulated", unit="deg", digits=3),
    (column="attitude_error_deg", label="attitude error vs flight", unit="deg", digits=3),
    (column="body_rate_1", label="body rate x", unit="rad/s", digits=6),
    (column="body_rate_2", label="body rate y", unit="rad/s", digits=6),
    (column="body_rate_3", label="body rate z", unit="rad/s", digits=6),
    (column="wheel_speed_rpm_1", label="wheel 1 speed", unit="rpm", digits=1),
    (column="wheel_speed_rpm_2", label="wheel 2 speed", unit="rpm", digits=1),
    (column="wheel_speed_rpm_3", label="wheel 3 speed", unit="rpm", digits=1),
    (column="wheel_torque_nm", label="wheel reaction torque", unit="N m", digits=7),
]

html = export_visualization(page_prefix; max_frames=2000, trail_s=window_s, texture_resolution="4k",
    frame=:inertial, ground_tracks=false, channels=channels,
    title="AGORA CYGNSS FM01 · the commanded slew, $(DATE_LABEL)",
    models=Dict(1 => CYGNSS_MODEL), model_scale=CYGNSS_MODEL_SCALE,
    model_rotation_deg=Dict(1 => CYGNSS_MODEL_ROTATION_DEG),
    references=references)
println("html: ", html, " ", filesize(html))

validity = full_hour === nothing ? "" : @sprintf(
    "The window ends where the maneuver does. Carried on open loop to the end of the hour the same run drifts to %.0f deg mean, because the external torques it does not carry accumulate and nothing corrects them. ",
    mean(full_hour.att[full_hour.t .> window_s]))
cdn = build_cdn_page(html, joinpath(OUTDIR, "artifact.html"), "AGORA CYGNSS FM01 Slew",
    "AGORA CYGNSS FM01 · the commanded slew, $(DATE_LABEL) $(Dates.format(epoch_dt, "HH:MM")) UTC",
    @sprintf("%.0f km, i %.2f°, %.1f min; EarthGGM05C 50x50, Sun + Moon, SRP, NRLMSISE-00; reaction-wheel momentum exchange, no external torque, dp8",
        (oe[1] - planet.Rp_e) / 1e3, rad2deg(oe[3]), 2pi * sqrt(oe[1]^3 / planet.μ) / 60),
    @sprintf("%.0f s, the commanded step at %+.0f s", window_s, COMMAND_STEP_S - T_CAL_S),
    "Private wheels-only diagnostic, not a validated flight reconstruction. " *
    @sprintf("The omitted/exchanged momentum ratio is %.2f; physical adequacy: %s. ", ledger.ratio, ledger.ratio < 1 ? "PASS" : "FAIL") *
    "Commanded-wheel replays and their qualifications are saved separately in reconstruction_summary.toml. " *
    "Their initial and changing wheel momentum use different scales, so the attitude scores are not a matched channel comparison. " *
    "Click the spacecraft for its channels; the LVLH pointing angle is the maneuver in one number. Drag to orbit, wheel to zoom, Space to pause, F to follow.")
cdn === nothing || println("cdn: ", cdn, " ", filesize(cdn))
