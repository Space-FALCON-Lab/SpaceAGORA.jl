##
# Shared setup for the CYGNSS attitude/torque reconstruction study: telemetry
# loading, custom SpaceAGORA effectors, spacecraft construction, and the
# `run_cygnss_case` simulation runner. `include`-d by every `run_*.jl` script
# in this folder -- see README.md for the full narrative and results.
##

include(joinpath(@__DIR__, "..", "..", "..", "examples", "common.jl"))
using SPICE
using StaticArrays
using LinearAlgebra
using DataFrames
using Arrow
using Interpolations
using Statistics

setup_gram_example!()

const STUDY_DIR = @__DIR__
const DATA_DIR = joinpath(STUDY_DIR, "data")
const PLOTS_DIR = joinpath(STUDY_DIR, "plots")
mkpath(PLOTS_DIR)

planet = Earth("", SPICE_PATH)

# Real calendar epoch of t_rel=0 in this telemetry file, derived from the raw
# GPS receiver telemetry (ADCS_BUS_FSW_INP.sensorData.gpsData.gpsWeek/gpsSec
# at the first sample: week 2386, gpsSec 521835.99999998), converted from GPS
# time to UTC (GPS is 18s ahead of UTC as of the 2016-12-31 leap second, still
# in effect through at least this epoch). An earlier, unrelated placeholder
# date (2018-12-15) was used here before this was derived; it made no
# measurable difference to propagated position error once corrected (ruled
# out as an explanation -- see "GPS position as ground truth" below for what
# actually explains the multi-km discrepancy), but the correct epoch is kept
# since it's simply more accurate.
const TELEMETRY_INITIAL_TIME = InitialTime(year=2025, month=10, day=4, hour=0, minute=56, second=58.0)
const _EPHEM_MODEL = SpiceEphemeridesModel()
const _ET_START = SimulationModel.SimulationCallbacks.ephemerides_time_seconds(TELEMETRY_INITIAL_TIME, _EPHEM_MODEL)

# ==============================================================================
# Telemetry loading
#
# Source: CYGNSS_slew_data.feather (repo root, HDF5-format flight telemetry,
# 1 hour @ ~1 Hz spanning a slew maneuver at t=901.75s), pre-extracted into
# data/. See README.md for the extraction provenance (raw HDF5 paths, dplCmd
# units decision, etc.).
# ==============================================================================

df = DataFrame(Arrow.Table(joinpath(DATA_DIR, "cygnss_slew_data_extracted.feather")))
sort!(df, :t_rel)

t_all = Float64.(df.t_rel)
t_range = range(t_all[1], t_all[end], length=length(t_all))

# Quaternion continuity unwrap (telemetry is scalar-first [qs, qx, qy, qz]).
qraw = Vector{Float64}.(eachrow(hcat(df.q_eci_0, df.q_eci_1, df.q_eci_2, df.q_eci_3)))
for i in 2:length(qraw)
    dot(qraw[i], qraw[i - 1]) < 0.0 && (qraw[i] = -qraw[i])
end
qm = reduce(hcat, qraw)'
q1_itp = cubic_spline_interpolation(t_range, Float64.(qm[:, 1])) # scalar
q2_itp = cubic_spline_interpolation(t_range, Float64.(qm[:, 2]))
q3_itp = cubic_spline_interpolation(t_range, Float64.(qm[:, 3]))
q4_itp = cubic_spline_interpolation(t_range, Float64.(qm[:, 4]))
q_truth_scalarfirst(t) = SVector{4, Float64}(q1_itp(t), q2_itp(t), q3_itp(t), q4_itp(t))
q_truth_scalarlast(t) = (qsf = q_truth_scalarfirst(t); SVector{4, Float64}(qsf[2], qsf[3], qsf[4], qsf[1]))

w1_itp = cubic_spline_interpolation(t_range, Float64.(df.w_eci_0))
w2_itp = cubic_spline_interpolation(t_range, Float64.(df.w_eci_1))
w3_itp = cubic_spline_interpolation(t_range, Float64.(df.w_eci_2))

Om1_itp = cubic_spline_interpolation(t_range, Float64.(df.Omega_rw_0))
Om2_itp = cubic_spline_interpolation(t_range, Float64.(df.Omega_rw_1))
Om3_itp = cubic_spline_interpolation(t_range, Float64.(df.Omega_rw_2))

dpl0_itp = cubic_spline_interpolation(t_range, Float64.(df.dplCmd_0))
dpl1_itp = cubic_spline_interpolation(t_range, Float64.(df.dplCmd_1))
dpl2_itp = cubic_spline_interpolation(t_range, Float64.(df.dplCmd_2))
m_body_itp(t) = SVector{3, Float64}(dpl0_itp(t), dpl1_itp(t), dpl2_itp(t)) # A*m^2, body frame, raw dplCmd

# ==============================================================================
# GPS position as ground truth.
#
# Originally used ADCS_BUS_FSW_TLM.adsAttFilter.pv_eci (the onboard attitude
# filter's own propagated/estimated position) as position ground truth.
# Direct cross-check against the raw GPS receiver measurement
# (ADCS_BUS_FSW_INP.sensorData.gpsData.PV_ecef, converted to ECI here via
# SpaceAGORA's own SPICE frame) at matching timestamps showed adsAttFilter's
# "ECI" disagrees with the raw GPS fix by 5-14 km with ZERO propagation
# involved -- a ground-truth frame/convention discrepancy, not a dynamics
# error. Sweeping the ECEF->ECI conversion epoch to find the best-fit timing
# offset only reduced the gap to ~5.5 km (not ~0), ruling out a pure clock/
# leap-second bug and pointing to a genuine frame-convention difference
# between full SPICE-grade Earth orientation (ITRF<->J2000, used here) and
# whatever simplified onboard sidereal-rotation model adsAttFilter uses
# internally. That ~5.5 km is essentially the same size as the "position
# error" previously being chased through gravity-model/IC fidelity in
# run_position_accuracy.jl, none of which could have addressed a ground-truth
# problem. Ground truth is now the raw GPS fix, converted to ECI ourselves
# (self-consistent with how the propagated trajectory is scored, since both
# use the same SPICE frame) instead of adsAttFilter.pv_eci.
#
# GPS position updates at ~1 Hz and is held between faster-cadence telemetry
# samples (adsAttFilter.time updates faster than gpsSec changes) -- dedupe to
# actual fix instants before building the interpolant.
# ==============================================================================

gps_raw = DataFrame(Arrow.Table(joinpath(DATA_DIR, "cygnss_slew_gps_raw_extracted.feather")))
sort!(gps_raw, :t_rel)

_gps_keep = trues(nrow(gps_raw))
for k in 2:nrow(gps_raw)
    if gps_raw.gps_sec[k] == gps_raw.gps_sec[k - 1] && gps_raw.gps_week[k] == gps_raw.gps_week[k - 1]
        _gps_keep[k] = false
    end
end
gps_eci = gps_raw[_gps_keep, :]

gps_eci.r_eci_0 = zeros(nrow(gps_eci))
gps_eci.r_eci_1 = zeros(nrow(gps_eci))
gps_eci.r_eci_2 = zeros(nrow(gps_eci))
gps_eci.v_eci_0 = zeros(nrow(gps_eci))
gps_eci.v_eci_1 = zeros(nrow(gps_eci))
gps_eci.v_eci_2 = zeros(nrow(gps_eci))
for k in 1:nrow(gps_eci)
    et = _ET_START + gps_eci.t_rel[k]
    l_pi = SimulationModel.planet_frame_lpi(planet, et, _EPHEM_MODEL)
    r_ecef = SVector{3, Float64}(gps_eci.r_ecef_0[k], gps_eci.r_ecef_1[k], gps_eci.r_ecef_2[k])
    v_ecef = SVector{3, Float64}(gps_eci.v_ecef_0[k], gps_eci.v_ecef_1[k], gps_eci.v_ecef_2[k])
    r_ii = l_pi' * r_ecef
    v_ii = l_pi' * v_ecef + cross(planet.ω, r_ii)
    gps_eci.r_eci_0[k], gps_eci.r_eci_1[k], gps_eci.r_eci_2[k] = r_ii
    gps_eci.v_eci_0[k], gps_eci.v_eci_1[k], gps_eci.v_eci_2[k] = v_ii
end

gps_t_range = range(gps_eci.t_rel[1], gps_eci.t_rel[end], length=nrow(gps_eci))
r1_itp = cubic_spline_interpolation(gps_t_range, gps_eci.r_eci_0)
r2_itp = cubic_spline_interpolation(gps_t_range, gps_eci.r_eci_1)
r3_itp = cubic_spline_interpolation(gps_t_range, gps_eci.r_eci_2)
v1_itp = cubic_spline_interpolation(gps_t_range, gps_eci.v_eci_0)
v2_itp = cubic_spline_interpolation(gps_t_range, gps_eci.v_eci_1)
v3_itp = cubic_spline_interpolation(gps_t_range, gps_eci.v_eci_2)
r_truth(t) = SVector{3, Float64}(r1_itp(t), r2_itp(t), r3_itp(t))
v_truth(t) = SVector{3, Float64}(v1_itp(t), v2_itp(t), v3_itp(t))

# ==============================================================================
# Validated constants (from the maneuver-window conservation-law regression)
# ==============================================================================

const J_RW = SMatrix{3, 3, Float64}(
    -0.885246, 0.0251202, -0.14032,
    -1.18407, -0.0215877, 0.0715373,
    -0.963195, 0.0369281, -0.275273
)
const I_WHEEL = 18.0e-3 / (6000.0 * 2.0 * pi / 60.0) # kg*m^2
const RPM_TO_RAD_S = 2.0 * pi / 60.0

const INERTIA_OVERRIDE = SMatrix{3, 3, Float64}(
    1.4e6, -1.71e4, 8.08e3,
    -1.71e4, 8.19e5, -5.35e3,
    8.08e3, -5.35e3, 1.95e6
) .* 1e-6

Om_rads(t) = SVector{3, Float64}(Om1_itp(t), Om2_itp(t), Om3_itp(t)) .* RPM_TO_RAD_S
H_wheels(t) = J_RW * (I_WHEEL .* Om_rads(t))
Hdot_wheels(t) = J_RW * (I_WHEEL * RPM_TO_RAD_S .* SVector{3, Float64}(
    Interpolations.gradient(Om1_itp, t)[1],
    Interpolations.gradient(Om2_itp, t)[1],
    Interpolations.gradient(Om3_itp, t)[1],
))

# ==============================================================================
# Custom ControlModel: reaction-wheel momentum-exchange replay from telemetry.
#
# The engine integrates the standard torque-free rigid body Euler equation
# I*w_dot = tau_total - w x (I*w) (src/dynamics/rotational/rigid_body_dynamics.jl)
# with no built-in notion of wheel momentum. To replay the wheel-telemetry-
# based exchange as an external torque contribution, supply the full
# momentum-transfer correction: tau_control = -w x H_wheels(t) - Hdot_wheels(t).
# ==============================================================================

struct WheelMomentumReplayControlModel
    t_offset::Float64 # engine integration time is relative to sim start; telemetry interpolants are absolute
end

function SimulationModel.calcControlForceTorque(
    model::WheelMomentumReplayControlModel,
    u,
    p::ODEParams,
    i::Int64,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    ω = SVector{3, Float64}(u.ω)
    t_abs = t + model.t_offset
    Hw = H_wheels(t_abs)
    Hwdot = Hdot_wheels(t_abs)
    τ = -cross(ω, Hw) - Hwdot
    return SVector{3, Float64}(0.0, 0.0, 0.0), τ
end

function SimulationModel.calcControlEffect!(
    model::WheelMomentumReplayControlModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64,
)
    return nothing
end

# ==============================================================================
# Custom ControlModel: kinematic torque back-out from the measured attitude
# trajectory, in a fully torque-free (no gravity-gradient/aero/magnetic)
# environment.
#
# Modified Euler's equation with reaction-wheel momentum:
#   I*w_dot + w x (I*w + H_wheels) + Hdot_wheels = tau_ext
# In a torque-free environment (tau_ext = 0), the wheels are the only torque
# source, so the *combined* wheel reaction term (w x H_wheels + Hdot_wheels)
# must equal -(I*w_dot + w x (I*w)) evaluated on the MEASURED trajectory. This
# backs the net wheel-reaction torque directly out of the measured q(t)/w(t)
# kinematics (w and its derivative from the telemetry spline), without going
# through the Omega_rw/J_RW wheel-telemetry model at all. Applying the
# negative of that as the control torque (i.e. tau_control = -(w x H_wheels +
# Hdot_wheels) = I*w_dot_meas + w_meas x (I*w_meas)) and integrating forward
# reproduces the measured trajectory near-exactly (up to spline/
# differentiation and integration error) -- this is inverse dynamics, not an
# independent physical prediction. See README.md ("Kinematic torque back-out")
# for why this result matters diagnostically.
# ==============================================================================

ω_meas(t) = -1.0 .* SVector{3, Float64}(w1_itp(t), w2_itp(t), w3_itp(t))
ω̇_meas(t) = -1.0 .* SVector{3, Float64}(
    Interpolations.gradient(w1_itp, t)[1],
    Interpolations.gradient(w2_itp, t)[1],
    Interpolations.gradient(w3_itp, t)[1],
)
τ_kinematic_backout(t) = INERTIA_OVERRIDE * ω̇_meas(t) + cross(ω_meas(t), INERTIA_OVERRIDE * ω_meas(t))

struct KinematicTorqueReplayControlModel
    t_offset::Float64
end

function SimulationModel.calcControlForceTorque(
    model::KinematicTorqueReplayControlModel,
    u,
    p::ODEParams,
    i::Int64,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    τ = τ_kinematic_backout(t + model.t_offset)
    return SVector{3, Float64}(0.0, 0.0, 0.0), τ
end

function SimulationModel.calcControlEffect!(
    model::KinematicTorqueReplayControlModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64,
)
    return nothing
end

# ==============================================================================
# Custom dynamic effector: telemetry-driven magnetorquer torque.
#
# Reuses the exact same src/ physics wired into MagneticTorqueRodModel
# (get_magnetic_field_dipole, calculate_magnetic_torque) but drives the
# commanded dipole moment from a time interpolant of real dplCmd telemetry
# instead of a static per-spacecraft Magnet configuration.
# ==============================================================================

# `rot` (scalar-last quaternion -> DCM) is an internal, non-exported helper
# (src/core/numerics/quaternion_utils.jl, textually `include`-d into each
# consuming src/ module rather than exported). Reproduced verbatim here since
# it isn't reachable through `using .SimulationModel` from outside src/.
@inline function rot(q::AbstractVector{<:Float64})::SMatrix{3, 3, Float64}
    q1, q2, q3, q4 = q[1], q[2], q[3], q[4]
    q1q1 = q1 * q1; q2q2 = q2 * q2; q3q3 = q3 * q3; q4q4 = q4 * q4
    q1q2 = q1 * q2; q1q3 = q1 * q3; q1q4 = q1 * q4
    q2q3 = q2 * q3; q2q4 = q2 * q4; q3q4 = q3 * q4
    return SMatrix{3, 3, Float64}(
        q1q1 - q2q2 - q3q3 + q4q4, 2 * (q1q2 - q3q4), 2 * (q1q3 + q2q4),
        2 * (q1q2 + q3q4), -q1q1 + q2q2 - q3q3 + q4q4, 2 * (q2q3 - q1q4),
        2 * (q1q3 - q2q4), 2 * (q2q3 + q1q4), -q1q1 - q2q2 + q3q3 + q4q4
    )
end

struct TelemetryMagnetorquerModel{F} <: SimulationModel.AbstractForceTorqueModel
    m_itp::F
    t_offset::Float64 # engine integration time is relative to sim start; telemetry interpolants are absolute
end

@inline SimulationModel.environment_requirements(::TelemetryMagnetorquerModel) =
    SimulationModel.EffectorEnvironmentRequirements(planet_frame=true)

@inline function SimulationModel.wrench(
    model::TelemetryMagnetorquerModel,
    x::SimulationModel.StateSample,
    env::SimulationModel.EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    x.q_ib === nothing && return zero3, zero3
    m_body = model.m_itp(t + model.t_offset)
    iszero(m_body) && return zero3, zero3
    planet_frame = env.planet_frame
    planet_frame === nothing && throw(ArgumentError("TelemetryMagnetorquerModel wrench requires env.planet_frame."))
    B_ii = get_magnetic_field_dipole(planet_frame.pos_pp, MMatrix{3, 3, Float64}(planet_frame.l_pi))
    B_body = rot(x.q_ib) * B_ii
    τ = calculate_magnetic_torque(m_body, B_body)
    return zero3, τ
end

# ==============================================================================
# Spacecraft construction (CYGNSS dimensions/mass from extra_examples/CYGNSS_test.jl,
# custom inertia tensor override from the validated conservation-law fit).
# ==============================================================================

const BUS_DIMS = (0.202, 0.521, 0.641)
const PANEL_DIMS = (0.497, 0.0005, 0.521)
const BUS_MASS = 28.94
const PANEL_MASS = 1.0e-6
const PANEL_OFFSET_Y = 0.569

# ==============================================================================
# SMA-averaged initial condition (position accuracy fix).
#
# Replicates the fix in test/gmat_scenario_matrix.jl's
# _build_cygnss_48hr_reference exactly: the raw single-sample IC (r0, v0)
# carries real telemetry noise -- vis-viva SMA (an exact constant of
# unperturbed 2-body motion) jitters sample-to-sample, and because this
# orbit's along-track behavior is extremely sensitive to velocity-magnitude/
# SMA error, that noise alone is enough to explain multi-km along-track drift
# over the propagation window. The fix keeps the raw sample's position AND
# velocity DIRECTION unchanged and corrects only the velocity MAGNITUDE so
# vis-viva SMA matches a short-window (N=5 raw telemetry rows) average of
# per-sample SMAs computed directly from telemetry's own reported (r, v)
# pairs -- no differentiation or curve-fitting, just noise averaging. N=5 was
# swept directly against propagated position RMSE in the other study
# (N=1..3 under-corrects, N>=12 over-corrects into real dynamics; N=4..8
# forms a consistent ~1.2-1.6 km plateau there), not re-swept here. Operates
# on the GPS-derived ECI series (gps_eci) so the IC stays self-consistent
# with the ground truth it's later scored against.
# ==============================================================================

function _sma_corrected_r0_v0(t_cal::Float64; n_sma_avg::Int=5)
    idx0 = searchsortedfirst(gps_eci.t_rel, t_cal)
    idx0 = clamp(idx0, 1, nrow(gps_eci))
    n = min(n_sma_avg, nrow(gps_eci) - idx0 + 1)

    a_samples = Float64[]
    for k in idx0:(idx0 + n - 1)
        rk = SVector{3, Float64}(gps_eci.r_eci_0[k], gps_eci.r_eci_1[k], gps_eci.r_eci_2[k])
        vk = SVector{3, Float64}(gps_eci.v_eci_0[k], gps_eci.v_eci_1[k], gps_eci.v_eci_2[k])
        push!(a_samples, 1.0 / (2.0 / norm(rk) - dot(vk, vk) / planet.μ))
    end
    a_target = mean(a_samples)

    r0 = SVector{3, Float64}(gps_eci.r_eci_0[idx0], gps_eci.r_eci_1[idx0], gps_eci.r_eci_2[idx0])
    v0_raw = SVector{3, Float64}(gps_eci.v_eci_0[idx0], gps_eci.v_eci_1[idx0], gps_eci.v_eci_2[idx0])
    r0_mag = norm(r0)
    v0_mag = norm(v0_raw)
    v_target = sqrt(planet.μ * (2.0 / r0_mag - 1.0 / a_target))
    v_scale = v_target / v0_mag

    return r0, v0_raw .* v_scale
end

function build_cygnss_spacecraft(t_cal::Float64; id::Int64=1, sma_correct_ic::Bool=false)
    q0 = q_truth_scalarlast(t_cal)
    w0 = -1.0 .* SVector{3, Float64}(w1_itp(t_cal), w2_itp(t_cal), w3_itp(t_cal))
    r0, v0 = if sma_correct_ic
        _sma_corrected_r0_v0(t_cal)
    else
        SVector{3, Float64}(r1_itp(t_cal), r2_itp(t_cal), r3_itp(t_cal)), SVector{3, Float64}(v1_itp(t_cal), v2_itp(t_cal), v3_itp(t_cal))
    end
    ic = CartesianInitialCondition(r0, v0; q=q0, ang_vel=w0)

    main_bus = SimulationModel.Link{0}(
        root=true, m=BUS_MASS,
        dims=MVector{3, Float64}(BUS_DIMS...),
        ref_area=BUS_DIMS[1] * BUS_DIMS[3],
    )
    left_panel = SimulationModel.Link{0}(
        root=false, m=PANEL_MASS,
        dims=MVector{3, Float64}(PANEL_DIMS...),
        ref_area=PANEL_DIMS[2] * PANEL_DIMS[3],
        r=MVector{3, Float64}(0.0, -PANEL_OFFSET_Y, 0.0),
    )
    right_panel = SimulationModel.Link{0}(
        root=false, m=PANEL_MASS,
        dims=MVector{3, Float64}(PANEL_DIMS...),
        ref_area=PANEL_DIMS[2] * PANEL_DIMS[3],
        r=MVector{3, Float64}(0.0, PANEL_OFFSET_Y, 0.0),
    )

    return SimulationModel.SpacecraftModel(
        SimulationModel.Joint[],
        [main_bus, left_panel, right_panel],
        main_bus,
        true,
        main_bus.m + left_panel.m + right_panel.m,
        0.0,
        INERTIA_OVERRIDE,
        0,
        0,
        ic,
        id,
    )
end

# ==============================================================================
# Simulation runner
# ==============================================================================

function run_cygnss_case(;
    t_cal::Float64,
    t_end::Float64,
    dynamic_effectors_factory::Function, # t_offset::Float64 -> Tuple of dynamic effectors
    label::String,
    dt_max::Float64=0.5,
    density_model=NoAtmosphereModel(), # pass GRAMAtmosphereModel(planet_name="earth") when testing aero
    control_effectors_factory::Function=(t_offset) -> (WheelMomentumReplayControlModel(t_offset),),
    sma_correct_ic::Bool=false, # N=5 vis-viva SMA averaging + velocity-magnitude rescale (see build_cygnss_spacecraft)
)
    spacecraft = build_cygnss_spacecraft(t_cal; sma_correct_ic=sma_correct_ic)
    mission_time = t_end - t_cal
    dynamic_effectors = dynamic_effectors_factory(t_cal)
    control_effectors = control_effectors_factory(t_cal)

    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            results_directory=joinpath(STUDY_DIR, "output"),
            normalize=false,
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=mission_time,
            orientation_sim=true,
            num_steps_to_save=2000,
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false,
        ),
        dynamics_model=DynamicsModel([spacecraft], dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=fill(0.1, length(control_effectors))),
        initial_time=TELEMETRY_INITIAL_TIME,
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-10,
            abstol_orbit=1e-10,
            reltol_quaternion=1e-10,
            abstol_quaternion=1e-12,
            reltol_angular_rate=1e-10,
            abstol_angular_rate=1e-12,
            dt_max_orbit=dt_max,
        ),
    )

    sol = run_simulation(args; return_solution=true)

    ts = range(t_cal, t_end, length=300)
    errs = Float64[]
    pos_errs_km = Float64[]
    for t in ts
        u = sol(t - t_cal)
        q_sim = SVector{4, Float64}(u.sc[1].q[1], u.sc[1].q[2], u.sc[1].q[3], u.sc[1].q[4]) # scalar-last
        q_sim_sf = SVector{4, Float64}(q_sim[4], q_sim[1], q_sim[2], q_sim[3])
        q_gt = q_truth_scalarfirst(t)
        push!(errs, rad2deg(2 * acos(clamp(abs(dot(q_sim_sf, q_gt)), 0.0, 1.0))))

        pos_sim = SVector{3, Float64}(u.sc[1].pos)
        pos_gt = r_truth(t)
        push!(pos_errs_km, norm(pos_sim - pos_gt) / 1000.0)
    end

    println("[$label] t=$(t_cal)-$(t_end)s  attitude mean=$(mean(errs)) deg  max=$(maximum(errs)) deg  final=$(errs[end]) deg  |  position mean=$(mean(pos_errs_km)) km  max=$(maximum(pos_errs_km)) km  final=$(pos_errs_km[end]) km")
    return (
        label=label,
        mean=mean(errs), max=maximum(errs), final=errs[end], errs=errs,
        pos_mean_km=mean(pos_errs_km), pos_max_km=maximum(pos_errs_km), pos_final_km=pos_errs_km[end], pos_errs_km=pos_errs_km,
        sol=sol,
    )
end
