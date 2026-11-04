"""
CYGNSS private telemetry conventions and loader. See
`docs/spaceagora_cygnss_reconstruction_record.md` for provenance and limits.

The flight counter is read as ET, converted to UTC through SPICE. Quaternion
mapping is scalar-first telemetry to scalar-last `(-q1,-q2,-q3,q0)`; body rates
use the positive recorded convention. The signed wheel mapping and inertia are
supplied separately. These inputs carry unresolved scale/polarity qualifications;
retaining a mapping is not an independent physical validation of it.

Run the convention and physical-adequacy checks through `scripts/dev/run.jl`.
Keep their measured results and all input-derived artifacts private.
"""
module CygnssSlewTelemetry

using Arrow
using DataFrames
using Dates
using LinearAlgebra
using StaticArrays
using Statistics
using TOML
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SPICE = SpaceAGORA.TelemetryVerification.SPICE

export SLEW_TIME_ORIGIN, SLEW_COMMAND_STEP_S, SLEW_RPM_TO_RAD_S
export slew_epoch_et, slew_epoch_utc, slew_quaternion_to_spaceagora, slew_scalar_first_attitude_matrix
export slew_body_rate, slew_wheel_axes
export attitude_angle_deg, lvlh_pointing_angle_deg, momentum_ledger
export load_slew_telemetry, load_slew_constants

"""
The instant the export's absolute time column counts seconds from: the J2000
epoch, 2000-01-01T12:00:00 TDB (SPICE's ET origin), which is
2000-01-01T11:58:55.816 UTC. The calendar value here names the instant; the
counter is ET, so it is converted with the leap-second kernel, never by adding
seconds to this DateTime. See the module docstring for the evidence.
"""
const SLEW_TIME_ORIGIN = DateTime(2000, 1, 1, 12, 0, 0)

"""
Where `acsLvlhPoint.qwTgt.q`, the flight software's LVLH pointing target, steps
from identity to the commanded rotation, in seconds from the first sample of the
export.
"""
const SLEW_COMMAND_STEP_S = 901.75

const SLEW_RPM_TO_RAD_S = 2pi / 60

"""
    slew_epoch_et(t_abs_s) -> Float64

SPICE ephemeris time (TDB seconds past J2000) of an export time stamp. The
export's counter already is that quantity, so this is the identity, kept as a
named function so every consumer states the reading it relies on.
"""
slew_epoch_et(t_abs_s::Real)::Float64 = Float64(t_abs_s)

"""
    slew_epoch_utc(t_abs_s) -> DateTime

UTC of an export time stamp, to the millisecond, through SPICE's `et2utc` and
the loaded leap-second kernel. Requires a leap-second kernel to be furnished
(the scenario and the checks load the starter pack through `Earth`); without
one SPICE raises, which is the intended failure rather than a silent
69-second offset.
"""
function slew_epoch_utc(t_abs_s::Real)::DateTime
    stamp = lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
        SPICE.et2utc(slew_epoch_et(t_abs_s), "ISOC", 3)
    end
    return DateTime(stamp, dateformat"yyyy-mm-ddTHH:MM:SS.sss")
end

"""
    slew_quaternion_to_spaceagora(q_telemetry) -> SVector{4, Float64}

The export's scalar-first quaternion `(t1, t2, t3, t4)` as the scalar-last
quaternion SpaceAGORA integrates: `(-t2, -t3, -t4, t1)`.
"""
@inline slew_quaternion_to_spaceagora(q::AbstractVector{<:Real})::SVector{4, Float64} =
    SVector{4, Float64}(-Float64(q[2]), -Float64(q[3]), -Float64(q[4]), Float64(q[1]))

"""
    slew_scalar_first_attitude_matrix(q) -> SMatrix{3, 3, Float64}

The attitude matrix of a scalar-first quaternion under the standard convention,
`(qs^2 - |qv|^2) I + 2 qv qv' - 2 qs [qv x]`. It exists so a test can state the
identity the mapping above relies on: `rot(slew_quaternion_to_spaceagora(q))` is
the TRANSPOSE of this, because negating the vector part conjugates the
quaternion and conjugating transposes the matrix. The export's own frame
direction is therefore the transposed one, which is what the nadir check
measures.
"""
function slew_scalar_first_attitude_matrix(q::AbstractVector{<:Real})::SMatrix{3, 3, Float64}
    qs = Float64(q[1])
    qv = SVector{3, Float64}(Float64(q[2]), Float64(q[3]), Float64(q[4]))
    skew = SMatrix{3, 3, Float64}(0.0, qv[3], -qv[2], -qv[3], 0.0, qv[1], qv[2], -qv[1], 0.0)
    return (qs^2 - dot(qv, qv)) * SMatrix{3, 3, Float64}(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0) +
        2 * qv * qv' - 2 * qs * skew
end

"""
    slew_body_rate(w_eci) -> SVector{3, Float64}

The body-frame angular velocity SpaceAGORA integrates, from the export's
`adsAttFilter.qw_eci.w` columns. It is `+w_eci`; the sign lives here, once, with
the evidence in the module docstring.
"""
@inline slew_body_rate(w::AbstractVector{<:Real})::SVector{3, Float64} =
    SVector{3, Float64}(Float64(w[1]), Float64(w[2]), Float64(w[3]))

"""
    slew_wheel_axes(matrix)

Apply the retained constants-file sign convention. This preserves the original
comparison configuration; provenance and unresolved export-polarity questions
belong to the private reconstruction record. Do not infer a new wheel polarity
from a better fit on the same arc.
"""
slew_wheel_axes(axes_from_file::AbstractMatrix{<:Real})::Matrix{Float64} = -Matrix{Float64}(axes_from_file)

"""
    attitude_angle_deg(qa, qb) -> Float64

The rotation angle between two scalar-last attitudes, in degrees, insensitive to
the sign a quaternion is written with.
"""
@inline function attitude_angle_deg(qa::AbstractVector{<:Real}, qb::AbstractVector{<:Real})::Float64
    a = SVector{4, Float64}(qa[1], qa[2], qa[3], qa[4])
    b = SVector{4, Float64}(qb[1], qb[2], qb[3], qb[4])
    return rad2deg(2.0 * acos(clamp(abs(dot(a / norm(a), b / norm(b))), 0.0, 1.0)))
end

"""
    lvlh_pointing_angle_deg(q, r, v) -> Float64

The angle between the vehicle's attitude and the local-vertical/local-horizontal
attitude it holds when it is on its pointing target: body +z on nadir, body +y
on the negative orbit normal, body +x completing the triad along track.

The triad is built from the state itself, so the
function needs no convention beyond the attitude one: `-r` is nadir and `r x v`
is the orbit normal.
"""
function lvlh_pointing_angle_deg(q::AbstractVector{<:Real}, r::AbstractVector{<:Real}, v::AbstractVector{<:Real})::Float64
    rr = SVector{3, Float64}(r[1], r[2], r[3])
    vv = SVector{3, Float64}(v[1], v[2], v[3])
    z_lvlh = -rr / norm(rr)
    h = cross(rr, vv)
    y_lvlh = -h / norm(h)
    x_lvlh = cross(y_lvlh, z_lvlh)
    qq = SVector{4, Float64}(q[1], q[2], q[3], q[4])
    # The pointing error is the rotation carrying the LVLH triad onto the body
    # triad; its matrix is (inertial -> body) applied to the LVLH axes written
    # as columns in inertial, and its angle follows from the trace.
    c_err = SM.rot(qq / norm(qq)) * hcat(x_lvlh, y_lvlh, z_lvlh)
    return rad2deg(acos(clamp((tr(c_err) - 1.0) / 2.0, -1.0, 1.0)))
end

"""
    load_slew_telemetry(adcs_path, pv_path) -> NamedTuple

Load the matching private packet-time tables. Return t_rel, absolute t_abs,
sign-unwrapped q/q_telemetry, body omega, N-by-wheel speeds_rad_s, and the
separate navigation series pv_t_rel, pv_index, pos_m, vel_mps. Repeated fixes
are dropped and retimed using chord/speed cadence and median packet phase.
That interpolation does not resolve absolute along-track phase uncertainty.
Reject inconsistent packet grids, irregular fixes and a failed nadir convention
check. The data and computed spread remain private.
"""
function load_slew_telemetry(adcs_path::AbstractString, pv_path::AbstractString)
    isfile(adcs_path) || throw(ArgumentError("no FM01 ADCS export at $(adcs_path)."))
    isfile(pv_path) || throw(ArgumentError("no FM01 position/velocity export at $(pv_path)."))
    adcs = DataFrame(Arrow.Table(String(adcs_path)))
    pv = DataFrame(Arrow.Table(String(pv_path)))
    t = Float64.(adcs.t_rel)
    issorted(t) || throw(ArgumentError("$(adcs_path): t_rel is not sorted."))
    Float64.(pv.t_rel) == t || throw(ArgumentError("$(adcs_path) and $(pv_path) are not on the same time base."))

    n = length(t)
    q_tel = Matrix{Float64}(undef, 4, n)
    q_sa = Matrix{Float64}(undef, 4, n)
    @inbounds for i in 1:n
        qi = SVector{4, Float64}(adcs.q_eci_0[i], adcs.q_eci_1[i], adcs.q_eci_2[i], adcs.q_eci_3[i])
        # q and -q are the same attitude and an onboard filter flips sign
        # freely; unwrap so an interpolated curve sees no false discontinuity.
        if i > 1 && dot(qi, SVector{4, Float64}(q_tel[1, i - 1], q_tel[2, i - 1], q_tel[3, i - 1], q_tel[4, i - 1])) < 0.0
            qi = -qi
        end
        q_tel[:, i] .= qi
        q_sa[:, i] .= slew_quaternion_to_spaceagora(qi)
    end
    omega = Matrix{Float64}(undef, 3, n)
    @inbounds for i in 1:n
        omega[:, i] .= slew_body_rate(SVector{3, Float64}(adcs.w_eci_0[i], adcs.w_eci_1[i], adcs.w_eci_2[i]))
    end
    pos_all = permutedims(hcat(Float64.(pv.r_eci_0), Float64.(pv.r_eci_1), Float64.(pv.r_eci_2)))
    vel_all = permutedims(hcat(Float64.(pv.v_eci_0), Float64.(pv.v_eci_1), Float64.(pv.v_eci_2)))
    speeds = hcat(Float64.(adcs.Omega_rw_0), Float64.(adcs.Omega_rw_1), Float64.(adcs.Omega_rw_2)) .* SLEW_RPM_TO_RAD_S

    # Drop the repeated rows of the 1 Hz state channel; see the docstring.
    keep = Int[1]
    @inbounds for i in 2:n
        if pos_all[1, i] != pos_all[1, i - 1] || pos_all[2, i] != pos_all[2, i - 1] || pos_all[3, i] != pos_all[3, i - 1]
            push!(keep, i)
        end
    end
    length(keep) >= 3 || throw(ArgumentError("$(pv_path): fewer than three distinct navigation fixes."))
    pos = pos_all[:, keep]
    vel = vel_all[:, keep]
    pv_t = _retime_navigation_fixes(t[keep], pos, vel, pv_path)

    # The nadir check is evaluated on the distinct fixes, where the position is
    # the one the filter actually reported.
    nadir_body = Matrix{Float64}(undef, 3, length(keep))
    @inbounds for (j, i) in enumerate(keep)
        r = SVector{3, Float64}(pos_all[1, i], pos_all[2, i], pos_all[3, i])
        nadir_body[:, j] .= SM.rot(SVector{4, Float64}(q_sa[:, i])) * (-r / norm(r))
    end
    spread = maximum(std(nadir_body[k, :]) for k in 1:3)
    spread < 0.15 || throw(ErrorException(
        "the telemetry quaternion mapping does not hold for $(adcs_path): the body-frame nadir direction varies " *
        "by $(spread), which is not a nadir-pointing spacecraft."))

    return (t_rel=t, t_abs=Float64.(adcs.t), q=q_sa, q_telemetry=q_tel, omega=omega,
        pv_t_rel=pv_t, pv_index=keep, pos_m=pos, vel_mps=vel, speeds_rad_s=speeds, nadir_spread=spread)
end

"""
    _retime_navigation_fixes(t_raw, pos, vel, source) -> Vector{Float64}

Put the navigation fixes back on the uniform clock they were taken on. See
[`load_slew_telemetry`](@ref) for why the raw stamps cannot be used.

The interval is measured from the data rather than read off the stamps: over one
fix the chord `|r_{k+1} - r_k|` divided by the mean speed is the elapsed time to
about a fifth of a millisecond, since position and velocity are both reported to
better than a metre and a metre per second. The fixes then go on a uniform grid
of that interval, phased by the median of the raw stamps.

Raises when any measured interval is more than two percent from the median one,
which would mean a fix is missing and a uniform grid is the wrong repair, or
when the phased grid ends up more than one interval from a raw stamp.
"""
function _retime_navigation_fixes(
    t_raw::AbstractVector{<:Real},
    pos::AbstractMatrix{<:Real},
    vel::AbstractMatrix{<:Real},
    source::AbstractString,
)::Vector{Float64}
    t = Vector{Float64}(t_raw)
    n = length(t)
    n >= 3 || return t
    intervals = Vector{Float64}(undef, n - 1)
    @inbounds for k in 1:(n - 1)
        chord = norm(SVector{3, Float64}(pos[1, k + 1] - pos[1, k], pos[2, k + 1] - pos[2, k], pos[3, k + 1] - pos[3, k]))
        speed = 0.5 * (norm(SVector{3, Float64}(vel[:, k])) + norm(SVector{3, Float64}(vel[:, k + 1])))
        speed > 0.0 || throw(ArgumentError("$(source): a navigation fix reports zero speed."))
        intervals[k] = chord / speed
    end
    cadence = median(intervals)
    cadence > 0.0 || throw(ArgumentError("$(source): the navigation fixes have a non-positive cadence."))
    worst_interval = maximum(abs.(intervals .- cadence)) / cadence
    worst_interval <= 0.02 || throw(ArgumentError(
        "$(source): a measured fix interval is $(round(100 * worst_interval; digits=2))% from the median " *
        "$(round(cadence; digits=6)) s; the fixes are not on one uniform clock and cannot be re-timed onto a grid."))
    grid = [(k - 1) * cadence for k in 1:n]
    phase = median(t .- grid)
    retimed = grid .+ phase
    worst_stamp = maximum(abs.(t .- retimed))
    worst_stamp <= cadence || throw(ArgumentError(
        "$(source): a raw time stamp is $(round(worst_stamp; digits=4)) s from its grid point, more than one " *
        "$(round(cadence; digits=6)) s interval; the stamps and the fixes disagree too much to anchor the grid."))
    return retimed
end

"""
    momentum_ledger(tel, constants; window, rod_duty=nothing) -> NamedTuple

The external torque the flown vehicle felt over `window = (t_lo, t_hi)` in
`t_rel` seconds, measured rather than modeled. Form the total angular momentum
in inertial space from measured quantities only,

    H_total(t) = C_bi(t)' ( I omega(t) + A J_w Omega_rw(t) ),

at the navigation fixes; its linear drift over the window IS the net external
torque, gravity gradient and torque rods together, by definition. Returned
fields: `drift_nm` (3-vector, N m), `gravity_gradient_nm` (the mean
gravity-gradient torque over the window, 3-vector), `omitted_nms` (the drift's
magnitude times the window length: the external momentum a replay carrying no
external torque leaves out), `exchanged_nms` (the swing of the wheel-momentum
magnitude over the window: what such a replay does carry), `ratio` (omitted
over exchanged) and `rod_duty` (per-rod mean |duty| when `rod_duty`, an N x 3
matrix on `tel.t_rel`, is given; otherwise `nothing`). The ratio is the margin
of a wheels-only replay: above one, the omitted torque moves more momentum
over the window than the wheels do.
"""
function momentum_ledger(tel, constants; window::Tuple{Real, Real}, rod_duty=nothing)
    lo, hi = Float64(window[1]), Float64(window[2])
    lo < hi || throw(ArgumentError("momentum_ledger: the window must be increasing."))
    mu = 3.986004418e14
    idx = findall(x -> lo <= x <= hi, tel.pv_t_rel)
    length(idx) >= 3 || throw(ArgumentError("momentum_ledger: fewer than three navigation fixes in the window."))
    m = length(idx)
    h = Matrix{Float64}(undef, 3, m)
    tau = Matrix{Float64}(undef, 3, m)
    hw_norm = Vector{Float64}(undef, m)
    for (j, jj) in enumerate(idx)
        i = tel.pv_index[jj]
        c_bi = SM.rot(SVector{4, Float64}(tel.q[:, i]))
        omega = SVector{3, Float64}(tel.omega[:, i])
        hw = SVector{3, Float64}(constants.wheel_axes * (constants.wheel_inertia .* SVector{3, Float64}(tel.speeds_rad_s[i, :])))
        hw_norm[j] = norm(hw)
        h[:, j] .= c_bi' * (constants.inertia * omega + hw)
        r_body = c_bi * SVector{3, Float64}(tel.pos_m[:, jj])
        rr = norm(r_body); rhat = r_body / rr
        tau[:, j] .= c_bi' * (3 * mu / rr^3 * cross(rhat, constants.inertia * rhat))
    end
    tt = tel.pv_t_rel[idx]
    tc = tt .- mean(tt)
    drift = SVector{3, Float64}([sum(tc .* (h[k, :] .- mean(h[k, :]))) / sum(tc .^ 2) for k in 1:3])
    gg = SVector{3, Float64}([mean(tau[k, :]) for k in 1:3])
    omitted = norm(drift) * (hi - lo)
    exchanged = maximum(hw_norm) - minimum(hw_norm)
    duty = if rod_duty === nothing
        nothing
    else
        wi = findall(x -> lo <= x <= hi, tel.t_rel)
        SVector{3, Float64}([mean(abs.(Float64.(rod_duty[wi, k]))) for k in 1:3])
    end
    return (window=(lo, hi), drift_nm=drift, gravity_gradient_nm=gg, omitted_nms=omitted,
        exchanged_nms=exchanged, ratio=omitted / exchanged, rod_duty=duty, fixes=m)
end

"""
    load_slew_constants(path)

Read restricted body inertia and signed wheel axes (inner lists are columns).
Optional effective_inertia_scale factors reproduce the supplied configuration;
they are not validated by this loader. Later research questions those factors,
so report their provenance and keep the physical-adequacy result visible.
"""
function load_slew_constants(path::AbstractString)
    isfile(path) || throw(ArgumentError("no ADCS constants at $(path)."))
    cfg = TOML.parsefile(String(path))
    inertia = SMatrix{3, 3, Float64}(reduce(hcat, [Float64.(row) for row in cfg["body"]["inertia_kg_m2"]]))
    isapprox(inertia, inertia'; atol=1e-12) || throw(ArgumentError("$(path): the body inertia is not symmetric."))
    file_axes = Matrix{Float64}(reduce(hcat, [Float64.(row) for row in cfg["wheels"]["wheel_axes"]]))
    norms = [norm(file_axes[:, k]) for k in 1:size(file_axes, 2)]
    all(0.5 .< norms .< 2.0) || throw(ArgumentError(
        "$(path): the wheel-axis columns have norms $(norms). The file's inner lists are the matrix COLUMNS; " *
        "read as rows they give norms far from the unit spin axes expected."))
    # Optional per-wheel effective inertia, relative to the datasheet inertia
    # (a supplied reconstruction assumption, not validated here). It
    # scales the axis columns, so the momentum model's single wheel inertia
    # stays the datasheet value and H_w = axes * (inertia .* speeds) carries
    # the calibration per wheel. Absent means ones.
    scale = haskey(cfg["wheels"], "effective_inertia_scale") ?
        Float64.(cfg["wheels"]["effective_inertia_scale"]) : ones(size(file_axes, 2))
    length(scale) == size(file_axes, 2) && all(isfinite, scale) && all(>(0.0), scale) || throw(ArgumentError(
        "$(path): effective_inertia_scale must give one finite positive factor per wheel."))
    return (inertia=inertia,
        wheel_inertia=Float64(cfg["wheels"]["inertia_kg_m2"]),
        wheel_axes=slew_wheel_axes(file_axes) .* permutedims(scale),
        wheel_axes_from_file=file_axes,
        effective_inertia_scale=scale,
        momentum_rating_nms=Float64(cfg["wheels"]["momentum_rating_nms"]),
        speed_rating_rpm=Float64(cfg["wheels"]["speed_rating_rpm"]))
end

end # module CygnssSlewTelemetry
