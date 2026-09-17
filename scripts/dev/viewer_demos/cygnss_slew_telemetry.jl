"""
    CygnssSlewTelemetry

The FM01 commanded-slew telemetry: its epoch, its attitude conventions, its
loaders, and the two derived angles the scenario reports.

Everything here that is a CONVENTION is a pure function of its arguments, so
`test/unit/mission/cygnss_slew_tests.jl` can pin it without the telemetry
files, which are gitignored. The loaders are separate and skip cleanly when the
files are absent.

The conventions, each with the evidence that settled it. Re-run all of them with
`julia --project=. scripts/dev/viewer_demos/cygnss_slew_checks.jl`.

## The epoch

The export's absolute time column counts seconds from the J2000 epoch
2000-01-01T12:00:00, not from midnight of that day; the two readings are twelve
hours apart. CYGNSS FM01's catalogue element set of 2025-06-12
(`cygnss_historical.tle`), propagated by SGP4 and rotated TEME to J2000, is
compared with the orbit plane the telemetry itself reports at its first sample.
The midnight reading leaves the right ascension of the ascending node 3.55 deg
and the inclination 0.032 deg away from the telemetered plane; the J2000-epoch
reading leaves them 0.28 deg and 0.005 deg away, after a 114-day propagation.
An order of magnitude on both, and the J2000-epoch reading is the one that
agrees with `docs/spaceagora_cygnss_reconstruction_record.md` section 2, which
dates the FM01 arc to 2025-10-04.

ASSUMPTION: the counter is read as UTC seconds. A counter kept in TAI, GPS or
barycentric dynamical time would move the epoch by 18 to 69 s. Nothing here can
resolve that, and nothing on the page depends on it.

## The attitude quaternion

The export's quaternion is left-handed, scalar-first and frame-to-body. The
mapping into the scalar-last quaternion SpaceAGORA integrates is

    (x, y, z, w) = (-t2, -t3, -t4, t1)

(`docs/spaceagora_cygnss_reconstruction_record.md`, section 2.2). The evidence:
`rot(q)` of the mapped quaternion carries the telemetered nadir direction into a
nearly constant body-frame vector, as a nadir-pointing spacecraft requires
(per-component standard deviation 0.002 to 0.079 over the hour), while the
unmapped ordering gives 0.34 to 0.39.

## The body rate

`omega_body = +w_eci`, in SpaceAGORA's convention and no other. This is NOT the
sign the standalone reconstruction in `extra_examples/` uses: that file
integrates the telemetry quaternion in its own convention, whose kinematics
carry the opposite sign on the cross-product term, and `-w_eci` is the body rate
there. Differentiating the mapped quaternion numerically and solving the
SpaceAGORA kinematics for the rate gives `+w_eci` to better than one percent at
every sample tested, including through the maneuver.

## The wheel-axis matrix

`data/telemetry/CYGNSS/cyg01_adcs_constants.toml` carries a wheel-axis matrix
regressed from flight telemetry against the OTHER body-rate sign, so in
SpaceAGORA's convention it enters NEGATED (see [`slew_wheel_axes`](@ref)). The
evidence is a momentum closure that uses no model at all: form
`H_total = C_bi(t)' (I omega(t) + A J_w Omega_rw(t))` from measured quantities
only and look at how constant it is in inertial space. The negated matrix holds
it to a per-axis standard deviation of 0.9 to 1.9e-4 N m s over the maneuver
window; the unnegated one gives 2.1 to 9.7e-4, two to five times worse.
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

export SLEW_TIME_ORIGIN, SLEW_COMMAND_STEP_S, SLEW_RPM_TO_RAD_S
export slew_epoch_utc, slew_quaternion_to_spaceagora, slew_scalar_first_attitude_matrix
export slew_body_rate, slew_wheel_axes
export attitude_angle_deg, lvlh_pointing_angle_deg
export load_slew_telemetry, load_slew_constants

"""
The instant the export's absolute time column counts seconds from: the J2000
epoch, 2000-01-01T12:00:00, read as UTC. See the module docstring for the test
that settled this against the alternative reading of midnight that day.
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
    slew_epoch_utc(t_abs_s) -> DateTime

UTC of an export time stamp. Leap seconds are not applied: the counter is read
as UTC seconds elapsed since the calendar instant [`SLEW_TIME_ORIGIN`], which is
the reading the module docstring states and marks as an assumption.
"""
slew_epoch_utc(t_abs_s::Real)::DateTime = SLEW_TIME_ORIGIN + Millisecond(round(Int, Float64(t_abs_s) * 1000))

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
    slew_wheel_axes(axes_from_file) -> Matrix{Float64}

The wheel-axis matrix in SpaceAGORA's convention: the constants file's matrix,
negated.

The file's matrix was regressed against the standalone reconstruction's body
rate, which is minus SpaceAGORA's. The conservation law the regression fitted,
`I omega + A J_w Omega_rw = constant`, is invariant under negating BOTH omega
and the wheel momentum, so the regression could not see which of the two signs
it had recovered; the kinematics can, and does. The forward torque law
`-dH/dt - omega x H` is not invariant, so this sign is not cosmetic: with the
matrix unnegated the same run reproduces the maneuver to 11 deg instead of
0.5 deg.
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

This is the maneuver in one number. It reads 0.09 deg before the commanded step
and about 10 deg after it. The triad is built from the state itself, so the
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

The FM01 ADCS export and its companion position/velocity table on one time base,
with the quaternion sign-unwrapped and converted into SpaceAGORA's convention
and the body rate signed into it.

Fields: `t_rel` (seconds from the first sample), `t_abs` (the export's own
clock), `q` (4 x N, SpaceAGORA convention), `q_telemetry` (4 x N, as exported,
unwrapped), `omega` (3 x N, body frame), `speeds_rad_s` (N x 3 wheel speeds),
`nadir_spread` (the convention check's own number), and the state columns on
their own time base: `pv_t_rel` (M), `pos_m` and `vel_mps` (3 x M), with
`pv_index` giving the 4 Hz row each fix was taken from.

The state columns are on their own time base because `adsAttFilter.pv_eci`
updates at 1 Hz while the attitude columns update at 4 Hz: three of every four
rows of the position and velocity are a repeat of the row before, 75.0% of the
file. Interpolating a staircase gives a curve that oscillates by kilometres
between the real fixes, so the repeats are dropped here and the 3602 distinct
fixes are returned as their own series. Nothing downstream should use the 4 Hz
position column.

`pv_t_rel` is also RE-TIMED. A fix inherits the time stamp of the 4 Hz packet
it first appears in, and those stamps wander: about one in ten is a packet late,
a few are half a second out, and over the hour they accumulate 0.76 s of drift
against the fixes themselves. A 0.25 s error is 1.9 km at orbital speed, which a
spline through the raw stamps turns into a kilometre-scale oscillation on either
side of a mislabeled fix.

The fixes themselves are on a clean uniform clock, and the data say so without
reference to any stamp: the chord between consecutive fixes divided by the speed
over it gives 1.000004 s for every one of the 3601 intervals, inside 0.0007 s.
So the interval is measured that way, the fixes are laid on a uniform grid of
it, and the grid's phase is the median of the raw stamps. (The chord understates
the arc by one part in 2e7 at this orbital rate, 0.2 ms over the whole hour,
which is far below the metre-level quantization of the position and velocity
columns.)

What this cannot fix is the grid's ABSOLUTE phase, which is known only to the
resolution of the stamps that anchor it, a fraction of a second, or about a
kilometre along track. Radial and cross-track comparisons against these fixes
are meaningful at the metre level; an along-track one is not, below a couple of
kilometres.

Raises rather than returning a silently mis-converted attitude: the body-frame
nadir direction must be nearly constant, or the quaternion convention does not
hold for this file and nothing downstream is meaningful.
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
    load_slew_constants(path) -> NamedTuple

The body inertia, wheel inertia and wheel-axis matrix from the gitignored
constants file, with the matrix already carried into SpaceAGORA's convention by
[`slew_wheel_axes`](@ref).

The file writes the wheel-axis matrix as three inner lists, which are its
COLUMNS: that is what makes their norms the 0.90 / 1.19 / 1.00 the file reports,
where reading them as rows gives 1.76 / 0.05 / 0.32. The function checks the
norms so the wrong reading cannot pass silently.
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
    return (inertia=inertia,
        wheel_inertia=Float64(cfg["wheels"]["inertia_kg_m2"]),
        wheel_axes=slew_wheel_axes(file_axes),
        wheel_axes_from_file=file_axes,
        momentum_rating_nms=Float64(cfg["wheels"]["momentum_rating_nms"]),
        speed_rating_rpm=Float64(cfg["wheels"]["speed_rating_rpm"]))
end

end # module CygnssSlewTelemetry
