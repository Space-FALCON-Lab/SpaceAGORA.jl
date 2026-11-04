"""
    CygnssTracks

The flown CYGNSS constellation over the 96-hour telemetry window: fitting a
state out of a sampled position arc, comparing two states in the orbit frame,
and reading back the track table `build_cygnss_ics.jl` writes.

This module is the testable half of that script. Like `CygnssICs`, everything
here is a pure function of its arguments: no SPICE, no file layout knowledge
beyond the one loader at the bottom, so each piece can be exercised against a
closed-form case with neither kernels nor telemetry present.

Meters, seconds and degrees throughout.
"""
module CygnssTracks

using Arrow
using DataFrames
using LinearAlgebra
using StaticArrays

export fit_state_from_arc, rtn_offset, orbit_plane_angle_deg, along_track_time_s
export ConstellationTracks, load_constellation_tracks, track_names, track_state_at

"""
    reference_sample_indices(times, stop_s; max_samples=2500)

Choose evenly spaced indices in an ordered navigation track up to `stop_s`,
retaining both endpoints. Decimation must not shorten the reference coverage.
"""
function reference_sample_indices(times::AbstractVector, stop_s::Real; max_samples::Int=2500)
    max_samples >= 2 || throw(ArgumentError("reference budget must be at least two samples"))
    isfinite(stop_s) || throw(ArgumentError("reference end time must be finite"))
    last = searchsortedlast(times, stop_s)
    last == 0 && return Int[]
    return unique!(round.(Int, range(1, last; length=min(last, max_samples))))
end

# ---------------------------------------------------------------------------
# Fitting a state out of a sampled position arc
# ---------------------------------------------------------------------------

"""
    fit_state_from_arc(t_s, positions, t_eval; order=5) -> (r, v, residual_rms_m)

Least-squares polynomial fit of degree `order` to a sampled position arc,
evaluated at `t_eval`: the position is the polynomial's value and the velocity
its derivative.

`positions` is a vector of three-vectors sampled at the times `t_s`, which need
not be evenly spaced. Each axis is fitted independently; `residual_rms_m` is the
root-mean-square of the fit residual over all three axes and all samples, which
is the number that says whether the arc was smooth enough for the fit to mean
anything.

Why fit at all, rather than read a velocity column: the CYGNSS Level 1 position
and velocity variables are `Int32` in meters and meters per second, so a single
velocity sample is quantized at 1 m/s, while the same instant fitted from three
hundred position samples moves by only a few tenths of a meter per second as the
arc length is varied. `docs/spaceagora_cygnss_reconstruction_record.md`
section 2.1 reaches the same conclusion and fits its own initial states "from
position over an arc rather than taken from a single quantized velocity sample".
The improvement is a few times, not orders of magnitude: what limits the fit is
the choice of arc and order, not the quantization it averages down, and
`residual_rms_m` is how that choice is checked.

The times are shifted to `t_eval` and then scaled by the largest offset in the
arc before the fit is formed. The shift makes the constant term the answer and the
linear term the velocity; the scaling is what keeps a high-order fit usable, and
it is not optional. Over a 600-second arc in raw seconds the degree-7 Vandermonde
column `tau^7` spans nineteen decades, the normal equations lose every digit, and
the "fit" comes back megameters away from the data. In the scaled variable every
column is bounded by one.
"""
function fit_state_from_arc(
    t_s::AbstractVector{<:Real},
    positions::AbstractVector{<:AbstractVector{<:Real}},
    t_eval::Real;
    order::Int = 5,
)
    n = length(t_s)
    n == length(positions) ||
        throw(ArgumentError("got $n times and $(length(positions)) positions"))
    n > order ||
        throw(ArgumentError("a degree-$order fit needs more than $order samples, got $n"))

    tau = [float(t) - float(t_eval) for t in t_s]
    scale = maximum(abs, tau)
    scale > 0 || throw(ArgumentError("every sample sits at t_eval; the arc has no span"))
    u = tau ./ scale
    A = Matrix{Float64}(undef, n, order + 1)
    for p in 0:order, i in 1:n
        A[i, p + 1] = u[i]^p
    end

    r = zeros(3)
    v = zeros(3)
    squares = 0.0
    for c in 1:3
        y = [Float64(p[c]) for p in positions]
        coeffs = A \ y
        r[c] = coeffs[1]                               # value at t_eval
        v[c] = order >= 1 ? coeffs[2] / scale : 0.0    # d/dt at t_eval, unscaled
        squares += sum(abs2, A * coeffs .- y)
    end
    return (SVector{3, Float64}(r), SVector{3, Float64}(v), sqrt(squares / (3 * n)))
end

# ---------------------------------------------------------------------------
# Comparing two states
# ---------------------------------------------------------------------------

"""
    rtn_offset(r_ref, v_ref, dr) -> (radial_m, along_m, cross_m)

A position difference resolved in the reference state's own orbit frame: radial
outward, along-track in the direction of motion, cross-track along the orbit
normal.

The along-track axis is `cross(normal, radial)` rather than the velocity
direction itself. For a near-circular orbit the two differ by only the flight
path angle, but the triad built this way is exactly orthonormal, so the three
components add up in quadrature to the total and no error leaks between axes.
"""
function rtn_offset(r_ref::AbstractVector, v_ref::AbstractVector, dr::AbstractVector)
    rv = SVector{3, Float64}(r_ref)
    vv = SVector{3, Float64}(v_ref)
    d = SVector{3, Float64}(dr)
    rhat = normalize(rv)
    nhat = normalize(cross(rv, vv))
    that = cross(nhat, rhat)
    return (radial_m = dot(d, rhat), along_m = dot(d, that), cross_m = dot(d, nhat))
end

"""
    orbit_plane_angle_deg(r_a, v_a, r_b, v_b) -> degrees

Angle between two states' orbit planes, as the angle between their
angular-momentum vectors. This is the single number that covers both an
inclination difference and a difference in right ascension of the ascending
node, and it is what a cross-track error measures.

The answer is an arccosine of a dot product that is close to one, so it carries
about half the available digits: the resolution floor is around a millionth of a
degree, and two identical states return that rather than exactly zero. Angles
below that are noise; the catalog-versus-flight plane differences this is used
for are thousandths of a degree, four orders above it.
"""
function orbit_plane_angle_deg(r_a, v_a, r_b, v_b)
    ha = normalize(cross(SVector{3, Float64}(r_a), SVector{3, Float64}(v_a)))
    hb = normalize(cross(SVector{3, Float64}(r_b), SVector{3, Float64}(v_b)))
    return acosd(clamp(dot(ha, hb), -1.0, 1.0))
end

"""
    along_track_time_s(along_m, speed_m_s) -> seconds

An along-track offset read as a timing error: how far ahead or behind its own
track a spacecraft is. For a catalog state this is usually the honest way to
state the error, because the error is a phase error and not a shape error.
"""
along_track_time_s(along_m::Real, speed_m_s::Real) = float(along_m) / float(speed_m_s)

# ---------------------------------------------------------------------------
# The track table
# ---------------------------------------------------------------------------

"""
The flight solution of the whole constellation over one window: the common
epoch, and one table row per spacecraft per sample. Positions and velocities are
J2000 Earth-centered inertial in meters, `t_s` is seconds from the epoch.
"""
struct ConstellationTracks
    epoch_utc::String
    frame::String
    table::DataFrame
end

const TRACK_COLUMNS = (
    "name", "norad_id", "time",
    "pos_ii_1", "pos_ii_2", "pos_ii_3",
    "vel_ii_1", "vel_ii_2", "vel_ii_3",
)

"""
    load_constellation_tracks(path; epoch_utc, frame) -> ConstellationTracks

Read the track table and check that it carries the columns a reference ghost
needs, that each spacecraft's samples are sorted in time, and that nothing is
missing. A half-valid table would draw a ghost that silently jumps, which is the
one failure a viewer cannot show you.

The table `build_cygnss_ics.jl` writes carries more than the required columns:
the Earth-fixed `sc_pos_*_pvt_m` and `sc_vel_*_pvt_mps` the archive reported, as
`Int32` because that is the type and the quantization the archive really has, and
a Boolean `arc_consistent`. The last one is false on the roughly two samples in a
thousand whose position disagrees with its own immediate neighbors by more than
40 m, which is several times the 9 m a one-second second difference should show
at this altitude. Those are the product's bad navigation fixes; a page drawing a
ghost should skip them rather than connect them.
"""
function load_constellation_tracks(
    path::AbstractString;
    epoch_utc::AbstractString = "2025-06-06T00:00:00Z",
    frame::AbstractString = "J2000 Earth-centered inertial, meters and meters per second",
)
    isfile(path) || throw(ArgumentError("track table not found: $path"))
    df = DataFrame(Arrow.Table(String(path)))
    for col in TRACK_COLUMNS
        hasproperty(df, Symbol(col)) ||
            throw(ArgumentError("track table is missing column \"$col\": $path"))
    end
    isempty(df) && throw(ArgumentError("track table is empty: $path"))
    for name in unique(df.name)
        t = df.time[df.name .== name]
        issorted(t) || throw(ArgumentError("track for \"$name\" is not sorted in time: $path"))
        any(isnan, t) && throw(ArgumentError("track for \"$name\" has a missing time: $path"))
    end
    return ConstellationTracks(String(epoch_utc), String(frame), df)
end

"The spacecraft in a track table, in the order they first appear."
track_names(tracks::ConstellationTracks) = unique(tracks.table.name)

"""
    track_state_at(tracks, name, t_s) -> (r, v)

The flown state of one spacecraft at an arbitrary time, by cubic Hermite
interpolation between the two bracketing samples using the samples' own
velocities.

Hermite rather than linear: the samples are one second apart and the orbit turns
through about 0.06 degrees in that second, so linear interpolation would cut the
corner by roughly 2 centimeters at the midpoint. That is small, but the
interpolation costs nothing more and the result is then smooth in velocity as
well as in position, which is what a ghost drawn at a simulation's own irregular
save times needs.

Times outside the track are an error rather than a clamp: a ghost that quietly
freezes at the end of its data looks exactly like a ghost that agrees.
"""
function track_state_at(tracks::ConstellationTracks, name::AbstractString, t_s::Real)
    df = tracks.table
    rows = findall(==(String(name)), df.name)
    isempty(rows) && throw(ArgumentError("no track for \"$name\""))
    t = float(t_s)
    times = @view df.time[rows]
    (t >= first(times) && t <= last(times)) ||
        throw(ArgumentError("t = $t s is outside the track for \"$name\" " *
                            "($(first(times)) .. $(last(times)) s)"))
    j = searchsortedlast(times, t)
    j == length(times) && (j -= 1)
    i0, i1 = rows[j], rows[j + 1]
    t0, t1 = df.time[i0], df.time[i1]
    h = t1 - t0
    r0 = SVector(df.pos_ii_1[i0], df.pos_ii_2[i0], df.pos_ii_3[i0])
    r1 = SVector(df.pos_ii_1[i1], df.pos_ii_2[i1], df.pos_ii_3[i1])
    v0 = SVector(df.vel_ii_1[i0], df.vel_ii_2[i0], df.vel_ii_3[i0])
    v1 = SVector(df.vel_ii_1[i1], df.vel_ii_2[i1], df.vel_ii_3[i1])
    h == 0 && return (r0, v0)
    s = (t - t0) / h
    h00 = 2s^3 - 3s^2 + 1
    h10 = s^3 - 2s^2 + s
    h01 = -2s^3 + 3s^2
    h11 = s^3 - s^2
    r = h00 .* r0 .+ (h10 * h) .* v0 .+ h01 .* r1 .+ (h11 * h) .* v1
    # Derivative of the same polynomial, so the returned velocity is consistent
    # with the returned position rather than a second, independent guess.
    g00 = (6s^2 - 6s) / h
    g10 = 3s^2 - 4s + 1
    g01 = (-6s^2 + 6s) / h
    g11 = 3s^2 - 2s
    v = g00 .* r0 .+ g10 .* v0 .+ g01 .* r1 .+ g11 .* v1
    return (r, v)
end

end # module
