"""
    CygnssICs

Initial conditions for the CYGNSS constellation at a single epoch: the frame
conversions that produce them, the classical elements that sanity-check them,
and the loader for the JSON file `build_cygnss_ics.jl` writes.

The frame conversions here are deliberately pure linear algebra. They take a
rotation (and, for the Earth-fixed case, the angular velocity that goes with
it) as an argument rather than reaching for SPICE themselves, so they can be
exercised against closed-form cases in a unit test that needs no kernels and no
telemetry. `build_cygnss_ics.jl` supplies the SPICE-derived rotations.

Frames, named once:

  * **ITRF93** — the Earth-fixed frame the CYGNSS Level-1 PVT columns are
    reported in. `docs/spaceagora_cygnss_reconstruction_record.md`, section 2.1:
    "The PVT positions are reported in WGS84 Earth-fixed coordinates. They must
    be converted from ITRF93 to J2000 before they can be compared with a
    SpaceAGORA trajectory."
  * **TEME** — true equator, mean equinox of date, the frame SGP4 outputs. An
    element set propagated by SGP4 is in TEME by construction, never in J2000.
  * **J2000** — the mean equator and mean equinox of J2000.0, which is the
    inertial frame SpaceAGORA integrates in and the frame the viewer draws.

Everything this module returns is in meters and meters per second.
"""
module CygnssICs

using JSON
using LinearAlgebra
using StaticArrays

export CYGNSS_EPOCH_UTC, CYGNSS_MU_EARTH, CYGNSS_R_EARTH_EQ
export SpacecraftIC, ConstellationICs
export earth_fixed_to_inertial, rotating_to_inertial, teme_to_j2000
export classical_elements, angle_spread_deg, wrap_360
export load_constellation_ics, parse_tle_blocks, TleRecord

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

"The common epoch every state in the initial-conditions file is expressed at."
const CYGNSS_EPOCH_UTC = "2025-06-06T00:00:00Z"

"""
Earth's gravitational parameter, m^3/s^2. DE421 value, the same one
`SpaceAGORA`'s `Earth()` planet carries (`src/environment/ephemerides/planets.jl`),
so elements computed here agree with elements computed from a simulation state.
Used only for the sanity-check elements, never for the states themselves.
"""
const CYGNSS_MU_EARTH = 3.98600436233e14

"Earth's equatorial radius, m, matching `Earth().Rp_e`. Used for altitude only."
const CYGNSS_R_EARTH_EQ = 6.3781366e6

# ---------------------------------------------------------------------------
# Frame conversions
# ---------------------------------------------------------------------------

"""
    rotating_to_inertial(l_pi, omega_fixed, r_fixed, v_fixed) -> (r_inertial, v_inertial)

State in a uniformly rotating body-fixed frame to the inertial frame, given
`l_pi`, the inertial-to-body-fixed rotation (so `l_pi'` maps body-fixed to
inertial), and `omega_fixed`, the body's angular velocity **expressed in the
body-fixed frame**.

The velocity term is `l_pi' * (v_fixed + omega_fixed x r_fixed)`, that is, the
transport term is formed in the body-fixed frame and rotated out with the
position. Forming it instead as `omega x r_inertial` with the same numerical
`omega` treats a body-fixed vector as an inertial one; for Earth in 2025 the
two spin axes differ by the accumulated precession since J2000 (about 0.35
degrees), which at CYGNSS altitude is worth about 3 m/s. That is not a rounding
difference, so the distinction is kept explicit here.

This is the uniform-rotation approximation. `earth_fixed_to_inertial` below
takes the exact rotation derivative instead and is what the generator uses.
"""
function rotating_to_inertial(
    l_pi::AbstractMatrix,
    omega_fixed::AbstractVector,
    r_fixed::AbstractVector,
    v_fixed::AbstractVector,
)
    rf = SVector{3, Float64}(r_fixed)
    vf = SVector{3, Float64}(v_fixed)
    w = SVector{3, Float64}(omega_fixed)
    R = SMatrix{3, 3, Float64}(l_pi)
    r_i = R' * rf
    v_i = R' * (vf + cross(w, rf))
    return (r_i, v_i)
end

"""
    earth_fixed_to_inertial(l_pi, dl_pi_dt, r_fixed, v_fixed) -> (r_inertial, v_inertial)

State in an Earth-fixed frame to the inertial frame using the exact rotation
derivative, which is what a full state transform (SPICE `sxform`) applies:

    r_i = L' r_f
    v_i = L' v_f + (dL/dt)' r_f

`l_pi` is the inertial-to-fixed rotation and `dl_pi_dt` its time derivative in
the same sense. Carrying the derivative rather than a constant spin rate picks
up polar motion and length-of-day, which the uniform-rotation form above misses.
"""
function earth_fixed_to_inertial(
    l_pi::AbstractMatrix,
    dl_pi_dt::AbstractMatrix,
    r_fixed::AbstractVector,
    v_fixed::AbstractVector,
)
    rf = SVector{3, Float64}(r_fixed)
    vf = SVector{3, Float64}(v_fixed)
    R = SMatrix{3, 3, Float64}(l_pi)
    dR = SMatrix{3, 3, Float64}(dl_pi_dt)
    return (R' * rf, R' * vf + dR' * rf)
end

"""
    teme_to_j2000(r_teme_to_j2000, r_teme, v_teme) -> (r_j2000, v_j2000)

Rotate an SGP4 output state from TEME into J2000. `r_teme_to_j2000` is the
rotation matrix for that pair at the state's epoch.

Both frames are quasi-inertial, so a single rotation carries position and
velocity alike: there is no transport term. The residual rate of the
TEME-to-J2000 rotation is precession and nutation, of order 5e-12 rad/s, which
over the seconds-scale validity of any one state is far below the kilometers of
error an element set already carries.
"""
function teme_to_j2000(r_teme_to_j2000::AbstractMatrix, r_teme::AbstractVector, v_teme::AbstractVector)
    R = SMatrix{3, 3, Float64}(r_teme_to_j2000)
    return (R * SVector{3, Float64}(r_teme), R * SVector{3, Float64}(v_teme))
end

# ---------------------------------------------------------------------------
# Classical elements, for sanity checks and for the constellation geometry
# ---------------------------------------------------------------------------

wrap_360(x::Real) = mod(float(x), 360.0)

"""
    classical_elements(r, v; mu=CYGNSS_MU_EARTH, r_body=CYGNSS_R_EARTH_EQ) -> NamedTuple

Classical elements of an inertial state, in meters, seconds and degrees:
`(a_m, e, inclination_deg, raan_deg, arg_perigee_deg, true_anomaly_deg,
arg_latitude_deg, period_s, altitude_m, apogee_altitude_m, perigee_altitude_m)`.

`arg_latitude_deg` is the angle from the ascending node to the spacecraft,
measured in the orbit plane. It is computed directly from the node direction
rather than as `arg_perigee + true_anomaly`, because CYGNSS eccentricity is
about 0.001: at that eccentricity the perigee direction is poorly determined
and the sum of two noisy angles is much worse conditioned than the one angle
that actually matters for where a spacecraft sits in its plane.

`altitude_m` is above the equatorial radius, so it is a scalar summary and not
a geodetic altitude.
"""
function classical_elements(
    r::AbstractVector,
    v::AbstractVector;
    mu::Real = CYGNSS_MU_EARTH,
    r_body::Real = CYGNSS_R_EARTH_EQ,
)
    rv = SVector{3, Float64}(r)
    vv = SVector{3, Float64}(v)
    rn = norm(rv)
    vn = norm(vv)
    h = cross(rv, vv)
    hn = norm(h)
    # Node vector: z-hat x h, pointing at the ascending node.
    n = SVector{3, Float64}(-h[2], h[1], 0.0)
    nn = norm(n)

    e_vec = (cross(vv, h) ./ mu) .- (rv ./ rn)
    e = norm(e_vec)

    a = 1.0 / (2.0 / rn - vn^2 / mu)
    incl = acosd(clamp(h[3] / hn, -1.0, 1.0))

    raan = nn > 0 ? wrap_360(atand(n[2], n[1])) : 0.0

    argp = if nn > 0 && e > 0
        w = acosd(clamp(dot(n, e_vec) / (nn * e), -1.0, 1.0))
        e_vec[3] < 0 ? 360.0 - w : w
    else
        0.0
    end

    nu = if e > 0
        f = acosd(clamp(dot(e_vec, rv) / (e * rn), -1.0, 1.0))
        dot(rv, vv) < 0 ? 360.0 - f : f
    else
        0.0
    end

    # Argument of latitude straight from the node, valid at any eccentricity.
    u = if nn > 0
        nh = n ./ nn
        wh = cross(h ./ hn, nh)          # in-plane, 90 deg ahead of the node
        wrap_360(atand(dot(rv, wh), dot(rv, nh)))
    else
        0.0
    end

    period = 2π * sqrt(a^3 / mu)
    return (
        a_m = a,
        e = e,
        inclination_deg = incl,
        raan_deg = raan,
        arg_perigee_deg = argp,
        true_anomaly_deg = nu,
        arg_latitude_deg = u,
        period_s = period,
        altitude_m = rn - r_body,
        apogee_altitude_m = a * (1 + e) - r_body,
        perigee_altitude_m = a * (1 - e) - r_body,
    )
end

"""
    angle_spread_deg(angles_deg) -> NamedTuple

Spread of a set of angles on the circle: `(min_deg, max_deg, span_deg,
largest_gap_deg)`. `span_deg` is the width of the smallest arc containing every
angle, found by sorting and taking the complement of the largest gap between
neighbors; a plain `maximum - minimum` would report 359 degrees for two
spacecraft one degree apart across the wrap.
"""
function angle_spread_deg(angles_deg::AbstractVector{<:Real})
    a = sort(wrap_360.(collect(float.(angles_deg))))
    n = length(a)
    n == 0 && return (min_deg = NaN, max_deg = NaN, span_deg = NaN, largest_gap_deg = NaN)
    n == 1 && return (min_deg = a[1], max_deg = a[1], span_deg = 0.0, largest_gap_deg = 360.0)
    gaps = [a[i + 1] - a[i] for i in 1:(n - 1)]
    push!(gaps, a[1] + 360.0 - a[end])
    biggest = maximum(gaps)
    return (
        min_deg = a[1],
        max_deg = a[end],
        span_deg = 360.0 - biggest,
        largest_gap_deg = biggest,
    )
end

# ---------------------------------------------------------------------------
# The initial-conditions file
# ---------------------------------------------------------------------------

"One spacecraft's entry in the initial-conditions file."
struct SpacecraftIC
    name::String
    norad_id::Union{Int, Nothing}
    r_ii_m::SVector{3, Float64}
    v_ii_m_s::SVector{3, Float64}
    provenance::String
    source::String
    epoch_offset_s::Float64
end

"The whole file: a common epoch, a named frame, and the spacecraft in it."
struct ConstellationICs
    epoch_utc::String
    frame::String
    spacecraft::Vector{SpacecraftIC}
    notes::String
end

# The British spelling of "catalogue" is the interface, fixed by the schema the
# viewer agent reads; it stays as written so the two sides agree on the literal.
# Prose in this repository uses American English (CLAUDE.md).
const VALID_PROVENANCE = ("telemetry", "catalogue", "nominal")

"""
    load_constellation_ics(path) -> ConstellationICs

Read and validate the initial-conditions JSON. Raises rather than returning a
half-built object: a missing key, a state that is not three numbers, or a
`provenance` outside `telemetry` / `catalogue` / `nominal` is an error, because
a viewer page that silently drew a spacecraft with no provenance would be
exactly the thing the file exists to prevent.
"""
function load_constellation_ics(path::AbstractString)
    isfile(path) || throw(ArgumentError("initial-conditions file not found: $path"))
    raw = JSON.parsefile(String(path))
    for key in ("epoch_utc", "frame", "spacecraft")
        haskey(raw, key) || throw(ArgumentError("initial-conditions file is missing \"$key\": $path"))
    end
    entries = raw["spacecraft"]
    entries isa AbstractVector || throw(ArgumentError("\"spacecraft\" must be an array: $path"))
    isempty(entries) && throw(ArgumentError("\"spacecraft\" is empty: $path"))

    sats = SpacecraftIC[]
    for (i, e) in enumerate(entries)
        for key in ("name", "r_ii_m", "v_ii_m_s", "provenance", "source", "epoch_offset_s")
            haskey(e, key) || throw(ArgumentError("spacecraft $i is missing \"$key\": $path"))
        end
        r = e["r_ii_m"]
        v = e["v_ii_m_s"]
        (r isa AbstractVector && length(r) == 3) ||
            throw(ArgumentError("spacecraft $i \"r_ii_m\" must be three numbers: $path"))
        (v isa AbstractVector && length(v) == 3) ||
            throw(ArgumentError("spacecraft $i \"v_ii_m_s\" must be three numbers: $path"))
        prov = String(e["provenance"])
        prov in VALID_PROVENANCE ||
            throw(ArgumentError("spacecraft $i has provenance \"$prov\", expected one of $(VALID_PROVENANCE): $path"))
        nid = get(e, "norad_id", nothing)
        push!(sats, SpacecraftIC(
            String(e["name"]),
            nid === nothing ? nothing : Int(nid),
            SVector{3, Float64}(Float64.(r)),
            SVector{3, Float64}(Float64.(v)),
            prov,
            String(e["source"]),
            Float64(e["epoch_offset_s"]),
        ))
    end
    return ConstellationICs(
        String(raw["epoch_utc"]),
        String(raw["frame"]),
        sats,
        String(get(raw, "notes", "")),
    )
end

# ---------------------------------------------------------------------------
# Two-line element sets
# ---------------------------------------------------------------------------

"One element set out of a snapshot file, with the snapshot it came from."
struct TleRecord
    block::String        # which capture the set came from
    name::String
    line1::String
    line2::String
end

"""
    parse_tle_blocks(path) -> Dict{String, Vector{TleRecord}}

Read a TLE file that is divided into named captures by `BLOCK <name> ...`
lines, ignoring `#` comments and blank lines. Returns each block's element sets
keyed by the block label.

The block structure exists because more than one archived capture brackets the
target epoch, and the generator scores them against telemetry rather than
picking one on faith.
"""
function parse_tle_blocks(path::AbstractString)
    isfile(path) || throw(ArgumentError("TLE file not found: $path"))
    blocks = Dict{String, Vector{TleRecord}}()
    current = "default"
    pending = String[]
    for raw in eachline(String(path))
        line = rstrip(raw)
        (isempty(line) || startswith(lstrip(line), "#")) && continue
        if startswith(line, "BLOCK ")
            current = String(split(line)[2])
            get!(blocks, current, TleRecord[])
            empty!(pending)
            continue
        end
        push!(pending, String(line))
        if length(pending) == 3
            startswith(pending[2], "1 ") && startswith(pending[3], "2 ") ||
                throw(ArgumentError("malformed three-line element set near \"$(pending[1])\" in $path"))
            push!(get!(blocks, current, TleRecord[]),
                  TleRecord(current, strip(pending[1]), pending[2], pending[3]))
            empty!(pending)
        end
    end
    isempty(pending) ||
        throw(ArgumentError("trailing partial element set ($(length(pending)) lines) in $path"))
    return blocks
end

end # module
