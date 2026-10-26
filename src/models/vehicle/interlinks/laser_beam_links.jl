# The laser beam link: an open-cavity laser link with opt-in physics (line of
# sight past the planet, divergence beyond the beam waist, pointing within each
# terminal's field of regard, and recoil at the mounted terminals). With every
# option at its default it is the laser momentum link. Included by
# `InterLinkModels` (interlink_graph.jl), after links.jl and power_links.jl.

"""
    AbstractLaserBeam

Supertype of the beam models by which a `LaserBeamLink` scales its force. A beam
implements `beam_force_factor`. The built-in model is `GaussianBeam`.
"""
abstract type AbstractLaserBeam end

"""
    GaussianBeam(; wavelength_m=1.064e-6, waist_radius_m, aperture_radius_m)

A Gaussian beam for a `LaserBeamLink`, with its waist at the transmitting
terminal and a circular receiving aperture centred on the beam at the other
terminal.

- **The beam.** It has wavelength `wavelength_m`, and a waist of 1/e² radius
  `waist_radius_m`. Its spot radius at range `d` is `w0 √(1 + (d / z_R)²)`,
  with the Rayleigh range `z_R = π w0² / λ` (`beam_spot_radius_m`). Its
  divergence half-angle is `λ / (π w0)`.
- **The capture.** The receiving aperture has radius `aperture_radius_m`. It
  captures `1 − exp(−2 r² / w²)` of the beam (`beam_capture_fraction`), which is
  the beam's force factor.
- **Limits.** All three lengths must be finite and positive.
"""
struct GaussianBeam <: AbstractLaserBeam
    wavelength_m::Float64
    waist_radius_m::Float64
    aperture_radius_m::Float64

    function GaussianBeam(wavelength_m::Real, waist_radius_m::Real, aperture_radius_m::Real)
        isfinite(wavelength_m) && wavelength_m > 0 || throw(ArgumentError("wavelength_m must be finite and positive."))
        isfinite(waist_radius_m) && waist_radius_m > 0 || throw(ArgumentError("waist_radius_m must be finite and positive."))
        isfinite(aperture_radius_m) && aperture_radius_m > 0 ||
            throw(ArgumentError("aperture_radius_m must be finite and positive."))
        new(Float64(wavelength_m), Float64(waist_radius_m), Float64(aperture_radius_m))
    end
end

GaussianBeam(; wavelength_m::Real=1.064e-6, waist_radius_m::Real, aperture_radius_m::Real) =
    GaussianBeam(wavelength_m, waist_radius_m, aperture_radius_m)

function beam_spot_radius_m(beam::GaussianBeam, distance_m::Real)::Float64
    rayleigh_m = pi * beam.waist_radius_m^2 / beam.wavelength_m
    return beam.waist_radius_m * sqrt(1.0 + (Float64(distance_m) / rayleigh_m)^2)
end

function beam_capture_fraction(beam::GaussianBeam, radius_m::Real, distance_m::Real)::Float64
    spot = beam_spot_radius_m(beam, distance_m)
    return -expm1(-2.0 * Float64(radius_m)^2 / spot^2)
end

# The far-field divergence half-angle (rad) of a Gaussian beam: λ / (π w0).
_beam_divergence_rad(beam::GaussianBeam) = beam.wavelength_m / (pi * beam.waist_radius_m)

"""
    beam_force_factor(beam, link, terminal_a, terminal_b, distance_m) -> Float64

The share of a `LaserBeamLink`'s near-field force, `magnification × power_w / c`,
that `link` applies at a separation of `distance_m` (m) with `beam`, between
the terminal hardware `terminal_a` and `terminal_b` (`nothing` for a virtual
terminal).

- The engine calls it at every right-hand-side evaluation for each selected
  link with a beam, and when the orbital scheduling policies score one.
- A `GaussianBeam` gives the capture fraction at the far aperture,
  `beam_capture_fraction(beam, beam.aperture_radius_m, distance_m)`.
- A beam model of your own implements it: a recycling cavity, a measured
  curve, a non-Gaussian beam.
- Default: none. Every beam type implements it.
"""
function beam_force_factor end

beam_force_factor(beam::GaussianBeam, link, terminal_a, terminal_b, distance_m) =
    beam_capture_fraction(beam, beam.aperture_radius_m, distance_m)

"""
    TerminalPointing(; field_of_regard_rad, error_rad=0.0)

How a `LaserBeamLink`'s terminals point.

- **Field of regard.** Each terminal can steer its beam anywhere within
  `field_of_regard_rad` (at most π) of its boresight. It points there ideally,
  with no slew or gimbal state.
- **Pointing error.** A static error `error_rad` multiplies the force by the
  far-field Gaussian loss `exp(−2 (error / θ)²)`, with `θ = λ / (π w0)` the
  beam's divergence. It needs a `GaussianBeam`. Jitter, a random error, is not
  modelled.
"""
struct TerminalPointing
    field_of_regard_rad::Float64
    error_rad::Float64

    function TerminalPointing(field_of_regard_rad::Real, error_rad::Real)
        isfinite(field_of_regard_rad) && 0 < field_of_regard_rad <= pi ||
            throw(ArgumentError("field_of_regard_rad must be in (0, π]."))
        isfinite(error_rad) && error_rad >= 0 || throw(ArgumentError("error_rad must be finite and nonnegative."))
        new(Float64(field_of_regard_rad), Float64(error_rad))
    end
end

TerminalPointing(; field_of_regard_rad::Real, error_rad::Real=0.0) = TerminalPointing(field_of_regard_rad, error_rad)

"""
    LaserBeamLink(; power_w, magnification, range_m, beam=nothing, occultation=false, pointing=nothing, recoil=:central)
    LaserBeamLink(terminal::LaserTerminal; beam=nothing, occultation=false, pointing=nothing, recoil=:central)

An open-cavity laser link that pushes its two spacecraft apart, as
`LaserMomentumLink` does, with opt-in physics.

- **Defaults.** With every option at its default, it is `LaserMomentumLink`,
  bit for bit: a central force of `magnification * power_w / c` on each
  endpoint, directed away from the other, equal and opposite, at any separation
  up to `range_m`.
- **Requirements.** `power_w` (W) must be nonnegative, and `magnification` and
  `range_m` (m) positive. The terminal form takes them from a `LaserTerminal`:
  its `max_power_w`, `magnification` and `range_m`.
- **`occultation=true`.** The link is unavailable while the segment between the
  two spacecraft meets the planet, a sphere of its equatorial radius.
  `occultation=(; min_ray_altitude_m=h)` raises that sphere by `h` metres, so
  that a beam grazing the atmosphere counts as blocked. Attenuation is not
  modelled.
- **`beam`.** The force is scaled by `beam_force_factor(beam, link, terminal_a,
  terminal_b, d)` at separation `d`.
  - With a `GaussianBeam`, the factor is the capture fraction at the far
    aperture. It is close to 1 within the Rayleigh range, and far beyond it
    the force falls as `1/d²`.
  - This is an engineering model, not a cavity model: it takes the
    magnification as the near-field value and scales it by one pass's
    geometric capture.
- **`pointing=TerminalPointing(...)`.** The link is available only while each
  end's line of sight to the other lies within the field of regard of its
  terminal's boresight, in its attitude. A pointing error scales the force.
  This needs terminal hardware at both ends and attitude (`orientation_sim`).
- **`recoil=:mounted`.** Each end's force acts at its terminal's mounted
  location, with the torque of its lever arm, as the terminal's own recoil
  would (`actuator_wrench`). The default, `:central`, applies it with no
  torque. Mounted recoil needs terminal hardware at both ends and attitude.

With pointing or mounted recoil, range, occultation, divergence and force
use the ray between the mounted terminals, transformed by their attitudes.
Other configurations use the spacecraft body-reference positions. The
terminal-aware hooks are required for pointing or mounted recoil; the
four-argument `link_force` and `link_available` cannot supply that geometry.

The two forces stay equal and opposite, so the run's linear momentum is
conserved. While the link is active, the terminals' own recoil applies nothing:
the link's force is their recoil.
"""
struct LaserBeamLink{B <: Union{Nothing, AbstractLaserBeam}} <: AbstractInterlink
    power_w::Float64
    magnification::Float64
    range_m::Float64
    beam::B
    occultation::Bool
    min_ray_altitude_m::Float64
    pointing::Union{Nothing, TerminalPointing}
    recoil::Symbol

    function LaserBeamLink(power_w::Real, magnification::Real, range_m::Real, beam::B, occultation::Bool,
                           min_ray_altitude_m::Real, pointing::Union{Nothing, TerminalPointing},
                           recoil::Symbol) where {B <: Union{Nothing, AbstractLaserBeam}}
        isfinite(power_w) && power_w >= 0 || throw(ArgumentError("power_w must be finite and nonnegative."))
        isfinite(magnification) && magnification > 0 || throw(ArgumentError("magnification must be finite and positive."))
        isfinite(range_m) && range_m > 0 || throw(ArgumentError("range_m must be finite and positive."))
        isfinite(min_ray_altitude_m) && min_ray_altitude_m >= 0 ||
            throw(ArgumentError("min_ray_altitude_m must be finite and nonnegative."))
        recoil in (:central, :mounted) || throw(ArgumentError("recoil must be :central or :mounted; got $(repr(recoil))."))
        pointing === nothing || pointing.error_rad == 0 || beam isa GaussianBeam || throw(ArgumentError(
            "A pointing error scales the force by a Gaussian beam's divergence; give the link a GaussianBeam."))
        new{B}(Float64(power_w), Float64(magnification), Float64(range_m), beam, occultation,
            Float64(min_ray_altitude_m), pointing, recoil)
    end
end

function LaserBeamLink(; power_w::Real, magnification::Real, range_m::Real, beam=nothing, occultation=false,
                       pointing=nothing, recoil::Symbol=:central)
    blocked, altitude = _occultation_option(occultation)
    return LaserBeamLink(power_w, magnification, range_m, beam, blocked, altitude, pointing, recoil)
end

LaserBeamLink(terminal::LaserTerminal; kwargs...) =
    LaserBeamLink(; power_w=terminal.max_power_w, magnification=terminal.magnification, range_m=terminal.range_m, kwargs...)

_occultation_option(occultation::Bool) = (occultation, 0.0)
function _occultation_option(occultation::NamedTuple)
    keys(occultation) == (:min_ray_altitude_m,) || throw(ArgumentError(
        "occultation is true, false or (; min_ray_altitude_m=h); got $(repr(occultation))."))
    return (true, Float64(occultation.min_ray_altitude_m))
end
_occultation_option(occultation) = throw(ArgumentError(
    "occultation is true, false or (; min_ray_altitude_m=h); got $(repr(occultation))."))

link_range_m(link::LaserBeamLink) = link.range_m
link_terminal_levels(link::LaserBeamLink) = (link.power_w, link.power_w)
terminal_supports(::LaserTerminal, ::LaserBeamLink) = true
link_applies_torque(link::LaserBeamLink) = link.recoil === :mounted
_link_needs_terminals(link::LaserBeamLink) = link.pointing !== nothing || link.recoil === :mounted
_link_needs_attitude(link::LaserBeamLink) = link.pointing !== nothing || link.recoil === :mounted
_needs_terminals_reason(::LaserBeamLink) = " (its pointing or mounted recoil reads them)"

# The occultation test of power links, against the planet's sphere raised by
# `margin_m`: a beam that passes within `margin_m` of the surface is blocked.
_segment_occulted(p, q, planet, margin_m::Float64) =
    planet !== nothing && hasproperty(planet, :Rp_e) && _segment_occulted(p, q, (Rp_e=Float64(planet.Rp_e) + margin_m,))

# Pointing and mounted recoil share one physical ray. Other modes retain the
# body-reference geometry, including the bitwise-compatible all-default path.
function _beam_geometry(link::LaserBeamLink, terminal_a, terminal_b, a, b)
    separation = a.pos_ii - b.pos_ii
    _link_needs_attitude(link) || return (a.pos_ii, b.pos_ii, separation)
    offset_a = _beam_terminal_offset(a, terminal_a)
    offset_b = _beam_terminal_offset(b, terminal_b)
    # Subtract the body positions before adding small mount offsets, to avoid
    # losing their precision at large inertial coordinates.
    return (a.pos_ii + offset_a, b.pos_ii + offset_b, separation + (offset_a - offset_b))
end

function _beam_terminal_offset(sample, terminal)
    terminal === nothing && throw(ArgumentError(
        "A LaserBeamLink with pointing or mounted recoil needs terminal hardware at both ends; spacecraft $(sample.spacecraft_index)'s is virtual."))
    sample.q_ib === nothing && throw(ArgumentError(
        "A LaserBeamLink with pointing or mounted recoil needs the endpoints' attitude (orientation_sim)."))
    _, location = terminal_pose(sample.spacecraft, terminal)
    return transpose(rot(sample.q_ib)) * location
end

function link_available(link::LaserBeamLink, terminal_a, terminal_b, a, b, ctx)
    position_a, position_b, ray = _beam_geometry(link, terminal_a, terminal_b, a, b)
    separation = norm(ray)
    0 < separation <= link.range_m || return false
    link.occultation && _segment_occulted(position_a, position_b, ctx.planet, link.min_ray_altitude_m) && return false
    link.pointing === nothing && return true
    return _within_field_of_regard(link.pointing.field_of_regard_rad, a, terminal_a, -ray) &&
        _within_field_of_regard(link.pointing.field_of_regard_rad, b, terminal_b, ray)
end
link_available(link::LaserBeamLink, a, b, ctx) = link_available(link, nothing, nothing, a, b, ctx)

function interlink_wrench(link::LaserBeamLink, terminal_a, terminal_b, a, b, ctx)
    force_a, force_b = _beam_forces(link, terminal_a, terminal_b, a, b)
    link.recoil === :mounted || return (force_a, _ZERO_FORCE, force_b, _ZERO_FORCE)
    return (force_a, _mounted_torque(a, terminal_a, force_a), force_b, _mounted_torque(b, terminal_b, force_b))
end
link_force(link::LaserBeamLink, a, b, ctx) = _beam_forces(link, nothing, nothing, a, b)

# The two forces: the laser momentum link's, scaled by the beam's factor and the
# pointing loss when the link has them.
function _beam_forces(link::LaserBeamLink, terminal_a, terminal_b, a, b)
    _, _, separation = _beam_geometry(link, terminal_a, terminal_b, a, b)
    distance = norm(separation)
    distance == 0 && return (zero(separation), zero(separation))
    scale = link.magnification * link.power_w / SPEED_OF_LIGHT_MPS
    if link.beam !== nothing
        scale *= Float64(beam_force_factor(link.beam, link, terminal_a, terminal_b, distance))
        link.pointing === nothing || link.pointing.error_rad == 0 ||
            (scale *= exp(-2.0 * (link.pointing.error_rad / _beam_divergence_rad(link.beam))^2))
    end
    force = scale * (separation / distance)
    return (force, -force)
end

# Whether the terminal-to-terminal direction `sight_ii` lies within
# `field_of_regard` (rad) of its terminal's boresight, in its body frame.
function _within_field_of_regard(field_of_regard::Float64, sample, terminal, sight_ii)
    terminal === nothing && throw(ArgumentError(
        "A LaserBeamLink with pointing needs terminal hardware at both ends; spacecraft $(sample.spacecraft_index)'s is virtual."))
    sample.q_ib === nothing && throw(ArgumentError(
        "A LaserBeamLink with pointing needs the endpoints' attitude (orientation_sim)."))
    boresight, _ = terminal_pose(sample.spacecraft, terminal)
    sight = rot(sample.q_ib) * sight_ii
    return atan(norm(cross(sight, boresight)), dot(sight, boresight)) <= field_of_regard
end

# The body torque (N·m) of the inertial force `force_ii` applied at endpoint
# `sample`'s terminal: the terminal's own recoil lever arm (`terminal_pose`).
function _mounted_torque(sample, terminal, force_ii)
    terminal === nothing && throw(ArgumentError(
        "A LaserBeamLink with recoil=:mounted needs terminal hardware at both ends; spacecraft $(sample.spacecraft_index)'s is virtual."))
    sample.q_ib === nothing && throw(ArgumentError(
        "A LaserBeamLink with recoil=:mounted needs the endpoints' attitude (orientation_sim)."))
    _, location = terminal_pose(sample.spacecraft, terminal)
    return cross(location, rot(sample.q_ib) * force_ii)
end
