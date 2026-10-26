# Power links: a beam from a transmitter terminal to a `PowerReceiver`, whose
# delivered power is generation in the receiver's power budget. Included by
# `InterLinkModels` (interlink_graph.jl), after links.jl.

"""
    link_received_power_w(link, receiver, a::StateSample, b::StateSample, ctx::LinkContext) -> Float64

The electrical power (W) that `receiver`, the terminal hardware at `link`'s
receiving end `b`, delivers to its spacecraft's bus while `link` is active, with
`a` the transmitting end. It is generation in that spacecraft's power budget.
The engine asks for it at every right-hand-side evaluation, for each active
directed link that ends on terminal hardware, and `PowerDeliveryPolicy` and
`EnergyBalancingPolicy` score candidate links by it. Default: 0.
"""
link_received_power_w(link, receiver, a, b, ctx) = 0.0

"""
    PowerBeamLink(; power_w, wavelength_m=1.064e-6, waist_radius_m=0.05, range_m)

A power beam from a transmitter terminal (a `LaserTerminal`) to a
`PowerReceiver`: a directed link that carries power and applies no force.

- **The beam.** The transmitter runs at `power_w` optical watts while the link
  is active, and draws its own electrical power for it (`power_w / efficiency`
  for a `LaserTerminal`), as a tier-1 laser load. The beam is Gaussian, with
  its waist, of 1/e² radius `waist_radius_m`, at the transmitter, and
  wavelength `wavelength_m`. Its spot radius at range `d` is
  `beam_spot_radius_m(link, d)`.
- **The receiver.** It captures `beam_capture_fraction(link, radius, d)` of the
  beam through its aperture of radius `radius`. It delivers that times its
  `efficiency`, up to its `max_power_w`.
- **Availability.** The link can close within `range_m`, while the segment
  between the two spacecraft clears the planet (a sphere of its equatorial
  radius).
- **Deferred.** It does not model pointing, the receiver's incidence angle, the
  atmosphere, beam quality or the beam's radiation pressure.
"""
struct PowerBeamLink <: AbstractInterlink
    power_w::Float64
    wavelength_m::Float64
    waist_radius_m::Float64
    range_m::Float64

    function PowerBeamLink(power_w::Real, wavelength_m::Real, waist_radius_m::Real, range_m::Real)
        isfinite(power_w) && power_w >= 0 || throw(ArgumentError("power_w must be finite and nonnegative."))
        isfinite(wavelength_m) && wavelength_m > 0 || throw(ArgumentError("wavelength_m must be finite and positive."))
        isfinite(waist_radius_m) && waist_radius_m > 0 || throw(ArgumentError("waist_radius_m must be finite and positive."))
        isfinite(range_m) && range_m > 0 || throw(ArgumentError("range_m must be finite and positive."))
        new(Float64(power_w), Float64(wavelength_m), Float64(waist_radius_m), Float64(range_m))
    end
end

PowerBeamLink(; power_w::Real, wavelength_m::Real=1.064e-6, waist_radius_m::Real=0.05, range_m::Real) =
    PowerBeamLink(power_w, wavelength_m, waist_radius_m, range_m)

link_range_m(link::PowerBeamLink) = link.range_m
link_directed(::PowerBeamLink) = true
link_terminal_levels(link::PowerBeamLink) = (link.power_w, 0.0)

terminal_supports(::LaserTerminal, ::PowerBeamLink) = true
terminal_supports(::LaserTerminal, ::PowerBeamLink, side::Symbol) = side === :transmitter
terminal_supports(::PowerReceiver, link) = false
terminal_supports(::PowerReceiver, ::PowerBeamLink) = true
terminal_supports(::PowerReceiver, ::PowerBeamLink, side::Symbol) = side === :receiver
link_terminal(::PowerReceiver) = true

"""
    beam_spot_radius_m(link::PowerBeamLink, distance_m) -> Float64
    beam_spot_radius_m(beam::GaussianBeam, distance_m) -> Float64

The 1/e² radius (m) of `link`'s beam (or of `beam`) at `distance_m` from its waist:
`w0 √(1 + (d / z_R)²)`, with the Rayleigh range `z_R = π w0² / λ`.
"""
function beam_spot_radius_m(link::PowerBeamLink, distance_m::Real)::Float64
    rayleigh_m = pi * link.waist_radius_m^2 / link.wavelength_m
    return link.waist_radius_m * sqrt(1.0 + (Float64(distance_m) / rayleigh_m)^2)
end

"""
    beam_capture_fraction(link::PowerBeamLink, radius_m, distance_m) -> Float64
    beam_capture_fraction(beam::GaussianBeam, radius_m, distance_m) -> Float64

The fraction of `link`'s beam power (or `beam`'s) within `radius_m` of its axis at
`distance_m`, the share a centred circular aperture of that radius captures:
`1 − exp(−2 r² / w²)`, with `w` the spot radius there.
"""
function beam_capture_fraction(link::PowerBeamLink, radius_m::Real, distance_m::Real)::Float64
    spot = beam_spot_radius_m(link, distance_m)
    return -expm1(-2.0 * Float64(radius_m)^2 / spot^2)
end

function link_available(link::PowerBeamLink, a, b, ctx)
    separation = norm(a.pos_ii - b.pos_ii)
    0 < separation <= link.range_m || return false
    return !_segment_occulted(a.pos_ii, b.pos_ii, ctx.planet)
end

function link_received_power_w(link::PowerBeamLink, receiver::PowerReceiver, a, b, ctx)
    distance = norm(a.pos_ii - b.pos_ii)
    captured = link.power_w * beam_capture_fraction(link, receiver.radius_m, distance)
    return min(receiver.efficiency * captured, receiver.max_power_w)
end

# Whether the segment from `p` to `q` (inertial, from the planet's centre) meets
# the planet, a sphere of its equatorial radius; a segment that only touches it
# counts as blocked. Without a planet, nothing blocks.
function _segment_occulted(p, q, planet)::Bool
    (planet === nothing || !hasproperty(planet, :Rp_e)) && return false
    direction = q - p
    span = dot(direction, direction)
    t = span == 0 ? 0.0 : clamp(-dot(p, direction) / span, 0.0, 1.0)
    return norm(p + t * direction) <= Float64(planet.Rp_e)
end
