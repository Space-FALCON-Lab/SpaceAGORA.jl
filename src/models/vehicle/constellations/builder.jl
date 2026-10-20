# Constellation layouts: the orbital planes of a constellation, the Walker
# patterns, and the builder that fills each slot with a fresh spacecraft.
# A constellation is spacecraft laid out in orbits together; links between
# them are interlinks (src/models/vehicle/interlinks/), a separate model.
# Included by `SpacecraftModels` (model.jl), after `SpacecraftModel`.

"""
    OrbitalPlane(; a, e=0.0, i, Ω, ω=0.0, M)

One orbital plane of a constellation: the semimajor axis `a` (m), eccentricity
`e`, inclination `i`, RAAN `Ω` and argument of periapsis `ω` (degrees) its
spacecraft share, and the mean anomaly (degrees) of each of its slots, in slot
order, as `M`. Lay out a constellation from planes with `build_constellation`;
`walker_planes` gives the planes of a Walker pattern.
"""
struct OrbitalPlane
    a::Float64
    e::Float64
    i::Float64
    Ω::Float64
    ω::Float64
    M::Vector{Float64}

    function OrbitalPlane(a::Real, e::Real, i::Real, Ω::Real, ω::Real, M)
        isfinite(a) && a > 0 || throw(ArgumentError("An orbital plane's semimajor axis a must be finite and positive; got $a."))
        isfinite(e) && 0 <= e < 1 || throw(ArgumentError("An orbital plane's eccentricity e must be in [0, 1); got $e."))
        all(isfinite, (i, Ω, ω)) || throw(ArgumentError("An orbital plane's angles i, Ω and ω must be finite."))
        anomalies = collect(Float64, M)
        isempty(anomalies) && throw(ArgumentError("An orbital plane needs at least one slot: M lists the mean anomaly of each."))
        all(isfinite, anomalies) || throw(ArgumentError("An orbital plane's mean anomalies M must be finite."))
        return new(Float64(a), Float64(e), Float64(i), Float64(Ω), Float64(ω), anomalies)
    end
end

OrbitalPlane(; a::Real, e::Real=0.0, i::Real, Ω::Real, ω::Real=0.0, M) = OrbitalPlane(a, e, i, Ω, ω, M)

Base.:(==)(x::OrbitalPlane, y::OrbitalPlane) =
    (x.a, x.e, x.i, x.Ω, x.ω) == (y.a, y.e, y.i, y.Ω, y.ω) && x.M == y.M

"""
    walker_planes(; total, planes, phasing, a=nothing, altitude=nothing, planet=nothing,
                  e=0.0, i, ω=0.0, raan0=0.0, M0=0.0, pattern=:delta) -> Vector{OrbitalPlane}

The planes of the Walker constellation `i: total/planes/phasing`: `total`
spacecraft in `planes` equally populated planes, all at inclination `i`
(degrees).

- Plane `j` (from 0) has RAAN `raan0 + j * spread / planes`, where the spread is
  360° for `pattern=:delta` and 180° for `pattern=:star`.
- Slot `k` (from 0) of plane `j` has mean anomaly
  `M0 + k * 360 / s + j * phasing * 360 / total` (degrees, reduced to
  [0, 360)), with `s = total / planes` slots per plane. `phasing` is an integer
  from 0 to `planes - 1`.
- Give the semimajor axis `a` (m), or the `altitude` (m) above the equatorial
  radius `Rp_e` of `planet`.

`e` and `ω` (degrees) are shared by every plane.
"""
function walker_planes(; total::Integer, planes::Integer, phasing::Integer,
                       a::Union{Nothing, Real}=nothing, altitude::Union{Nothing, Real}=nothing, planet=nothing,
                       e::Real=0.0, i::Real, ω::Real=0.0, raan0::Real=0.0, M0::Real=0.0, pattern::Symbol=:delta)
    total > 0 && planes > 0 || throw(ArgumentError("A Walker constellation needs positive total and planes; got $total and $planes."))
    total % planes == 0 || throw(ArgumentError(
        "A Walker constellation spreads its spacecraft equally over its planes: planes ($planes) must divide total ($total)."))
    0 <= phasing < planes || throw(ArgumentError("A Walker phasing is an integer from 0 to planes - 1 ($(planes - 1)); got $phasing."))
    pattern in (:delta, :star) || throw(ArgumentError("A Walker pattern is :delta or :star; got :$pattern."))
    semimajor = _walker_semimajor_axis(a, altitude, planet)
    spread = pattern === :delta ? 360.0 : 180.0
    per_plane = total ÷ planes
    return [OrbitalPlane(semimajor, e, i, raan0 + j * spread / planes, ω,
                [mod(M0 + k * 360.0 / per_plane + j * phasing * 360.0 / total, 360.0) for k in 0:per_plane - 1])
            for j in 0:planes - 1]
end

function _walker_semimajor_axis(a, altitude, planet)::Float64
    (a === nothing) != (altitude === nothing) || throw(ArgumentError(
        "A Walker constellation takes the semimajor axis a or the altitude (with its planet), not both."))
    a === nothing || return Float64(a)
    planet !== nothing && hasproperty(planet, :Rp_e) || throw(ArgumentError(
        "A Walker altitude is measured above the planet's equatorial radius Rp_e: pass the planet too."))
    return Float64(planet.Rp_e) + Float64(altitude)
end

"""
    build_constellation(factory, planes; first_id=1) -> Vector{SpacecraftModel}

The spacecraft of a constellation laid out in `planes` (`OrbitalPlane`s), plane
by plane and slot by slot. For each slot it calls `factory(plane, slot, index)`,
with the plane and slot numbered from 1 and `index` the spacecraft's position in
the returned vector. The factory builds and returns a new `SpacecraftModel`; the
builder then sets on it:

- the orbit of its slot: an `InitialCondition` with the plane's elements and the
  slot's mean anomaly, keeping the attitude `q` and body rate `ang_vel` of the
  factory's own initial condition;
- `id = first_id + index - 1`;
- `constellation_slot = ConstellationSlot(plane, slot)`.

Spacecraft, links, mounted devices and many GNC effectors hold state, so the
factory must build them anew on every call. A factory that returns the same
spacecraft twice, or shares a joint, link, mounted device, attachment or
mutable GNC effector between two spacecraft, is refused.
"""
function build_constellation(factory, planes::AbstractVector{OrbitalPlane}; first_id::Integer=1)
    isempty(planes) && throw(ArgumentError("build_constellation needs at least one orbital plane."))
    spacecraft = SpacecraftModel[]
    slots = ConstellationSlot[]
    owners = IdDict{Any, Int}()
    for (p, plane) in enumerate(planes), (s, mean_anomaly) in enumerate(plane.M)
        index = length(spacecraft) + 1
        model = factory(p, s, index)
        model isa SpacecraftModel || throw(ArgumentError(
            "build_constellation: the factory returned a $(typeof(model)) for plane $p, slot $s; it must return a SpacecraftModel."))
        slot = ConstellationSlot(p, s)
        _claim_spacecraft_parts!(owners, model, index, slot, slots)
        q, ang_vel = _initial_attitude_of(model.initial_condition)
        model.initial_condition = InitialCondition(a=plane.a, e=plane.e, i=plane.i, Ω=plane.Ω, ω=plane.ω,
            M=mean_anomaly, q=q, ang_vel=ang_vel)
        model.id = first_id + index - 1
        model.constellation_slot = slot
        push!(spacecraft, model)
        push!(slots, slot)
    end
    return spacecraft
end

"""
    walker_constellation(factory; total, planes, phasing, first_id=1, kwargs...) -> Vector{SpacecraftModel}

The spacecraft of a Walker constellation: `build_constellation(factory,
walker_planes(; total, planes, phasing, kwargs...); first_id)`. For example, a
Walker delta 53°: 24/3/1 at 550 km:

```julia
factory(plane, slot, index) = begin
    bus = Link(root=true, m=150.0)
    SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia)
end
spacecraft = walker_constellation(factory; total=24, planes=3, phasing=1,
    altitude=550e3, planet=planet, i=53)
```
"""
walker_constellation(factory; first_id::Integer=1, kwargs...) =
    build_constellation(factory, walker_planes(; kwargs...); first_id)

_initial_attitude_of(ic) = hasproperty(ic, :q) && hasproperty(ic, :ang_vel) ?
    (ic.q, ic.ang_vel) : (DEFAULT_INITIAL_CONDITION_Q, DEFAULT_INITIAL_CONDITION_ANG_VEL)

# Record every mutable part of `model` as held by spacecraft `index`, refusing a
# part another spacecraft already holds: the spacecraft, its joints, its links
# with their devices and device lists, its attachments and its GNC effectors.
# Immutable values hold no state of their own and may be shared.
function _claim_spacecraft_parts!(owners::IdDict{Any, Int}, model::SpacecraftModel, index::Int,
                                  slot::ConstellationSlot, slots::Vector{ConstellationSlot})
    claim!(part) = ismutable(part) && begin
        holder = get(owners, part, 0)
        if holder != 0 && holder != index
            throw(ArgumentError(
                "build_constellation: the spacecraft for plane $(slot.plane), slot $(slot.slot) shares a " *
                "$(nameof(typeof(part))) with the one for plane $(slots[holder].plane), slot $(slots[holder].slot). " *
                "Spacecraft hold their own state, so the factory must build every part anew on each call."))
        end
        owners[part] = index
    end
    claim!(model)
    foreach(claim!, model.joints)
    for link in unique(objectid, [model.root; model.links])
        claim!(link)
        for name in fieldnames(Link)
            value = getfield(link, name)
            value isa StaticArray && continue
            claim!(value)
            value isa AbstractVector && foreach(claim!, value)
        end
    end
    foreach(claim!, model.attachments)
    for effectors in (model.guidance.guidance_effectors, model.navigation.navigation_effectors,
                      model.control.control_effectors)
        foreach(claim!, effectors)
    end
    return owners
end
