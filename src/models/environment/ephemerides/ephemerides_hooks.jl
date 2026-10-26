# The ephemerides hooks a model implements, model first.

"""
    ephemerides_frame_lpi(model, planet, et) -> SMatrix{3, 3, Float64}

The rotation from the inertial (J2000) frame to `planet`'s planet-fixed frame at
ephemeris time `et`, so that `r_pf = L * r_ii`. The engine calls it wherever it
needs the planet-fixed frame, for example for the atmosphere, gravity-harmonics
and magnetic-field samples, guidance and the scene export. It returns the
rotation only; a caller that needs planet-relative velocity adds the transport
term with the planet's spin vector.

`et` is in seconds past the J2000 epoch on the model's own timeline: the run's
start, [`ephemerides_epoch_seconds`](@ref), plus the elapsed time.
`SpiceEphemeridesModel` asks SPICE for the transform; `SimpleEphemeridesModel`
turns the planet about its spin axis, by default Earth by its Greenwich mean
sidereal time.
"""
function ephemerides_frame_lpi end

"""
    ephemerides_sun_ii(model, planet, et) -> Union{Nothing, SVector{3, Float64}}

Unit vector from the planet's center to the Sun, in the inertial (J2000) frame
the integrated state uses, or `nothing` when this ephemeris backend cannot
resolve the Sun for this body. `et` is as for
[`ephemerides_frame_lpi`](@ref). The power callbacks and the saved Sun
direction call it.

`SpiceEphemeridesModel` asks SPICE for the Sun's position relative to the
planet, so it works for every body whose SPK is furnished. `SimpleEphemeridesModel`
carries no planetary ephemeris, so it answers only for Earth, from a
low-precision analytic Sun; it returns `nothing` for any other body rather than
inventing a direction.
"""
function ephemerides_sun_ii end

"""
    ephemerides_epoch_seconds(model, initial_time) -> Float64

The run's start epoch, `initial_time` (an `InitialTime`, UTC), in seconds past
the J2000 epoch on the model's own timeline. The engine calls it before the run
and adds the elapsed time to it for each ephemeris time `et` it passes to
[`ephemerides_frame_lpi`](@ref) and [`ephemerides_sun_ii`](@ref).

`SpiceEphemeridesModel` returns SPICE's ephemeris time (TDB);
`SimpleEphemeridesModel` returns leap-second-free UTC seconds past J2000 plus
its `reference_epoch_seconds`.
"""
function ephemerides_epoch_seconds end

"""
    ephemerides_requires_spice(model) -> Bool

Whether the model reads SPICE kernels. When it does, the planet-fixed position
and velocity transforms ask SPICE for the full state transform rather than
[`ephemerides_frame_lpi`](@ref), and the run counts the model's SPICE frame
calls.
"""
function ephemerides_requires_spice end

"""
    ephemerides_cache_key(model) -> Tuple

A key for the model's timeline and frames: the engine may reuse an ephemeris
cache built for one run in a later run whose model has an equal key. Two models
with equal keys must give the same frames and times.
"""
function ephemerides_cache_key end
