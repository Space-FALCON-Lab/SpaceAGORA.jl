struct SpiceEphemeridesModel <: AbstractEphemeridesModel
end

@kwdef struct SimpleEphemeridesModel <: AbstractEphemeridesModel
    reference_epoch_seconds::Float64 = 0.0
    # NaN is a sentinel meaning "planet-true prime meridian": Earth uses GMST
    # (see planet_frame_lpi below), other planets keep a zero angle at the
    # reference epoch. An explicit finite value selects the legacy linear model
    # θ = pm + ω₃·(et − reference_epoch_seconds) exactly as given.
    prime_meridian_at_reference_rad::Float64 = NaN
end

const _J2000_UTC = DateTime(2000, 1, 1, 12, 0, 0)
const _SPICE_POSITION_KM_TO_M = 1.0e3
const _EARTH_HIGH_PREC_BODY_FIXED_FRAME = "ITRF93"
const _EARTH_FALLBACK_BODY_FIXED_FRAME = "IAU_EARTH"

@inline function _spice_position_j2000_m_unlocked(
    target::AbstractString,
    et::Float64,
    observer::AbstractString
)::SVector{3, Float64}
    # NAIF SPK position outputs are in km; convert once at the SPICE boundary.
    return SVector{3, Float64}(spkpos(target, et, "J2000", "none", observer)[1]) * _SPICE_POSITION_KM_TO_M
end

@inline function spice_position_j2000_m(
    target::AbstractString,
    et::Float64,
    observer::AbstractString
)::SVector{3, Float64}
    return lock(tracked_lock(:spice_ephemeris)) do
        _spice_position_j2000_m_unlocked(target, et, observer)
    end
end

@inline ephemerides_requires_spice(::SpiceEphemeridesModel)::Bool = true
@inline ephemerides_requires_spice(::SimpleEphemeridesModel)::Bool = false

@inline function _initial_time_datetime(initial_time)::DateTime
    base = DateTime(
        Int(initial_time.year),
        Int(initial_time.month),
        Int(initial_time.day),
        Int(initial_time.hour),
        Int(initial_time.minute),
        0
    )
    return base + Millisecond(round(Int, 1000 * Float64(initial_time.second)))
end

@inline function ephemerides_time_seconds(initial_time, ::SpiceEphemeridesModel)::Float64
    start_epoch = from_utc(_initial_time_datetime(initial_time))
    return lock(tracked_lock(:spice_body)) do
        utc2et(to_utc(start_epoch))
    end
end

@inline function ephemerides_time_seconds(initial_time, model::SimpleEphemeridesModel)::Float64
    start_time = _initial_time_datetime(initial_time)
    elapsed_ms = Dates.value(start_time - _J2000_UTC)
    return model.reference_epoch_seconds + elapsed_ms / 1000.0
end

@inline function _rotation_about_spin_axis(θ::Float64)::SMatrix{3, 3, Float64}
    c = cos(θ)
    s = sin(θ)
    return @SMatrix [
        c  s  0.0
       -s  c  0.0
        0.0  0.0  1.0
    ]
end

@inline function _spice_body_fixed_frame(planet)::String
    return planet.name == "Moon" ? "MOON_PA_DE421" : (planet.name == "Earth" ? _EARTH_HIGH_PREC_BODY_FIXED_FRAME : "IAU_" * uppercase(planet.name))
end

function _spice_planet_frame_lpi(planet, et::Float64)::SMatrix{3, 3, Float64}
    return lock(tracked_lock(:spice_frame)) do
        if planet.name == "Earth"
            try
                return SMatrix{3, 3, Float64}(pxform("J2000", _EARTH_HIGH_PREC_BODY_FIXED_FRAME, et))
            catch
                return SMatrix{3, 3, Float64}(pxform("J2000", _EARTH_FALLBACK_BODY_FIXED_FRAME, et))
            end
        end
        return SMatrix{3, 3, Float64}(pxform("J2000", _spice_body_fixed_frame(planet), et))
    end
end

@inline function planet_frame_lpi(planet, et::Float64, ::SpiceEphemeridesModel)::SMatrix{3, 3, Float64}
    return _spice_planet_frame_lpi(planet, et)
end

# IAU-82 GMST as a function of the simple model's own timeline: leap-second-free
# UTC seconds past the J2000 epoch (2000-01-01T12:00 UTC), used directly as UT1
# (|UT1 − UTC| < 0.9 s, under 4e-3 deg of Earth rotation).
@inline function _earth_gmst_iau82_rad(ut1_seconds_past_j2000::Float64)::Float64
    Tu = ut1_seconds_past_j2000 / (36525.0 * 86400.0)
    gmst_s = @evalpoly(Tu, 67310.54841, 3.164400184812866e9, 0.093104, -6.2e-6)
    return mod2pi(rem(gmst_s, 86400.0) * (2π / 86400.0))
end

@inline function planet_frame_lpi(planet, et::Float64, model::SimpleEphemeridesModel)::SMatrix{3, 3, Float64}
    elapsed = et - model.reference_epoch_seconds
    pm = model.prime_meridian_at_reference_rad
    θ = if isnan(pm)
        # Planet-true default. The J2000 x-axis is the vernal equinox, not the
        # Greenwich meridian: Earth's rotation angle relative to it is GMST
        # (280.46 deg at the J2000 epoch), so a zero angle here would misplace
        # every geographic (lat/lon-keyed) model evaluation in longitude.
        planet.name == "Earth" ? _earth_gmst_iau82_rad(elapsed) : planet.ω[3] * elapsed
    else
        pm + planet.ω[3] * elapsed
    end
    # With J2000 as the internal inertial frame, the spin-axis rotation is the
    # direct J2000→PCPF transform for the simple ephemeris model.
    return _rotation_about_spin_axis(θ)
end

@inline ephemerides_cache_key(::SpiceEphemeridesModel) = (:spice,)

@inline function ephemerides_cache_key(model::SimpleEphemeridesModel)
    pm = model.prime_meridian_at_reference_rad
    return (
        :simple,
        round(Int64, model.reference_epoch_seconds * 1e6),
        # NaN is the planet-true-default sentinel (round would throw on it);
        # typemin is unreachable from any physically sane explicit override.
        isnan(pm) ? typemin(Int64) : round(Int64, pm * 1e12)
    )
end
