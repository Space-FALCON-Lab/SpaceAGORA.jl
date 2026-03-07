struct SpiceEphemeridesModel <: AbstractEphemeridesModel
end

@kwdef struct SimpleEphemeridesModel <: AbstractEphemeridesModel
    reference_epoch_seconds::Float64 = 0.0
    prime_meridian_at_reference_rad::Float64 = 0.0
end

const _J2000_UTC = DateTime(2000, 1, 1, 12, 0, 0)

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
    return lock(SPICE_LOCK) do
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

@inline function planet_frame_lpi(planet, et::Float64, ::SpiceEphemeridesModel)::SMatrix{3, 3, Float64}
    return lock(SPICE_LOCK) do
        SMatrix{3, 3, Float64}(pxform("J2000", "IAU_$(planet.name)", et)) * planet.J2000_to_pci'
    end
end

@inline function planet_frame_lpi(planet, et::Float64, model::SimpleEphemeridesModel)::SMatrix{3, 3, Float64}
    θ = model.prime_meridian_at_reference_rad + planet.ω[3] * (et - model.reference_epoch_seconds)
    return _rotation_about_spin_axis(θ)
end

@inline ephemerides_cache_key(::SpiceEphemeridesModel) = (:spice,)

@inline function ephemerides_cache_key(model::SimpleEphemeridesModel)
    return (
        :simple,
        round(Int64, model.reference_epoch_seconds * 1e6),
        round(Int64, model.prime_meridian_at_reference_rad * 1e12)
    )
end
