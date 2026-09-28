# Central body description and the sampled J2000 -> body-fixed rotation table.
#
# The rotation comes from the same `planet_frame_lpi` the runtime uses for
# latitude/longitude, for both the SPICE and simple ephemerides paths, so a
# ground track drawn by the viewer lands where the CSV says it does. Sampling
# it here means the viewer needs no SPICE kernels; it interpolates between
# samples and extrapolates past the table with `spin_rad_s`.

const ROTATION_MAX_SAMPLES = 4096

"""
    rotation_sample_times(mission_time_s, data_rate_s; max_samples=4096) -> Vector{Float64}

Uniform sample times from 0 to `mission_time_s` inclusive, at `data_rate_s` or
coarser so the table never exceeds `max_samples` entries.
"""
function rotation_sample_times(mission_time_s::Real, data_rate_s::Real; max_samples::Integer=ROTATION_MAX_SAMPLES)::Vector{Float64}
    mission_time = Float64(mission_time_s)
    max_samples >= 2 || throw(ArgumentError("max_samples must be at least 2, got $(max_samples)."))
    mission_time > 0.0 || return [0.0]
    step = max(Float64(data_rate_s), mission_time / (max_samples - 1))
    step > 0.0 || throw(ArgumentError("data_rate_s must be positive, got $(data_rate_s)."))
    times = collect(0.0:step:mission_time)
    if isempty(times) || times[end] < mission_time
        push!(times, mission_time)
    end
    return times
end

@inline function _continuous_quaternion(q::SVector{4, Float64}, previous::Union{Nothing, SVector{4, Float64}})::SVector{4, Float64}
    # q and -q are the same rotation; keep the sign continuous so a viewer
    # slerp between neighbours takes the short arc.
    previous === nothing && return q
    return dot(q, previous) < 0.0 ? -q : q
end

"""
    planet_rotation_table(planet, ephemerides_model, et_start_s, t_s) -> Vector{SVector{4, Float64}}

J2000 -> body-fixed quaternion (scalar-last) at each elapsed time in `t_s`,
sign-continuous along the table.
"""
function planet_rotation_table(
    planet::AbstractPlanet,
    ephemerides_model::AbstractEphemeridesModel,
    et_start_s::Real,
    t_s::AbstractVector{<:Real}
)::Vector{SVector{4, Float64}}
    table = Vector{SVector{4, Float64}}(undef, length(t_s))
    previous = nothing
    for (k, t) in enumerate(t_s)
        dcm = planet_frame_lpi(planet, Float64(et_start_s) + Float64(t), ephemerides_model)
        q = _continuous_quaternion(_svec4(dcm_to_quaternion(SMatrix{3, 3, Float64}(dcm))), previous)
        table[k] = q
        previous = q
    end
    return table
end

@inline function _planet_inertial_frame(planet::AbstractPlanet)::String
    return hasproperty(planet, :inertial_frame) ? String(string(getproperty(planet, :inertial_frame))) : "J2000"
end

"""
    planet_spec(planet, ephemerides_model, et_start_s, t_s) -> PlanetSpec

Planet shape, spin, texture key (`lowercase(planet.name)`) and the rotation
table sampled at `t_s`.
"""
function planet_spec(
    planet::AbstractPlanet,
    ephemerides_model::AbstractEphemeridesModel,
    et_start_s::Real,
    t_s::AbstractVector{<:Real}
)::PlanetSpec
    times = Float64[Float64(t) for t in t_s]
    return PlanetSpec(
        String(planet.name),
        Float64(planet.Rp_e),
        Float64(planet.Rp_p),
        _svec3(planet.ω),
        _planet_inertial_frame(planet),
        lowercase(String(planet.name)),
        times,
        planet_rotation_table(planet, ephemerides_model, et_start_s, times)
    )
end
