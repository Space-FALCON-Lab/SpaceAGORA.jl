module EphemeridesModels

using Dates
using AstroTime
using StaticArrays
using LinearAlgebra
using SPICE

using ...RuntimeServices: SPICE_LOCK, tracked_lock
using ..AbstractTypes: AbstractEphemeridesModel

export SpiceEphemeridesModel, SimpleEphemeridesModel
export ephemerides_time_seconds, planet_frame_lpi, ephemerides_requires_spice, ephemerides_cache_key
export ephemerides_sun_direction_ii, ephemerides_body_direction_ii

include(joinpath(@__DIR__, "simple_ephemerides.jl"))

end # module EphemeridesModels
