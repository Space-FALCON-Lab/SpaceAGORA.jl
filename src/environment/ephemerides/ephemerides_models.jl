module EphemeridesModels

using Dates
using AstroTime
using StaticArrays
using SPICE

using ...RuntimeServices: SPICE_LOCK
using ..AbstractTypes: AbstractEphemeridesModel

export SpiceEphemeridesModel, SimpleEphemeridesModel
export ephemerides_time_seconds, planet_frame_lpi, ephemerides_requires_spice, ephemerides_cache_key

include(joinpath(@__DIR__, "simple_ephemerides.jl"))

end # module EphemeridesModels
