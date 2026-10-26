module EphemeridesModels

using Dates
using StaticArrays
using LinearAlgebra
using SPICE

using ...RuntimeServices: SPICE_LOCK, tracked_lock
using ..AbstractTypes: AbstractEphemeridesModel

export SpiceEphemeridesModel, SimpleEphemeridesModel
export ephemerides_frame_lpi, ephemerides_sun_ii, ephemerides_epoch_seconds
export ephemerides_requires_spice, ephemerides_cache_key

include(joinpath(@__DIR__, "ephemerides_hooks.jl"))
include(joinpath(@__DIR__, "simple_ephemerides.jl"))

end # module EphemeridesModels
