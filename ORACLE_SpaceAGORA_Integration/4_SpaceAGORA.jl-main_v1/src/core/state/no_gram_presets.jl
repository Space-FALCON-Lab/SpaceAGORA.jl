module NoGramPresets

using ..AbstractTypes: AbstractPlanet, AbstractDensityModel
using ..Planets: Earth, Mars, Venus
using ..EnvironmentModels: NoAtmosphereModel, ExponentialAtmosphereModel
using ..EphemeridesModels: SimpleEphemeridesModel
using ..VehicleThermalModels: MaxwellianHeat
using ..SimConfig: EnvironmentModel

export make_no_gram_planet, make_no_gram_density_model, make_no_gram_environment

"""
    make_no_gram_planet(planet)

Return a baseline planet object for the no-GRAM onboarding mode without loading
SPICE kernels. Supported keys are `:earth`, `:mars`, and `:venus`.
"""
@inline function make_no_gram_planet(planet::AbstractPlanet)
    return planet
end

function make_no_gram_planet(planet::Symbol)::AbstractPlanet
    key = Symbol(lowercase(String(planet)))
    if key === :earth
        return Earth()
    elseif key === :mars
        return Mars()
    elseif key === :venus
        return Venus()
    end
    throw(ArgumentError("Unsupported no-GRAM planet $(repr(planet)). Supported planets: :earth, :mars, :venus."))
end

@inline make_no_gram_planet(planet::AbstractString) = make_no_gram_planet(Symbol(strip(lowercase(planet))))

"""
    make_no_gram_density_model(planet, atmosphere)

Construct the documented first-class no-GRAM atmosphere model for a baseline
run. The minimal quickstart baseline uses `atmosphere=:none`, which yields
`NoAtmosphereModel()`. `atmosphere=:none` yields `NoAtmosphereModel()` and
`atmosphere=:exponential` yields an `ExponentialAtmosphereModel(planet)`. Pass a
prebuilt `AbstractDensityModel` instance directly to use a custom analytic
model such as `PiecewiseExponentialAtmosphereModel(...)`.
"""
@inline function make_no_gram_density_model(planet::AbstractPlanet, density_model::AbstractDensityModel)
    return density_model
end

function make_no_gram_density_model(planet::AbstractPlanet, atmosphere::Symbol)::AbstractDensityModel
    key = Symbol(lowercase(String(atmosphere)))
    if key === :none
        return NoAtmosphereModel()
    elseif key === :exponential
        return ExponentialAtmosphereModel(planet)
    end
    throw(ArgumentError("Unsupported no-GRAM atmosphere $(repr(atmosphere)). Supported options: :none, :exponential."))
end

@inline make_no_gram_density_model(planet::AbstractPlanet, atmosphere::AbstractString) =
    make_no_gram_density_model(planet, Symbol(strip(lowercase(atmosphere))))

"""
    make_no_gram_environment(; planet=:earth, atmosphere=:none, EI_km=120.0, wind=false, topography=false)

Build a first-class no-GRAM environment configuration that does not depend on
GRAM assets or SPICE kernels. The default is the minimal quickstart baseline:
Earth, `NoAtmosphereModel()`, and `SimpleEphemeridesModel()`. Use explicit
arguments to opt into other documented fallback atmospheres such as
`ExponentialAtmosphereModel(planet)`.
"""
function make_no_gram_environment(;
    planet::Union{AbstractPlanet, Symbol, AbstractString}=:earth,
    atmosphere::Union{AbstractDensityModel, Symbol, AbstractString}=:none,
    EI_km::Real=120.0,
    wind::Bool=false,
    topography::Bool=false,
    topo_degree::Integer=0,
    topo_order::Integer=0,
    thermal_model=nothing
)
    planet_model = make_no_gram_planet(planet)
    density_model = make_no_gram_density_model(planet_model, atmosphere)
    thermal_model_resolved = isnothing(thermal_model) ?
        MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet_model) :
        thermal_model

    return EnvironmentModel(
        planet=planet_model,
        EI=Float64(EI_km),
        density_model=density_model,
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=thermal_model_resolved,
        topography=topography,
        topo_degree=Int(topo_degree),
        topo_order=Int(topo_order),
        wind=wind
    )
end

end # module NoGramPresets
