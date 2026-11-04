# Compose selected environment models without evaluating their physics.
# Included inside SimConfig by configuration.jl.

# 2.7. Environment Model
# TODO: Convert all the strings to abstract types to avoid needing if-else statements in complete passage and other functions. This will also make it easier to add new models in the future without needing to change the main code.
# i) Struct
@kwdef struct EnvironmentModel{P <: AbstractPlanet, D <: AbstractDensityModel, E <: AbstractEphemeridesModel, T <: AbstractThermalModel}
    # Physical environment model
    planet::P # Planet for which to run the simulation (used for gravity model, atmospheric model, etc.)
    EI::Float64 # Entry Interface altitude in km (used for determining when to start applying atmospheric effects)
    density_model::D # Atmospheric model to use (Constant, Exponential, GRAM, NRLMSISE-00)
    ephemerides_model::E = SpiceEphemeridesModel() # Planet-frame/ephemeris backend (SPICE for high fidelity, simplified analytic mode for onboarding)
    topography::Bool = false # Whether to include topography in the simulation for altitude calculation
    topo_degree::Int = 90 # Maximum degree of spherical harmonics for topography
    topo_order::Int = 90 # Maximum order of spherical harmonics for topography
    wind::Bool = true # Whether to include wind in the simulation for atmospheric effects
    thermal_model::T # Thermal model to use (Maxwellian heat transfer, Convective and Radiative)

    # ii) Inner constructor with data type conversion and validation
    function EnvironmentModel(
        planet::P,
        EI::Real,
        density_model::D,
        ephemerides_model::E,
        topography::Bool,
        topo_degree::Integer,
        topo_order::Integer,
        wind::Bool,
        thermal_model::T
    ) where {P <: AbstractPlanet, D <: AbstractDensityModel, E <: AbstractEphemeridesModel, T <: AbstractThermalModel}
        EI >= 0 || throw(ArgumentError("EnvironmentModel.EI must be >= 0 km; got $EI."))
        topo_degree >= 0 || throw(ArgumentError("EnvironmentModel.topo_degree must be >= 0; got $topo_degree."))
        topo_order >= 0 || throw(ArgumentError("EnvironmentModel.topo_order must be >= 0; got $topo_order."))
        return new{P, D, E, T}(
            planet,
            Float64(EI),
            density_model,
            ephemerides_model,
            topography,
            Int(topo_degree),
            Int(topo_order),
            wind,
            thermal_model
        )
    end
end # struct EnvironmentModel
