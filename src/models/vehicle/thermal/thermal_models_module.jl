module VehicleThermalModels
    using ..AbstractTypes: AbstractThermalModel, AbstractPlanet

    export MaxwellianHeat, SuttonGravesHeat, TabularHeat, getHeatRate

    """
        getHeatRate(model, S, T, ρ, v, α) -> Float64

    The heat rate of one thermal link, which every `AbstractThermalModel`
    implements for its own type. The built-in models return W/cm², the unit of
    the saved `heat_rate` columns, and the engine integrates the rate into each
    link's `heat_loads` state without converting it.

    - `S`: molecular speed ratio, `sqrt(γ / 2)` times the Mach number (no unit).
    - `T`: free-stream temperature, K.
    - `ρ`: free-stream density, kg/m³.
    - `v`: speed relative to the atmosphere (planet rotation and wind removed),
      m/s.
    - `α`: incidence angle of the link to the flow, rad.

    Whenever the engine samples heating where the density is positive, it calls
    this for each link of each spacecraft, with all arguments `Float64`, and
    stores a non-finite or negative result as zero. Before integration starts,
    `run_simulation` checks that the configured thermal model has this method
    and stops with an `ArgumentError` otherwise.
    """
    function getHeatRate end

    include(joinpath(@__DIR__, "thermal_models.jl"))
end
