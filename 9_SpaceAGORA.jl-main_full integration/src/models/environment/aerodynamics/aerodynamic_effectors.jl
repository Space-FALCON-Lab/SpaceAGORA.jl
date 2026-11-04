# Compatibility facade only. Atmospheric sampling belongs to environment;
# the spacecraft aerodynamic wrench models are implemented under dynamics.
module AerodynamicEffectors

using ..DynamicEffectors.AerodynamicEffectors: aerodynamic_coefficient_fM

export aerodynamic_coefficient_fM

end # module AerodynamicEffectors
