module DynamicsTranslational

using StaticArrays

include(joinpath(@__DIR__, "position_kinematics.jl"))
include(joinpath(@__DIR__, "point_mass_dynamics.jl"))

export position_derivative
export zero_position_derivative
export acceleration_from_force
export mass_derivative
export assign_full_translational_rhs!
export assign_slow_translational_rhs!
export assign_control_only_translational_rhs!
export assign_force_only_translational_rhs!

end # module DynamicsTranslational
