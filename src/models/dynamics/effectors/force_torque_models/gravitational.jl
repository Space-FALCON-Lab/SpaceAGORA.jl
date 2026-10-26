# Which dynamic effectors model gravity. Included at the end of
# `DynamicEffectors` (force_torque_models.jl), where the gravity models are in
# scope.

"""
    is_gravitational(effector) -> Bool

Whether a dynamic effector models gravity: the central body's field and its
harmonics, third bodies, and the gravity-gradient torque. An accelerometer
does not feel these effectors' forces (see `Accelerometer`). A user-defined
gravity model should add a method returning true.
"""
is_gravitational(::Any) = false
is_gravitational(::Union{ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel,
                         NBodyGravityModel, GravitationalHarmonicsModel, GravityGradientTorqueModel}) = true
