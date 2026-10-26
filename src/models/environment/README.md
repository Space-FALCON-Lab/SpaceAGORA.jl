# Environmental models and field kernels

Environmental models describe the surroundings sampled by a simulation:
planet properties, atmosphere, ephemerides, radiation inputs and physical
fields. They supply environmental values; spacecraft-specific force and torque
adapters live in `../dynamics/effectors/`.

The distinction is about responsibility, not a new Julia namespace boundary.
The small `aerodynamics/aerodynamic_effectors.jl` and
`gravity/gravity_effectors.jl` modules are compatibility facades that forward
existing dynamic-effector names. They do not contain independent force-model
implementations.

Analytic gravity kernels such as central acceleration, Encke differencing and
J2 secular rates live in `gravity/`. Runtime state access, model types and
spacecraft wrench adapters live in `../dynamics/effectors/gravity/`. The ordered
`gravity/gravity_models.jl` include manifest preserves the original loading
sequence and evaluates both sets of definitions in the existing
`DynamicEffectors.GravityEffectors` module. That deliberate composition keeps
method/type identities and numerical expressions intact.
