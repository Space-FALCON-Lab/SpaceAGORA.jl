# Coupled force and torque models

This directory adapts physical models to a spacecraft state and the sampled
environment. Its methods supply forces and torques to the dynamics engine.
Atmosphere, ephemeris and field sampling remain environmental responsibilities;
applying those samples to spacecraft mass, attitude, geometry and hardware is a
dynamic-effector responsibility.

`force_torque_models.jl` composes the existing `DynamicEffectors` namespace;
`force_torque_models/gravitational.jl` marks which of its effectors model
gravity (`is_gravitational`, read by the accelerometer).
`perturbations.jl` is an ordered include manifest for gravity, radiation,
magnetics and continuous attitude-control effectors. The separate files are
evaluated in the existing `PerturbationEffectors` module. In particular, the
LVLH cascade remains a continuous body-torque law evaluated by the RHS. It is
not a sampled GNC callback.

The analytic gravity adapters in `gravity/` are still defined in
`DynamicEffectors.GravityEffectors`. Their mathematical kernels live in
`../../environment/gravity/`; the compatibility composition file there keeps
the original definition order. Splitting these files changes neither model
identities nor the single compiled force kernels shared by scalar and batch
routes.

Thread eligibility and multibody worker selection belong to
`orchestration/parallel`. Historical aerodynamic helper bindings alias that
shared policy implementation so existing callers retain the same functions
and per-solve cache.
