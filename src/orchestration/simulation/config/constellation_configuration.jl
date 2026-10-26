# The existing spacecraft collection and selected dynamics effectors for one run.
# Included inside SpacecraftModels to preserve DynamicsModel type identity and public access.
# The spacecraft stay a plain vector; build_constellation (src/models/vehicle/constellations/)
# fills one from orbital planes.

import ...GNCRuntime
using ...GNCRuntime: SpacecraftScope, AllSpacecraft, BoundSpacecraft, CoupledSpacecraft, spacecraft_scope

"""
    ScopedEffector

A dynamic effector with the scope `scoped` gave it: the effector itself in
`effector`, and the spacecraft it acts on in `scope`. Build one with `scoped`.
"""
struct ScopedEffector{E, S <: SpacecraftScope}
    effector::E
    scope::S
end

"""
    scoped(effector, scope) -> ScopedEffector

A dynamic effector, as an entry of a `DynamicsModel`'s effector tuple, that acts
only on the spacecraft `scope` names: `BoundSpacecraft(positions)`, a
`SpacecraftSelection` (by position, id, constellation plane or slot, or a
predicate), or `AllSpacecraft()`, which is the same as no scope.

```julia
DynamicsModel(spacecraft, (
    InverseSquaredJ2GravityModel(),
    scoped(AerodynamicCoefficientConstant(), SpacecraftSelection(planes=1)),
))
```

Outside its scope the effector does not exist for that spacecraft: it is not
evaluated, samples no environment, adds no force or torque and writes no
diagnostics there, on every solver route. The run resolves the scope to run
positions once, when it starts, and its configuration then holds the effector
itself with the resolved positions (see `DynamicsModel`). A `CoupledSpacecraft`
scope is refused: a force model acts on one spacecraft at a time, and forces
between spacecraft are interlinks.
"""
function scoped(effector, scope)
    scope isa SpacecraftScope || throw(ArgumentError(
        "scoped takes a scope AllSpacecraft(), BoundSpacecraft(positions) or SpacecraftSelection(...), " *
        "not $(repr(scope))."))
    scope isa CoupledSpacecraft && throw(ArgumentError(
        "scoped does not take a CoupledSpacecraft: a force model acts on one spacecraft at a time, and forces " *
        "between spacecraft are interlinks (link_force). Give the spacecraft it acts on as BoundSpacecraft or " *
        "SpacecraftSelection."))
    effector isa ScopedEffector && throw(ArgumentError(
        "The effector is already scoped; give scoped the effector itself."))
    spacecraft_scope(effector) isa AllSpacecraft || throw(ArgumentError(
        "$(nameof(typeof(effector))) already declares its scope through spacecraft_scope; the run applies that " *
        "scope, so give the effector without scoped."))
    return ScopedEffector(effector, scope)
end

GNCRuntime.spacecraft_scope(e::ScopedEffector) = e.scope

"""
The spacecraft each dynamic effector of a run acts on, resolved to run positions
when the run starts (internal). Entry `j` describes effector `j` of the run's
tuple: `nothing` when it acts on every spacecraft, otherwise the sorted run
positions and a mask over the run's `n_spacecraft`. `atmosphere` marks the
spacecraft that some atmosphere-consuming effector acts on (`nothing` when one
acts on every spacecraft). A run with no scoped effector has no table.
"""
struct EffectorScopes
    positions::Vector{Union{Nothing, Vector{Int}}}
    masks::Vector{Union{Nothing, BitVector}}
    declared::Vector{Any}
    n_spacecraft::Int
    atmosphere::Union{Nothing, BitVector}
end

# Whether effector `eff_idx` of the run acts on spacecraft `sat_idx`.
@inline _effector_applies(::Nothing, eff_idx::Int, sat_idx::Int)::Bool = true
@inline function _effector_applies(scopes::EffectorScopes, eff_idx::Int, sat_idx::Int)::Bool
    @inbounds mask = scopes.masks[eff_idx]
    return mask === nothing || @inbounds mask[sat_idx]
end

"""
    DynamicsModel(spacecraft, dynamic_effectors)

Struct containing the information needed for the dynamics propagation of all spacecraft.

- `spacecraft`: vector of spacecraft models, one per propagated vehicle, each represented by its root link (main bus or core body). This is the vector over which the threads are parallelized, so each thread handles the dynamics propagation of one spacecraft and its attached links.
- `dynamic_effectors`: tuple of force and torque models (gravity, drag, etc.) applied during the dynamics propagation. An entry `scoped(effector, scope)` acts only on the spacecraft its scope names.

After `res = run_simulation(config; return_results=true)`, `res.configuration.dynamics_model.spacecraft[i]` is the spacecraft that ran.
In that configuration a scoped entry has been resolved: its `dynamic_effectors` hold the effector itself, and
`effector_scopes` the run positions each one acts on (`nothing` when no effector is scoped). Rebuilding a
`DynamicsModel` from those bare effectors drops the scopes; keep the tuple with the `scoped` entries to reuse them.
"""
struct DynamicsModel{T_Effectors<:Tuple}
    spacecraft::Vector{SpacecraftModel} # Vector of root links (main bus or core bodies)
    dynamic_effectors::T_Effectors # Tuple of dynamic effector models (gravity, drag, etc.)
    effector_scopes::Union{Nothing, EffectorScopes} # Resolved effector scopes of a run (internal), or nothing

    function DynamicsModel(roots::Vector{SpacecraftModel}, dynamic_effectors::T_Effectors,
                           effector_scopes::Union{Nothing, EffectorScopes}=nothing) where {T_Effectors<:Tuple}
        new{T_Effectors}(roots, dynamic_effectors, effector_scopes)
    end
end
