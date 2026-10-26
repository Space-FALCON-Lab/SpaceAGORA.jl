# Selection of the active path (heat rate, heat load, drag) and actuator (slew) limits.
"""
Active MPC limits: `heat_rate`, `heat_load`, `drag` and actuator `slew`.
Build one with [`mpc_constraints`](@ref).
"""
Base.@kwdef struct AerobrakingMPCConstraintSet
    heat_rate::Bool
    heat_load::Bool
    drag::Bool
    slew::Bool
end

const _MPC_CONSTRAINT_SYMBOLS = Set((:heat_rate, :heat_load, :drag, :slew))

"""
    mpc_constraints(; heat_rate, heat_load, drag, slew)
    mpc_constraints(names...)

Select the active MPC limits, by keyword or by name (`:heat_rate`,
`:heat_load`, `:drag`, `:slew`); returns an [`AerobrakingMPCConstraintSet`](@ref).
"""
function mpc_constraints(;
    heat_rate::Bool,
    heat_load::Bool,
    drag::Bool,
    slew::Bool,
)
    return AerobrakingMPCConstraintSet(
        heat_rate=heat_rate,
        heat_load=heat_load,
        drag=drag,
        slew=slew,
    )
end

function mpc_constraints(active::Union{Symbol, AbstractString}...)
    requested = Set(Symbol.(active))
    for item in requested
        item in _MPC_CONSTRAINT_SYMBOLS ||
            throw(ArgumentError("Unsupported MPC constraint $(repr(item)). Valid constraints are $(sort!(collect(_MPC_CONSTRAINT_SYMBOLS)))."))
    end
    return AerobrakingMPCConstraintSet(
        heat_rate=:heat_rate in requested,
        heat_load=:heat_load in requested,
        drag=:drag in requested,
        slew=:slew in requested,
    )
end

"""
    constraint_active(set, name::Symbol) -> Bool

Whether the limit `name` (`:heat_rate`, `:heat_load`, `:drag` or `:slew`) is active in `set`.
"""
function constraint_active(set::AerobrakingMPCConstraintSet, name::Symbol)::Bool
    name == :heat_rate && return set.heat_rate
    name == :heat_load && return set.heat_load
    name == :drag && return set.drag
    name == :slew && return set.slew
    throw(ArgumentError("Unsupported MPC constraint $(repr(name))."))
end

"""
    constraint_names(set) -> Vector{Symbol}

The active limits of an [`AerobrakingMPCConstraintSet`](@ref), in a fixed order.
"""
function constraint_names(set::AerobrakingMPCConstraintSet)
    names = Symbol[]
    set.heat_rate && push!(names, :heat_rate)
    set.heat_load && push!(names, :heat_load)
    set.drag && push!(names, :drag)
    set.slew && push!(names, :slew)
    return names
end

"""
    apply_constraints(config, set) -> AerobrakingMPCConfig

A copy of `config` with the active limits of `set`.
"""
function apply_constraints(config::AerobrakingMPCConfig, set::AerobrakingMPCConstraintSet)
    fields = NamedTuple{fieldnames(AerobrakingMPCConfig)}(
        ntuple(i -> getfield(config, i), fieldcount(AerobrakingMPCConfig)))
    return AerobrakingMPCConfig(; fields..., constraint_flags(set)...)
end

# The `AerobrakingMPCConfig` switches of a constraint set.
constraint_flags(set::AerobrakingMPCConstraintSet) = (
    use_constraints=any((set.heat_rate, set.heat_load, set.drag, set.slew)),
    use_slew_constraint=set.slew,
    use_qdot_constraint=set.heat_rate,
    use_heat_load_constraint=set.heat_load,
    use_drag_constraint=set.drag,
)
