#=
"""
    MPC constraint selection.

    The active constraints must be stated by the user or example file. This
    avoids a default case where a limit is active only because it was hidden in
    the source code.
"""
=#
Base.@kwdef struct AerobrakingMPCConstraintSet
    heat_rate::Bool
    heat_load::Bool
    drag::Bool
    slew::Bool
end

const _MPC_CONSTRAINT_SYMBOLS = Set((:heat_rate, :heat_load, :drag, :slew))

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

function constraint_active(set::AerobrakingMPCConstraintSet, name::Symbol)::Bool
    name == :heat_rate && return set.heat_rate
    name == :heat_load && return set.heat_load
    name == :drag && return set.drag
    name == :slew && return set.slew
    throw(ArgumentError("Unsupported MPC constraint $(repr(name))."))
end

function constraint_names(set::AerobrakingMPCConstraintSet)
    names = Symbol[]
    set.heat_rate && push!(names, :heat_rate)
    set.heat_load && push!(names, :heat_load)
    set.drag && push!(names, :drag)
    set.slew && push!(names, :slew)
    return names
end

function apply_constraints(config::AerobrakingMPCConfig, set::AerobrakingMPCConstraintSet)
    return AerobrakingMPCConfig(
        mode=config.mode,
        bus_reference_area_m2=config.bus_reference_area_m2,
        controllable_area_m2=config.controllable_area_m2,
        mass_kg=config.mass_kg,
        drag_coefficient=config.drag_coefficient,
        qdot_max_w_cm2=config.qdot_max_w_cm2,
        heat_load_max_j_cm2=config.heat_load_max_j_cm2,
        drag_max_n=config.drag_max_n,
        area_slew_max_m2_s=config.area_slew_max_m2_s,
        use_constraints=any((set.heat_rate, set.heat_load, set.drag, set.slew)),
        use_slew_constraint=set.slew,
        use_qdot_constraint=set.heat_rate,
        use_heat_load_constraint=set.heat_load,
        use_drag_constraint=set.drag,
        target_energy_mj_kg=config.target_energy_mj_kg,
        area_weight=config.area_weight,
        area_slew_weight=config.area_slew_weight,
        slack_weight=config.slack_weight,
        target_energy_weight=config.target_energy_weight,
        max_depletion_energy_weight=config.max_depletion_energy_weight,
        osqp_eps_abs=config.osqp_eps_abs,
        osqp_eps_rel=config.osqp_eps_rel,
        osqp_max_iter=config.osqp_max_iter,
    )
end
