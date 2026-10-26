# Construct MPC parameters from the active SpaceAGORA simulation objects.
@inline _mpc_scenario_args(obj) = hasproperty(obj, :args) ? getproperty(obj, :args) : obj

@inline function _maybe_property(obj, name::Symbol, default=nothing)
    return hasproperty(obj, name) ? getproperty(obj, name) : default
end

function _planet_mu(planet)::Float64
    value = _maybe_property(planet, :μ, _maybe_property(planet, :mu, nothing))
    value === nothing && throw(ArgumentError("Planet does not expose μ/mu."))
    return Float64(value)
end

function _planet_radius_m(planet)::Float64
    value = _maybe_property(planet, :Rp_e, _maybe_property(planet, :Re, nothing))
    value === nothing && throw(ArgumentError("Planet does not expose Rp_e/Re."))
    return Float64(value)
end

# Rotation rate about the pole: a scalar, or the z component of the planet's
# spin vector in its planet-fixed frame.
function _planet_rotation_rate_rad_s(planet)::Float64
    value = _maybe_property(planet, :ω, _maybe_property(planet, :Ω, 0.0))
    return Float64(value isa Number ? value : value[3])
end

"""
    mpc_params_from_spaceagora(args_or_params) -> KSPropagationParams

The planet constants of a `SimulationConfiguration` (or run parameters `p`).
"""
function mpc_params_from_spaceagora(obj)::KSPropagationParams
    args = _mpc_scenario_args(obj)
    planet = args.environment_model.planet
    return KSPropagationParams(
        Re=_planet_radius_m(planet),
        μ=_planet_mu(planet),
        J2=Float64(_maybe_property(planet, :J2, 0.0)),
        Ω=_planet_rotation_rate_rad_s(planet),
    )
end

"""
    spacecraft_mass_kg(spacecraft) -> Float64

Dry plus propellant mass (kg), from `dry_mass`, the root link, or the sum of the links.
"""
function spacecraft_mass_kg(spacecraft)::Float64
    dry = _maybe_property(spacecraft, :dry_mass, nothing)
    prop = _maybe_property(spacecraft, :prop_mass, 0.0)
    if dry !== nothing
        return Float64(dry) + Float64(prop)
    end
    root = _maybe_property(spacecraft, :root, nothing)
    if root !== nothing && hasproperty(root, :m)
        return Float64(root.m) + Float64(prop)
    end
    links = _maybe_property(spacecraft, :links, nothing)
    links === nothing && throw(ArgumentError("Cannot derive spacecraft mass from this object."))
    return sum(Float64(_maybe_property(link, :m, 0.0)) for link in links) + Float64(prop)
end

function spacecraft_reference_areas(spacecraft; controlled_panel_links)
    links = _maybe_property(spacecraft, :links, nothing)
    links === nothing && throw(ArgumentError("Cannot derive reference areas from spacecraft without links."))
    controlled = Set(Int.(controlled_panel_links))
    total = 0.0
    controlled_area = 0.0
    for (idx, link) in pairs(links)
        area = Float64(_maybe_property(link, :ref_area, 0.0))
        total += max(0.0, area)
        idx in controlled && (controlled_area += max(0.0, area))
    end
    return (
        total_reference_area_m2=total,
        bus_reference_area_m2=max(0.0, total - controlled_area),
        controllable_area_m2=controlled_area,
    )
end

"""
    mpc_config_from_spaceagora(args; mode, spacecraft_index, controlled_panel_links,
        constraints, drag_coefficient, qdot_max_w_cm2, heat_load_max_j_cm2,
        drag_max_n, area_slew_max_m2_s, target_energy_mj_kg, area_weight,
        area_slew_weight, slack_weight, target_energy_weight,
        max_depletion_energy_weight, osqp_eps_abs, osqp_eps_rel, osqp_max_iter)

An [`AerobrakingMPCConfig`](@ref) whose areas and mass come from spacecraft
`spacecraft_index` of `args`: the controlled panel links give the controllable
area and the other links the bus area.
"""
function mpc_config_from_spaceagora(
    obj;
    mode::AerobrakingMPCMode,
    spacecraft_index::Integer,
    controlled_panel_links,
    constraints::AerobrakingMPCConstraintSet,
    drag_coefficient::Real,
    qdot_max_w_cm2::Real,
    heat_load_max_j_cm2::Real,
    drag_max_n::Real,
    area_slew_max_m2_s::Real,
    target_energy_mj_kg::Real,
    area_weight::Real,
    area_slew_weight::Real,
    slack_weight::Real,
    target_energy_weight::Real,
    max_depletion_energy_weight::Real,
    osqp_eps_abs::Real,
    osqp_eps_rel::Real,
    osqp_max_iter::Integer,
)
    args = _mpc_scenario_args(obj)
    sc = args.dynamics_model.spacecraft[Int(spacecraft_index)]
    areas = spacecraft_reference_areas(sc; controlled_panel_links=controlled_panel_links)
    return AerobrakingMPCConfig(;
        mode=mode,
        bus_reference_area_m2=areas.bus_reference_area_m2,
        controllable_area_m2=areas.controllable_area_m2,
        mass_kg=spacecraft_mass_kg(sc),
        drag_coefficient=Float64(drag_coefficient),
        qdot_max_w_cm2=Float64(qdot_max_w_cm2),
        heat_load_max_j_cm2=Float64(heat_load_max_j_cm2),
        drag_max_n=Float64(drag_max_n),
        area_slew_max_m2_s=Float64(area_slew_max_m2_s),
        constraint_flags(constraints)...,
        target_energy_mj_kg=Float64(target_energy_mj_kg),
        area_weight=Float64(area_weight),
        area_slew_weight=Float64(area_slew_weight),
        slack_weight=Float64(slack_weight),
        target_energy_weight=Float64(target_energy_weight),
        max_depletion_energy_weight=Float64(max_depletion_energy_weight),
        osqp_eps_abs=Float64(osqp_eps_abs),
        osqp_eps_rel=Float64(osqp_eps_rel),
        osqp_max_iter=Int(osqp_max_iter),
    )
end
