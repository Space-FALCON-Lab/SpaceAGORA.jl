#=
"""
    SpaceAGORA case setup.

    This file reads planet constants, spacecraft mass, exposed area bounds, and
    density calls from SpaceAGORA objects. Limits and tuning are passed in by
    the example file.
"""
=#
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

function _planet_rotation_rate_rad_s(planet)::Float64
    value = _maybe_property(planet, :ω, _maybe_property(planet, :Ω, _maybe_property(planet, :omega, 0.0)))
    value isa Number && return Float64(value)
    try
        return norm(value)
    catch
        return 0.0
    end
end

function mpc_params_from_spaceagora(obj)::AerobrakingMPCParams
    args = _mpc_scenario_args(obj)
    planet = args.environment_model.planet
    return AerobrakingMPCParams(
        Re=_planet_radius_m(planet),
        μ=_planet_mu(planet),
        J2=Float64(_maybe_property(planet, :J2, 0.0)),
        Ω=_planet_rotation_rate_rad_s(planet),
    )
end

function mpc_prediction_gravity_model(obj)::Symbol
    args = _mpc_scenario_args(obj)
    effs = args.dynamics_model.dynamic_effectors
    names = Symbol.(nameof.(typeof.(effs)))
    :GravitationalHarmonicsModel in names && return :harmonics
    :InverseSquaredJ2GravityModel in names && return :inverse_square_j2
    :InverseSquaredGravityModel in names && return :inverse_square
    return :unknown
end

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
    cfg = AerobrakingMPCConfig(
        mode=mode,
        bus_reference_area_m2=areas.bus_reference_area_m2,
        controllable_area_m2=areas.controllable_area_m2,
        mass_kg=spacecraft_mass_kg(sc),
        drag_coefficient=Float64(drag_coefficient),
        qdot_max_w_cm2=Float64(qdot_max_w_cm2),
        heat_load_max_j_cm2=Float64(heat_load_max_j_cm2),
        drag_max_n=Float64(drag_max_n),
        area_slew_max_m2_s=Float64(area_slew_max_m2_s),
        use_constraints=any((constraints.heat_rate, constraints.heat_load, constraints.drag, constraints.slew)),
        use_slew_constraint=constraints.slew,
        use_qdot_constraint=constraints.heat_rate,
        use_heat_load_constraint=constraints.heat_load,
        use_drag_constraint=constraints.drag,
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
    return cfg
end

function _spaceagora_density_getter(density_getter)
    density_getter !== nothing && return density_getter
    if isdefined(Main, :SpaceAGORA)
        if isdefined(Main.SpaceAGORA, :getDensity)
            return Main.SpaceAGORA.getDensity
        end
        if isdefined(Main.SpaceAGORA, :SimulationModel) &&
                isdefined(Main.SpaceAGORA.SimulationModel, :getDensity)
            return Main.SpaceAGORA.SimulationModel.getDensity
        end
    end
    throw(ArgumentError("SpaceAGORA.getDensity is not loaded; pass density_getter explicitly."))
end

function density_and_gradient_from_spaceagora(
    density_model,
    altitude_m::Real,
    latitude::Real,
    longitude::Real,
    elapsed_time_s::Real,
    wind::Bool,
    p;
    dh_m::Real=10.0,
    density_getter=nothing,
)
    getter = _spaceagora_density_getter(density_getter)
    h = Float64(altitude_m)
    rho, temperature, wind_vec = getter(density_model, h, Float64(latitude), Float64(longitude), Float64(elapsed_time_s), wind, p)
    step = max(abs(Float64(dh_m)), 1.0)
    rho_plus = getter(density_model, h + step, Float64(latitude), Float64(longitude), Float64(elapsed_time_s), wind, p)[1]
    rho_minus = getter(density_model, h - step, Float64(latitude), Float64(longitude), Float64(elapsed_time_s), wind, p)[1]
    return (
        rho=Float64(rho),
        drho_dh=(Float64(rho_plus) - Float64(rho_minus)) / (2.0 * step),
        temperature=Float64(temperature),
        wind=wind_vec,
    )
end

function density_function_from_spaceagora(
    obj;
    latitude::Real=0.0,
    longitude::Real=0.0,
    wind::Bool=false,
    density_getter=nothing,
)
    args = _mpc_scenario_args(obj)
    density_context = hasproperty(obj, :args) ? obj : args
    density_model = args.environment_model.density_model
    getter = _spaceagora_density_getter(density_getter)
    return (altitude_m, elapsed_time_s) -> Float64(
        getter(
            density_model,
            Float64(altitude_m),
            Float64(latitude),
            Float64(longitude),
            Float64(elapsed_time_s),
            wind,
            density_context,
        )[1],
    )
end
