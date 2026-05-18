struct PanelComponent
    link
    side_sign::Float64
end

struct SphereConeReferenceBody
    ref_area::Float64
    nose_radius::Float64
    base_radius::Float64
    cone_half_angle_deg::Float64
    hypersonic_pressure_model::Symbol
end

struct FlatPlateReferenceBody
    ref_area::Float64
    skin_friction_coefficient::Float64
    zero_lift_drag_coefficient::Float64
end

struct DifferentialPanelCommand
    left_area_scale::Float64
    right_area_scale::Float64
    left_cant_deg::Float64
    right_cant_deg::Float64
    left_alpha_deg::Float64
    right_alpha_deg::Float64
end

struct AeroForceState
    label::String
    body_link
    panel_components::Vector{PanelComponent}
    trim_panel_components::Vector{PanelComponent}
    body_alpha_rad::Float64
    panel_alpha_rad::Float64
    trim_panel_alpha_rad::Float64
    panel_area_total_m2::Float64
    panel_area_each_m2::Float64
    trim_panel_area_total_m2::Float64
    trim_panel_area_each_m2::Float64
end

struct CalibratedAeroCase
    body_label::String
    deployed_label::String
    target_beta_high::Float64
    target_beta_low::Float64
    target_beta_ratio::Float64
    achieved_beta_high::Float64
    achieved_beta_low::Float64
    achieved_beta_ratio::Float64
    mass_kg::Float64
    reference_temperature_k::Float64
    reference_mach::Float64
    body_state::AeroForceState
    deployed_state::AeroForceState
    characteristic_length_m::Float64
    body_cda_ref_m2::Float64
    body_cla_ref_m2::Float64
    deployed_cda_ref_m2::Float64
    deployed_cla_ref_m2::Float64
    trim_only_cda_ref_m2::Float64
    trim_only_cla_ref_m2::Float64
end

@inline function _link_ref_area(link)
    return Float64(link.ref_area)
end

function build_body_link(config::PrelimConfig)
    geom = config.geometry
    area = π * geom.base_radius_m^2
    return SphereConeReferenceBody(
        area,
        geom.nose_radius_m,
        geom.base_radius_m,
        geom.cone_half_angle_deg,
        geom.hypersonic_pressure_model,
    )
end

function build_panel_links(panel_area_total_m2::Float64, config::PrelimConfig)
    geom = config.geometry
    panel_area_total_m2 <= 0.0 && return PanelComponent[]
    panel_area_each_m2 = panel_area_total_m2 / geom.panel_count
    width = sqrt(panel_area_each_m2 * geom.panel_aspect_ratio)
    height = panel_area_each_m2 / width
    components = PanelComponent[]
    side_signs = geom.panel_count == 2 ? [-1.0, 1.0] : collect(range(-1.0, 1.0; length=geom.panel_count))
    for side_sign in side_signs
        link = FlatPlateReferenceBody(
            panel_area_each_m2,
            geom.panel_skin_friction_coefficient,
            geom.panel_zero_lift_drag_coefficient,
        )
        push!(components, PanelComponent(link, Float64(side_sign)))
    end
    return components
end

function build_trim_panel_links(panel_area_total_m2::Float64, config::PrelimConfig)
    return build_panel_links(panel_area_total_m2, config)
end

@inline function _flat_plate_cl_cd(alpha_rad::Float64, plate::FlatPlateReferenceBody)
    sinα = sin(alpha_rad)
    cosα = cos(alpha_rad)
    Cf = plate.skin_friction_coefficient
    CD0 = plate.zero_lift_drag_coefficient
    CL = 2.0 * abs(sinα) * sinα * cosα - Cf * sinα
    CD = 2.0 * abs(sinα)^3 + Cf * cosα + CD0
    return CL, CD
end

@inline function _sphere_cone_cl_cd(alpha_rad::Float64, body::SphereConeReferenceBody)
    k = body.nose_radius / body.base_radius
    ratio_sq = k^2
    δ = deg2rad(body.cone_half_angle_deg)
    frustum_factor = 1.0 - ratio_sq * cos(δ)^2
    CN = frustum_factor * cos(δ)^2 * sin(2.0 * alpha_rad)
    CA = (1.0 - sin(δ)^4) * ratio_sq +
         (2.0 * sin(δ)^2 * cos(alpha_rad)^2 + cos(δ)^2 * sin(alpha_rad)^2) * frustum_factor
    CL = CN * cos(alpha_rad) - CA * sin(alpha_rad)
    CD = CA * cos(alpha_rad) + CN * sin(alpha_rad)
    return CL, CD
end

function _component_cl_cd(link::FlatPlateReferenceBody, alpha_rad::Float64, mach::Float64, gamma::Float64)
    return _flat_plate_cl_cd(alpha_rad, link)
end

function _component_cl_cd(link::SphereConeReferenceBody, alpha_rad::Float64, mach::Float64, gamma::Float64)
    return _sphere_cone_cl_cd(alpha_rad, link)
end

function _component_cla_cda(link, alpha_rad::Float64, mach::Float64, gamma::Float64)
    CL, CD = _component_cl_cd(link, alpha_rad, mach, gamma)
    area = _link_ref_area(link)
    return CL * area, CD * area
end

function _sum_component_cla_cda(links::AbstractVector, alpha_rad::Float64, mach::Float64, gamma::Float64)
    cla_total_m2 = 0.0
    cda_total_m2 = 0.0
    for link in links
        primitive_link = hasproperty(link, :link) ? getproperty(link, :link) : link
        cla_m2, cda_m2 = _component_cla_cda(primitive_link, alpha_rad, mach, gamma)
        cla_total_m2 += cla_m2
        cda_total_m2 += cda_m2
    end
    return cla_total_m2, cda_total_m2
end

function _component_cl_cd_area(link, alpha_rad::Float64, mach::Float64, gamma::Float64)
    CL, CD = _component_cl_cd(link, alpha_rad, mach, gamma)
    area = _link_ref_area(link)
    return CL, CD, area
end

function neutral_panel_command(config::PrelimConfig)
    alpha_deg = config.geometry.panel_alpha_deg
    return DifferentialPanelCommand(1.0, 1.0, 0.0, 0.0, alpha_deg, alpha_deg)
end

function make_differential_panel_command(
    config::PrelimConfig;
    favored_side_sign::Real=1.0,
    differential_fraction::Real,
    cant_deg::Real,
)
    frac = clamp(Float64(differential_fraction), 0.0, 1.0)
    cant = clamp(Float64(cant_deg), 0.0, 90.0)
    alpha_neutral = config.geometry.panel_alpha_deg
    alpha_deflected = max(alpha_neutral - frac * config.geometry.panel_lateral_deflection_deg, 0.0)
    if Float64(favored_side_sign) >= 0.0
        return DifferentialPanelCommand(
            max(1.0 - frac, 0.0),
            1.0,
            0.0,
            cant,
            alpha_neutral,
            alpha_deflected,
        )
    end
    return DifferentialPanelCommand(
        1.0,
        max(1.0 - frac, 0.0),
        -cant,
        0.0,
        alpha_deflected,
        alpha_neutral,
    )
end

function _panel_command(config::PrelimConfig, component::PanelComponent, lateral_command::Float64)
    η = clamp(lateral_command, -1.0, 1.0)
    abs(η) <= 1e-9 && return (alpha_deg=config.geometry.panel_alpha_deg, cant_deg=0.0, area_scale=1.0)
    if sign(η) != sign(component.side_sign)
        return (alpha_deg=config.geometry.panel_alpha_deg, cant_deg=0.0, area_scale=1.0)
    end
    return (
        alpha_deg=max(config.geometry.panel_alpha_deg - abs(η) * config.geometry.panel_lateral_deflection_deg, 0.0),
        cant_deg=90.0 * sign(component.side_sign),
        area_scale=1.0,
    )
end

function _panel_command(
    _config::PrelimConfig,
    component::PanelComponent,
    command::DifferentialPanelCommand,
)
    if component.side_sign < 0.0
        return (
            alpha_deg=command.left_alpha_deg,
            cant_deg=command.left_cant_deg,
            area_scale=command.left_area_scale,
        )
    end
    return (
        alpha_deg=command.right_alpha_deg,
        cant_deg=command.right_cant_deg,
        area_scale=command.right_area_scale,
    )
end

function aerodynamic_loads_3d(
    config::PrelimConfig,
    aero_case::CalibratedAeroCase,
    deployed::Bool,
    mach::Float64,
    gamma::Float64,
    lateral_command::Float64,
    panel_command_override::Union{Nothing, DifferentialPanelCommand}=nothing,
    trim_active::Bool=false,
)
    state = deployed ? aero_case.deployed_state : aero_case.body_state
    body_CL, body_CD, body_area = _component_cl_cd_area(state.body_link, state.body_alpha_rad, mach, gamma)
    CDA_total_m2 = body_CD * body_area
    CLA_vertical_m2 = body_CL * body_area
    CSA_m2 = 0.0
    if deployed && !isempty(state.panel_components)
        fixed_panel_cla_m2, fixed_panel_cda_m2 = _sum_component_cla_cda(state.panel_components, state.panel_alpha_rad, mach, gamma)
        CDA_total_m2 += fixed_panel_cda_m2
        CLA_vertical_m2 += fixed_panel_cla_m2
    end
    if deployed && trim_active && !isempty(state.trim_panel_components)
        command = panel_command_override === nothing ? neutral_panel_command(config) : panel_command_override
        for component in state.trim_panel_components
            cmd = _panel_command(config, component, command)
            α_panel_rad = deg2rad(cmd.alpha_deg)
            CL_panel, CD_panel, area_panel = _component_cl_cd_area(component.link, α_panel_rad, mach, gamma)
            effective_area = cmd.area_scale * area_panel
            cda_m2 = CD_panel * effective_area
            cla_m2 = CL_panel * effective_area
            cant_rad = deg2rad(cmd.cant_deg)
            cla_vertical_m2 = cla_m2 * cos(cant_rad)
            csa_m2 = cla_m2 * sin(cant_rad)
            CDA_total_m2 += cda_m2
            CLA_vertical_m2 += cla_vertical_m2
            CSA_m2 += csa_m2
        end
    end
    beta_eff = aero_case.mass_kg / CDA_total_m2
    return (
        state_label=state.label,
        CDA_m2=CDA_total_m2,
        CLA_m2=CLA_vertical_m2,
        CSA_m2=CSA_m2,
        beta_eff_kg_m2=beta_eff,
    )
end

function calibrate_aero_case(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    target_beta_high::Float64,
    target_beta_ratio::Float64,
)::CalibratedAeroCase
    body_link = build_body_link(config)
    body_alpha_rad = deg2rad(config.geometry.body_alpha_deg)
    panel_alpha_rad = deg2rad(config.geometry.panel_alpha_deg)
    trim_panel_alpha_rad = deg2rad(config.geometry.panel_alpha_deg)
    _, temperature_k = density_temperature(adapter, config.h0_m, config.theta0_rad, 0.0)
    reference_mach = mach_number(adapter, config.V0_mps, temperature_k)

    body_cla_ref_m2, body_cda_ref_m2 = _component_cla_cda(body_link, body_alpha_rad, reference_mach, config.planet.γ)
    mass_kg = target_beta_high * body_cda_ref_m2
    target_beta_low = target_beta_high / target_beta_ratio

    dummy_panel_links = build_panel_links(1.0, config)
    panel_cla_unit_m2, panel_cda_unit_m2 = isempty(dummy_panel_links) ?
        (0.0, 0.0) :
        _sum_component_cla_cda(dummy_panel_links, panel_alpha_rad, reference_mach, config.planet.γ)

    desired_deployed_cda_m2 = mass_kg / target_beta_low
    desired_panel_cda_m2 = max(desired_deployed_cda_m2 - body_cda_ref_m2, 0.0)
    panel_area_total_m2 = panel_cda_unit_m2 > 0.0 ? desired_panel_cda_m2 / panel_cda_unit_m2 : 0.0
    panel_links = build_panel_links(panel_area_total_m2, config)
    panel_cla_ref_m2, panel_cda_ref_m2 = isempty(panel_links) ?
        (0.0, 0.0) :
        _sum_component_cla_cda(panel_links, panel_alpha_rad, reference_mach, config.planet.γ)

    trim_panel_area_total_m2 = config.trim_panel_area_fraction_of_deployed * panel_area_total_m2
    trim_panel_links = build_trim_panel_links(trim_panel_area_total_m2, config)
    trim_panel_cla_ref_m2, trim_panel_cda_ref_m2 = isempty(trim_panel_links) ?
        (0.0, 0.0) :
        _sum_component_cla_cda(trim_panel_links, trim_panel_alpha_rad, reference_mach, config.planet.γ)

    deployed_cla_ref_m2 = body_cla_ref_m2 + panel_cla_ref_m2
    deployed_cda_ref_m2 = body_cda_ref_m2 + panel_cda_ref_m2

    achieved_beta_high = mass_kg / body_cda_ref_m2
    achieved_beta_low = mass_kg / deployed_cda_ref_m2
    achieved_beta_ratio = achieved_beta_high / achieved_beta_low
    panel_area_each_m2 = panel_area_total_m2 / config.geometry.panel_count
    trim_panel_area_each_m2 = trim_panel_area_total_m2 / config.geometry.panel_count
    body_state = AeroForceState(
        config.geometry.body_label,
        body_link,
        PanelComponent[],
        PanelComponent[],
        body_alpha_rad,
        panel_alpha_rad,
        trim_panel_alpha_rad,
        0.0,
        0.0,
        0.0,
        0.0,
    )
    deployed_state = AeroForceState(
        config.geometry.deployed_label,
        body_link,
        panel_links,
        trim_panel_links,
        body_alpha_rad,
        panel_alpha_rad,
        trim_panel_alpha_rad,
        panel_area_total_m2,
        panel_area_each_m2,
        trim_panel_area_total_m2,
        trim_panel_area_each_m2,
    )

    return CalibratedAeroCase(
        config.geometry.body_label,
        config.geometry.deployed_label,
        target_beta_high,
        target_beta_low,
        target_beta_ratio,
        achieved_beta_high,
        achieved_beta_low,
        achieved_beta_ratio,
        mass_kg,
        temperature_k,
        reference_mach,
        body_state,
        deployed_state,
        2.0 * config.geometry.base_radius_m,
        body_cda_ref_m2,
        body_cla_ref_m2,
        deployed_cda_ref_m2,
        deployed_cla_ref_m2,
        trim_panel_cda_ref_m2,
        trim_panel_cla_ref_m2,
    )
end

function derive_fixed_mass_beta_targets(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter;
    entry_mass_kg::Real,
    deployed_beta_target_kg_m2::Real,
)
    probe_case = calibrate_aero_case(config, adapter, 1.0, 1.0)
    body_cda_ref_m2 = probe_case.body_cda_ref_m2
    beta_high = Float64(entry_mass_kg) / body_cda_ref_m2
    beta_low = Float64(deployed_beta_target_kg_m2)
    beta_ratio = beta_high / beta_low
    return (
        beta_high_kg_m2=beta_high,
        beta_low_kg_m2=beta_low,
        beta_ratio=beta_ratio,
        body_cda_ref_m2=body_cda_ref_m2,
        probe_case=probe_case,
    )
end

function deployed_drag_skirt_equivalent_area(config::PrelimConfig)
    geom = config.geometry
    deployed_radius_m = 0.5 * geom.deployed_drag_surface_diameter_m
    stowed_radius_m = geom.base_radius_m
    skirt_height_m = geom.drag_skirt_height_m
    if !(deployed_radius_m > stowed_radius_m && skirt_height_m > 0.0)
        return 0.0
    end
    slant_height_m = hypot(skirt_height_m, deployed_radius_m - stowed_radius_m)
    return π * (stowed_radius_m + deployed_radius_m) * slant_height_m
end

function derive_fixed_mass_beta_targets_from_deployed_geometry(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter;
    entry_mass_kg::Real,
)
    body_link = build_body_link(config)
    body_alpha_rad = deg2rad(config.geometry.body_alpha_deg)
    panel_alpha_rad = deg2rad(config.geometry.panel_alpha_deg)
    _, temperature_k = density_temperature(adapter, config.h0_m, config.theta0_rad, 0.0)
    reference_mach = mach_number(adapter, config.V0_mps, temperature_k)

    body_cla_ref_m2, body_cda_ref_m2 = _component_cla_cda(body_link, body_alpha_rad, reference_mach, config.planet.γ)
    panel_area_total_m2 = deployed_drag_skirt_equivalent_area(config)
    panel_links = build_panel_links(panel_area_total_m2, config)
    panel_cla_ref_m2, panel_cda_ref_m2 = isempty(panel_links) ?
        (0.0, 0.0) :
        _sum_component_cla_cda(panel_links, panel_alpha_rad, reference_mach, config.planet.γ)

    deployed_cda_ref_m2 = body_cda_ref_m2 + panel_cda_ref_m2
    beta_high = Float64(entry_mass_kg) / body_cda_ref_m2
    beta_low = Float64(entry_mass_kg) / deployed_cda_ref_m2
    beta_ratio = beta_high / beta_low

    return (
        beta_high_kg_m2=beta_high,
        beta_low_kg_m2=beta_low,
        beta_ratio=beta_ratio,
        body_cda_ref_m2=body_cda_ref_m2,
        panel_cda_ref_m2=panel_cda_ref_m2,
        deployed_cda_ref_m2=deployed_cda_ref_m2,
        panel_area_total_m2=panel_area_total_m2,
        panel_area_each_m2=panel_area_total_m2 / config.geometry.panel_count,
    )
end

function aerodynamic_loads(
    aero_case::CalibratedAeroCase,
    deployed::Bool,
    mach::Float64,
    gamma::Float64,
    trim_active::Bool=false,
)
    state = deployed ? aero_case.deployed_state : aero_case.body_state
    body_cla_m2, body_cda_m2 = _component_cla_cda(state.body_link, state.body_alpha_rad, mach, gamma)
    panel_cla_m2 = 0.0
    panel_cda_m2 = 0.0
    if deployed && !isempty(state.panel_components)
        panel_cla_m2, panel_cda_m2 = _sum_component_cla_cda(state.panel_components, state.panel_alpha_rad, mach, gamma)
    end
    trim_cla_m2 = 0.0
    trim_cda_m2 = 0.0
    if deployed && trim_active && !isempty(state.trim_panel_components)
        trim_cla_m2, trim_cda_m2 = _sum_component_cla_cda(state.trim_panel_components, state.trim_panel_alpha_rad, mach, gamma)
    end
    cla_total_m2 = body_cla_m2 + panel_cla_m2 + trim_cla_m2
    cda_total_m2 = body_cda_m2 + panel_cda_m2 + trim_cda_m2
    beta_eff = aero_case.mass_kg / cda_total_m2
    return (; state_label=state.label, CLA_m2=cla_total_m2, CDA_m2=cda_total_m2, beta_eff_kg_m2=beta_eff)
end

function scale_aero_case_mass(
    aero_case::CalibratedAeroCase,
    mass_scale::Real,
)::CalibratedAeroCase
    scale = Float64(mass_scale)
    scale > 0.0 || throw(ArgumentError("mass_scale must be positive."))
    return CalibratedAeroCase(
        aero_case.body_label,
        aero_case.deployed_label,
        aero_case.target_beta_high * scale,
        aero_case.target_beta_low * scale,
        aero_case.target_beta_ratio,
        aero_case.achieved_beta_high * scale,
        aero_case.achieved_beta_low * scale,
        aero_case.achieved_beta_ratio,
        aero_case.mass_kg * scale,
        aero_case.reference_temperature_k,
        aero_case.reference_mach,
        aero_case.body_state,
        aero_case.deployed_state,
        aero_case.characteristic_length_m,
        aero_case.body_cda_ref_m2,
        aero_case.body_cla_ref_m2,
        aero_case.deployed_cda_ref_m2,
        aero_case.deployed_cla_ref_m2,
        aero_case.trim_only_cda_ref_m2,
        aero_case.trim_only_cla_ref_m2,
    )
end
