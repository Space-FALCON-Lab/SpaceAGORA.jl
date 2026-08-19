if !isdefined(@__MODULE__, :REPO_ROOT)
    const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
end

if !isdefined(@__MODULE__, :SimulationModel)
    include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
    using .SimulationModel
end

using Test
using ComponentArrays
using LinearAlgebra
using SpecialFunctions
using StaticArrays

const _EDG_SM = SimulationModel
const _ODYSSEY_LIMIT_QDOT_W_CM2 = 0.15
const _ODYSSEY_LIMIT_Q_J_CM2 = 30.0
const _ODYSSEY_LIMIT_DYN_PRESS_PA = 0.5
const _ODYSSEY_RA_INITIAL_M = 9_800e3
const _ODYSSEY_HP_INITIAL_M = 100e3
const _ODYSSEY_EPOCH = (year=2001, month=11, day=6, hour=19, minute=0, second=0.0)
const _VEX_RA_INITIAL_M = 72_649e3
const _VEX_HP_INITIAL_M = 136e3
const _VEX_TARGET_RA_M = 72_480e3
const _VEX_EPOCH = (year=2014, month=5, day=19, hour=4, minute=0, second=0.0)

function _edg_test_context(;
    guidance_modes=(:max_energy_depletion,),
    max_energy_submodes=(:heat_rate, :structural_load, :heat_load),
    heat_load_switch_solver=:closed_form,
    target_apoapsis_radius_m=4.84e6,
    heat_rate_limit_w_cm2=_ODYSSEY_LIMIT_QDOT_W_CM2,
    heat_load_limit_j_cm2=1e-4,
    structural_load_limit_pa=_ODYSSEY_LIMIT_DYN_PRESS_PA,
)
    planet = _EDG_SM.Mars()
    planet.L_PI .= @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    rp = planet.Rp_e + _ODYSSEY_HP_INITIAL_M
    semi_major = 0.5 * (_ODYSSEY_RA_INITIAL_M + rp)
    periapsis_speed = sqrt(planet.μ * (2.0 / rp - 1.0 / semi_major))
    pos = SVector{3, Float64}(rp, 0.0, 0.0)
    vel = SVector{3, Float64}(0.0, periapsis_speed, 0.0)
    ic = _EDG_SM.CartesianInitialCondition(pos, vel)

    root = _EDG_SM.Link{0}(
        root=true,
        m=391.0,
        dims=MVector{3, Float64}(2.2, 2.6, 1.7),
        ref_area=2.2 * 1.7,
        reflection_coefficient=0.9,
    )
    left = _EDG_SM.Link{0}(
        root=false,
        m=10.0,
        dims=MVector{3, Float64}(0.01, 3.76 / 2.0, 1.93),
        ref_area=(3.76 / 2.0) * 1.93,
        r=MVector{3, Float64}(0.0, -2.0, 0.0),
        reflection_coefficient=0.9,
    )
    right = _EDG_SM.Link{0}(
        root=false,
        m=10.0,
        dims=MVector{3, Float64}(0.01, 3.76 / 2.0, 1.93),
        ref_area=(3.76 / 2.0) * 1.93,
        r=MVector{3, Float64}(0.0, 2.0, 0.0),
        reflection_coefficient=0.9,
    )
    spacecraft = _EDG_SM.SpacecraftModel(
        _EDG_SM.Joint[],
        [root, left, right],
        root,
        true,
        411.0,
        50.0,
        root.inertia,
        0,
        0,
        ic,
        1,
    )

    state = _EDG_SM.AerobrakingEnergyDepletionState(num_sats=1)
    config = _EDG_SM.AerobrakingEnergyDepletionConfig(
        guidance_modes=guidance_modes,
        max_energy_submodes=max_energy_submodes,
        heat_load_switch_solver=heat_load_switch_solver,
        controlled_panel_links=(2, 3),
        target_apoapsis_radius_m=target_apoapsis_radius_m,
        max_alpha_rad=pi / 2,
        min_alpha_rad=1e-4,
        heat_rate_limit_w_cm2=heat_rate_limit_w_cm2,
        heat_load_limit_j_cm2=heat_load_limit_j_cm2,
        structural_load_limit_pa=structural_load_limit_pa,
        planning_horizon_s=5_000.0,
        switch_recompute_interval_s=1.0,
    )
    guidance = _EDG_SM.AerobrakingEnergyDepletionGuidanceModel(config, state)
    control = _EDG_SM.AerobrakingEnergyDepletionControlModel(config, state)
    args = _EDG_SM.SimulationConfiguration(
        environment_model=_EDG_SM.EnvironmentModel(
            planet=planet,
            EI=160.0,
            density_model=_EDG_SM.ExponentialAtmosphereModel(1e-4, 100e3, 20e3; temperature_k=150.0),
            ephemerides_model=_EDG_SM.SimpleEphemeridesModel(),
            topography=false,
            topo_degree=0,
            topo_order=0,
            wind=false,
            thermal_model=_EDG_SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        ),
        dynamics_model=_EDG_SM.DynamicsModel([spacecraft], (_EDG_SM.InverseSquaredGravityModel(),)),
        guidance_model=_EDG_SM.GuidanceModel((guidance,), [1.0]),
        navigation_model=_EDG_SM.NavigationModel((), Float64[]),
        control_model=_EDG_SM.ControlModel((control,), [1.0]),
        initial_time=_EDG_SM.InitialTime(
            year=_ODYSSEY_EPOCH.year,
            month=_ODYSSEY_EPOCH.month,
            day=_ODYSSEY_EPOCH.day,
            hour=_ODYSSEY_EPOCH.hour,
            minute=_ODYSSEY_EPOCH.minute,
            second=_ODYSSEY_EPOCH.second,
        ),
    )
    p = _EDG_SM.ODEParams{1}(args=args)
    p.shared_buffers.et_start[] = _EDG_SM.ephemerides_time_seconds(args.initial_time, args.environment_model.ephemerides_model)
    u = ComponentVector(sc=[(pos=collect(pos), vel=collect(vel), mass=461.0, heat_loads=zeros(3))])
    return (; config, state, guidance, control, args, p, u, spacecraft)
end

@testset "VEX N-body target-energy correction" begin
    hooks = _EDG_SM.ControlHooks
    spice_path = joinpath(REPO_ROOT, "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")
    planet = _EDG_SM.Venus("", spice_path)
    periapsis_radius = planet.Rp_e + _VEX_HP_INITIAL_M
    nominal_target = hooks._edg_target_energy_from_apoapsis(
        planet,
        _VEX_TARGET_RA_M,
        periapsis_radius,
    )
    @test hooks._edg_corrected_target_energy_from_apoapsis(
        planet,
        _VEX_TARGET_RA_M,
        periapsis_radius,
        125.0,
    ) == nominal_target - 125.0

    ctx = _edg_test_context()
    epoch = _EDG_SM.InitialTime(
        year=_VEX_EPOCH.year,
        month=_VEX_EPOCH.month,
        day=_VEX_EPOCH.day,
        hour=_VEX_EPOCH.hour,
        minute=_VEX_EPOCH.minute,
        second=_VEX_EPOCH.second,
    )
    ephemerides_model = _EDG_SM.SpiceEphemeridesModel()
    sun_gravity = _EDG_SM.NBodyGravityModel(
        body_names=("Sun",),
        primary_body_name="Venus",
        planet=planet,
    )
    environment = _EDG_SM.EnvironmentModel(
        planet=planet,
        EI=200.0,
        density_model=_EDG_SM.ExponentialAtmosphereModel(1e-9, _VEX_HP_INITIAL_M, 15e3),
        ephemerides_model=ephemerides_model,
        thermal_model=_EDG_SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        topo_degree=0,
        topo_order=0,
        wind=false,
    )
    args_without_nbody = _EDG_SM.SimulationConfiguration(
        environment_model=environment,
        dynamics_model=_EDG_SM.DynamicsModel(
            [ctx.spacecraft],
            (_EDG_SM.InverseSquaredJ2GravityModel(),),
        ),
        guidance_model=ctx.args.guidance_model,
        navigation_model=ctx.args.navigation_model,
        control_model=ctx.args.control_model,
        initial_time=epoch,
    )
    args_with_nbody = _EDG_SM.SimulationConfiguration(
        environment_model=environment,
        dynamics_model=_EDG_SM.DynamicsModel(
            [ctx.spacecraft],
            (_EDG_SM.InverseSquaredJ2GravityModel(), sun_gravity),
        ),
        guidance_model=ctx.args.guidance_model,
        navigation_model=ctx.args.navigation_model,
        control_model=ctx.args.control_model,
        initial_time=epoch,
    )
    p_without_nbody = _EDG_SM.ODEParams{1}(args=args_without_nbody)
    p_with_nbody = _EDG_SM.ODEParams{1}(args=args_with_nbody)
    et = _EDG_SM.ephemerides_time_seconds(epoch, ephemerides_model)
    p_without_nbody.shared_buffers.et_start[] = et
    p_with_nbody.shared_buffers.et_start[] = et

    position = SVector{3, Float64}(_VEX_RA_INITIAL_M, 0.0, 0.0)
    velocity = SVector{3, Float64}(0.0, 1_000.0, 0.0)
    mass = 650.0
    acceleration_without_nbody = hooks._edg_prediction_gravity_acceleration(
        p_without_nbody,
        position,
        velocity,
        mass,
        0.0,
    )
    acceleration_with_nbody = hooks._edg_prediction_gravity_acceleration(
        p_with_nbody,
        position,
        velocity,
        mass,
        0.0,
    )
    expected_nbody = _EDG_SM.nbody_acceleration_ii_at_epoch(
        sun_gravity,
        position,
        p_with_nbody,
        et,
    )
    @test norm(expected_nbody) > 0.0
    @test acceleration_with_nbody - acceleration_without_nbody ≈ expected_nbody rtol=1e-12
end

@testset "Energy depletion GNC config and guidance" begin
    @test _ODYSSEY_LIMIT_QDOT_W_CM2 == 0.15
    @test _ODYSSEY_LIMIT_Q_J_CM2 == 30.0
    @test _ODYSSEY_LIMIT_DYN_PRESS_PA == 0.5
    @test _ODYSSEY_RA_INITIAL_M == 9_800e3
    @test _ODYSSEY_HP_INITIAL_M == 100e3

    cfg = _EDG_SM.AerobrakingEnergyDepletionConfig(
        guidance_modes=(:max_energy_depletion,),
        max_energy_submodes=(:heat_rate,),
    )
    @test cfg.guidance_modes == (:max_energy_depletion,)
    @test cfg.max_energy_submodes == (:heat_rate,)
    @test cfg.second_switch_reevaluation
    @test cfg.heat_load_security_mode
    @test_throws ArgumentError _EDG_SM.AerobrakingEnergyDepletionConfig(guidance_modes=(:not_a_mode,))

    ctx = _edg_test_context(guidance_modes=(:max_energy_depletion,))
    _EDG_SM.calcGuidanceEffect!(ctx.guidance, ctx.u, ctx.p, 0.0, 1)
    @test ctx.state.selected_mode[1] == :max_energy_depletion
    @test ctx.state.energy_bracketing_count[1] == 0
    @test !ctx.state.energy_bracketing_evaluated[1]

    ctx_targeting = _edg_test_context(guidance_modes=(:targeting, :max_energy_depletion))
    _EDG_SM.calcGuidanceEffect!(ctx_targeting.guidance, ctx_targeting.u, ctx_targeting.p, 0.0, 1)
    @test ctx_targeting.state.energy_bracketing_count[1] == 1
    @test ctx_targeting.state.energy_bracketing_evaluated[1]
    sc = ctx_targeting.u.sc[1]
    pos = SVector{3, Float64}(sc.pos)
    vel = SVector{3, Float64}(sc.vel)
    low_drag, max_energy = _EDG_SM.ControlHooks._edg_targeting_bracket_outcomes(
        ctx_targeting.config,
        ctx_targeting.p,
        ctx_targeting.spacecraft,
        pos,
        vel,
        Float64(sc.mass),
        0.0;
        heat_load_j_cm2=0.0,
        heat_rate_control=true,
        structural_control=true,
    )
    bracket_min, bracket_max = extrema((low_drag.energy_jkg, max_energy.energy_jkg))
    @test ctx_targeting.state.bracket_min_energy_jkg[1] ≈ bracket_min
    @test ctx_targeting.state.bracket_max_energy_jkg[1] ≈ bracket_max
    @test all(isfinite, max_energy.heat_load_switches_s)
    @test max_energy.heat_load_switches_s[1] < max_energy.heat_load_switches_s[2]

    ctx_fallback = _edg_test_context(guidance_modes=(:targeting, :max_energy_depletion), target_apoapsis_radius_m=1.0e12)
    _EDG_SM.calcGuidanceEffect!(ctx_fallback.guidance, ctx_fallback.u, ctx_fallback.p, 0.0, 1)
    @test ctx_fallback.state.selected_mode[1] == :max_energy_depletion
    @test !ctx_fallback.state.targeting_active[1]
end

@testset "Energy depletion control switches and constraints" begin
    for solver in (:closed_form, :tpbvp_integration)
        ctx = _edg_test_context(heat_load_switch_solver=solver)
        _EDG_SM.calcGuidanceEffect!(ctx.guidance, ctx.u, ctx.p, 0.0, 1)
        _EDG_SM.calcControlEffect!(ctx.control, ctx.u, ctx.p, 0.0, 1)
        switches = ctx.state.heat_load_switches_s[1]
        @test length(switches) == 2
        @test all(isfinite, switches)
        @test switches[1] < switches[2]
    end

    ctx_before_drag = _edg_test_context(heat_load_switch_solver=:closed_form)
    planet = ctx_before_drag.args.environment_model.planet
    ctx_before_drag.u.sc[1].pos .= (planet.Rp_e + 300e3, 0.0, 0.0)
    _EDG_SM.calcGuidanceEffect!(ctx_before_drag.guidance, ctx_before_drag.u, ctx_before_drag.p, 0.0, 1)
    _EDG_SM.calcControlEffect!(ctx_before_drag.control, ctx_before_drag.u, ctx_before_drag.p, 0.0, 1)
    @test ctx_before_drag.state.heat_load_switches_s[1] == (Inf, Inf)
    @test !ctx_before_drag.state.heat_load_switch_solved[1]

    ctx_once = _edg_test_context(heat_load_switch_solver=:closed_form)
    _EDG_SM.calcGuidanceEffect!(ctx_once.guidance, ctx_once.u, ctx_once.p, 0.0, 1)
    _EDG_SM.calcControlEffect!(ctx_once.control, ctx_once.u, ctx_once.p, 0.0, 1)
    first_switches = ctx_once.state.heat_load_switches_s[1]
    ctx_once.u.sc[1].heat_loads .= 100.0
    _EDG_SM.calcControlEffect!(ctx_once.control, ctx_once.u, ctx_once.p, 10.0, 1)
    @test ctx_once.state.heat_load_switches_s[1] == first_switches
    @test ctx_once.state.heat_load_switch_solved[1]
    @test ctx_once.state.last_switch_solve_t[1] == 0.0

    ctx_targeting = _edg_test_context(guidance_modes=(:targeting,), max_energy_submodes=(:heat_rate, :structural_load))
    ctx_targeting.state.selected_mode[1] = :targeting
    ctx_targeting.state.targeting_active[1] = true
    ctx_targeting.state.bracket_min_energy_jkg[1] = -1.0e7
    ctx_targeting.state.bracket_max_energy_jkg[1] = -9.0e6
    ctx_targeting.state.target_energy_jkg[1] = -9.5e6
    _EDG_SM.calcControlEffect!(ctx_targeting.control, ctx_targeting.u, ctx_targeting.p, 0.0, 1)
    @test isfinite(ctx_targeting.state.targeting_switch_s[1])
    @test ctx_targeting.state.targeting_active[1]
    @test ctx_targeting.state.selected_mode[1] == :targeting
    @test ctx_targeting.state.heat_load_switches_s[1] == (Inf, Inf)

    ctx_constraints = _edg_test_context(
        max_energy_submodes=(:heat_rate, :structural_load),
        heat_rate_limit_w_cm2=Inf,
        structural_load_limit_pa=0.5,
    )
    ctx_constraints.state.selected_mode[1] = :max_energy_depletion
    _EDG_SM.calcControlEffect!(ctx_constraints.control, ctx_constraints.u, ctx_constraints.p, 0.0, 1)
    @test ctx_constraints.state.last_alpha_rad[1] ==
        min(ctx_constraints.state.last_alpha_heat_rate_rad[1], ctx_constraints.state.last_alpha_structural_rad[1])

    hooks = _EDG_SM.ControlHooks
    sc_state = ctx_constraints.u.sc[1]
    pos = SVector{3, Float64}(sc_state.pos)
    vel = SVector{3, Float64}(sc_state.vel)
    env = hooks._edg_targeting_prediction_environment(ctx_constraints.p, pos, vel, 0.0)
    track = hooks._edg_closed_form_heat_load_trajectory(
        ctx_constraints.config,
        ctx_constraints.p,
        ctx_constraints.spacecraft,
        pos,
        vel,
        Float64(sc_state.mass),
        0.0,
        env,
    )
    raw_high = fill(ctx_constraints.config.max_alpha_rad, length(track.time))
    alpha_heat = hooks._edg_constrained_heat_load_alpha_profile(
        ctx_constraints.config, ctx_constraints.p, ctx_constraints.spacecraft,
        ctx_constraints.config.controlled_panel_links, track, raw_high;
        heat_rate_control=true, structural_control=false,
    )
    alpha_structural = hooks._edg_constrained_heat_load_alpha_profile(
        ctx_constraints.config, ctx_constraints.p, ctx_constraints.spacecraft,
        ctx_constraints.config.controlled_panel_links, track, raw_high;
        heat_rate_control=false, structural_control=true,
    )
    alpha_high = hooks._edg_constrained_heat_load_alpha_profile(
        ctx_constraints.config, ctx_constraints.p, ctx_constraints.spacecraft,
        ctx_constraints.config.controlled_panel_links, track, raw_high;
        heat_rate_control=true, structural_control=true,
    )
    @test alpha_high ≈ min.(alpha_heat, alpha_structural)
    @test !isdefined(hooks, :_EDG_TPBVP_HEAT_LOAD_PREDICTION_FACTOR)

    shooting_env = hooks._edg_targeting_prediction_environment(ctx_constraints.p, pos, vel, 0.0)
    shooting_coeffs = hooks._edg_heat_load_coefficients(
        ctx_constraints.config, ctx_constraints.p, ctx_constraints.spacecraft, shooting_env,
    )
    shooting_scale_height = hooks._edg_heat_load_scale_height(ctx_constraints.p)
    shooting_costates = SVector{3, Float64}(1.0e9, 0.0, 0.0)
    low_altitude_state = SVector{9, Float64}(pos..., vel..., shooting_costates...)
    high_altitude_pos = normalize(pos) * (planet.Rp_e + 120e3)
    high_altitude_state = SVector{9, Float64}(high_altitude_pos..., vel..., shooting_costates...)
    low_altitude_derivative, _, _ = hooks._edg_heat_load_shooting_stage(
        ctx_constraints.config, ctx_constraints.p, ctx_constraints.spacecraft,
        Float64(sc_state.mass), 0.0, 0.0, low_altitude_state,
        shooting_coeffs, shooting_scale_height, 0.05, ctx_constraints.config.max_alpha_rad;
        heat_rate_control=false, structural_control=false,
    )
    high_altitude_derivative, _, _ = hooks._edg_heat_load_shooting_stage(
        ctx_constraints.config, ctx_constraints.p, ctx_constraints.spacecraft,
        Float64(sc_state.mass), 0.0, 0.0, high_altitude_state,
        shooting_coeffs, shooting_scale_height, 0.05, ctx_constraints.config.max_alpha_rad;
        heat_rate_control=false, structural_control=false,
    )
    low_gravity = hooks._edg_prediction_gravity_acceleration(
        ctx_constraints.p, pos, vel, Float64(sc_state.mass), 0.0,
    )
    high_gravity = hooks._edg_prediction_gravity_acceleration(
        ctx_constraints.p, high_altitude_pos, vel, Float64(sc_state.mass), 0.0,
    )
    low_aero = SVector{3, Float64}(low_altitude_derivative[4:6]) - low_gravity
    high_aero = SVector{3, Float64}(high_altitude_derivative[4:6]) - high_gravity
    @test norm(low_aero) > norm(high_aero)

    root_q = copy(ctx_constraints.spacecraft.links[1].q)
    @test ctx_constraints.spacecraft.links[2].α == ctx_constraints.state.last_alpha_rad[1]
    @test ctx_constraints.spacecraft.links[3].α == ctx_constraints.state.last_alpha_rad[1]
    @test ctx_constraints.spacecraft.links[1].q == root_q

    infeasible_env = (
        dynamic_pressure=1.0e6,
        temperature=150.0,
        molecular_speed_ratio=10.0,
    )
    reference_drag_area = _EDG_SM.ControlHooks._energy_depletion_struct_drag_area(
        ctx_constraints.spacecraft,
        infeasible_env.temperature,
        infeasible_env.molecular_speed_ratio,
        ctx_constraints.config.controlled_panel_links,
        ctx_constraints.config.max_alpha_rad,
        ctx_constraints.config,
    )
    minimum_drag = infeasible_env.dynamic_pressure * _EDG_SM.ControlHooks._energy_depletion_struct_drag_area(
        ctx_constraints.spacecraft,
        infeasible_env.temperature,
        infeasible_env.molecular_speed_ratio,
        ctx_constraints.config.controlled_panel_links,
        ctx_constraints.config.min_alpha_rad,
        ctx_constraints.config,
    )
    @test minimum_drag > ctx_constraints.config.structural_load_limit_pa * reference_drag_area
    @test _EDG_SM.ControlHooks._energy_depletion_struct_load_root_alpha(
        ctx_constraints.config,
        infeasible_env,
        ctx_constraints.spacecraft,
        ctx_constraints.config.controlled_panel_links,
        ctx_constraints.config.max_alpha_rad,
    ) == ctx_constraints.config.min_alpha_rad

    force, torque = _EDG_SM.calcControlForceTorque(ctx_constraints.control, ctx_constraints.u.sc[1], ctx_constraints.p, 1, 0.0)
    @test force == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test torque == SVector{3, Float64}(0.0, 0.0, 0.0)
end

@testset "EDG Phase 1 legacy equation parity" begin
    hooks = _EDG_SM.ControlHooks
    ctx = _edg_test_context(max_energy_submodes=(:heat_rate, :structural_load))

    taf = 1.0
    rho = 2.2e-7
    temperature = 150.0
    speed_ratio = 17.0
    alpha = 0.61
    gamma = ctx.args.environment_model.planet.γ
    gas_constant = ctx.args.environment_model.planet.R
    s_sin = speed_ratio * sin(alpha)
    legacy_heat_rate = rho * (1e-4 * taf * gas_constant * temperature * sqrt(gas_constant * temperature / (2pi))) * (
        (speed_ratio^2 + gamma / (gamma - 1.0) - (gamma + 1.0) / (2.0 * (gamma - 1.0))) *
        (exp(-s_sin^2) + sqrt(pi) * s_sin * (1.0 + erf(s_sin))) -
        0.5 * exp(-s_sin^2)
    )
    @test hooks._energy_depletion_heat_rate_calc(
        taf,
        rho,
        temperature,
        temperature,
        gas_constant,
        gamma,
        speed_ratio,
        alpha,
    ) ≈ legacy_heat_rate rtol=2e-14

    heat_limit = 0.15
    alpha_heat = hooks._energy_depletion_heatrate_root_alpha(
        taf=taf,
        rho=rho,
        T_p=temperature,
        R=gas_constant,
        gamma=gamma,
        S=speed_ratio,
        max_alpha=pi / 2,
        min_alpha=1e-4,
        heat_rate_limit=heat_limit,
        alpha_past=0.2,
    )
    @test hooks._energy_depletion_heat_rate_calc(
        taf,
        rho,
        temperature,
        temperature,
        gas_constant,
        gamma,
        speed_ratio,
        alpha_heat,
    ) ≈ heat_limit - 1e-5 atol=2e-10

    sigma = 0.9
    function legacy_plate(alpha)
        ss = speed_ratio * sin(alpha)
        exponential = exp(-ss^2)
        error_term = 1.0 + erf(ss)
        cn = ((((2.0 - sigma) / sqrt(pi)) * ss + sigma / 2.0) * exponential +
              ((2.0 - sigma) * (ss^2 + 0.5) + sigma * sqrt(pi) * ss / 2.0) * error_term) / speed_ratio^2
        ca = sigma * cos(alpha) * (exponential + sqrt(pi) * ss * error_term) / (sqrt(pi) * speed_ratio)
        return cn * cos(alpha) - ca * sin(alpha), ca * cos(alpha) + cn * sin(alpha)
    end
    panel_area = 3.76 * 1.93
    bus_area = 2.2 * 1.7
    cl_panel, cd_panel = legacy_plate(alpha)
    cl_bus, cd_bus = legacy_plate(pi / 2)
    expected_cl = (cl_panel * panel_area + cl_bus * bus_area) / (panel_area + bus_area)
    expected_cd = (cd_panel * panel_area + cd_bus * bus_area) / (panel_area + bus_area)
    actual_cl, actual_cd = hooks._edg_legacy_spacecraft_aero_coefficients(
        ctx.spacecraft,
        speed_ratio,
        (2, 3),
        alpha,
    )
    @test actual_cl ≈ expected_cl rtol=2e-14
    @test actual_cd ≈ expected_cd rtol=2e-14

    times = [0.0, 100.0]
    f1, f2 = hooks._edg_closed_form_f1_f2(
        ctx.args.environment_model.planet,
        4350.0,
        deg2rad(-5.25),
        times,
        100.0,
    )
    @test f1 ≈ 4.6815 atol=1e-12
    @test f2[1] == 0.0
    @test f2[2] > 0.0
    @test hooks._edg_closed_form_density(ctx.args.environment_model.planet, 175e3) > 0.0
end

@testset "EDG Phase 2 and 3 scheduling parity" begin
    hooks = _EDG_SM.ControlHooks

    ctx_before_targeting = _edg_test_context(
        guidance_modes=(:targeting, :max_energy_depletion),
    )
    planet = ctx_before_targeting.args.environment_model.planet
    ctx_before_targeting.u.sc[1].pos .= (planet.Rp_e + 300e3, 0.0, 0.0)
    _EDG_SM.calcGuidanceEffect!(ctx_before_targeting.guidance, ctx_before_targeting.u, ctx_before_targeting.p, 0.0, 1)
    @test !ctx_before_targeting.state.energy_bracketing_evaluated[1]
    @test ctx_before_targeting.state.energy_bracketing_count[1] == 0
    @test ctx_before_targeting.state.selected_mode[1] == :inactive

    targeting_grid = hooks._edg_targeting_prediction_time_grid(10.4)
    @test first(targeting_grid) == 0.0
    @test last(targeting_grid) == 10.4
    @test maximum(diff(targeting_grid)) <= 1.0

    fit = hooks._edg_closed_form_correction_fit(planet)
    @test fit[1:4] == (4350.0, 450.0, -5.25, 2.25)
    @test fit[5] == (4.560850, -1.715925, -0.375150, 0.376025, -0.254500)
    @test fit[6][1:4] == (-4.74628450, -3.69824500, -6.06348525, 1.07688200)

    rp = planet.Rp_e + _ODYSSEY_HP_INITIAL_M
    a = 0.5 * (_ODYSSEY_RA_INITIAL_M + rp)
    e = 1.0 - rp / a
    interface_radius = planet.Rp_e + 160e3
    nu_exit = acos((a * (1.0 - e^2) / interface_radius - 1.0) / e)
    nu_outbound = 0.5 * nu_exit
    semi_latus_rectum = a * (1.0 - e^2)
    radius_outbound = semi_latus_rectum / (1.0 + e * cos(nu_outbound))
    pos_outbound = SVector{3, Float64}(
        radius_outbound * cos(nu_outbound),
        radius_outbound * sin(nu_outbound),
        0.0,
    )
    velocity_scale = sqrt(planet.μ / semi_latus_rectum)
    vel_outbound = velocity_scale * SVector{3, Float64}(-sin(nu_outbound), e + cos(nu_outbound), 0.0)
    outbound_duration = hooks._edg_drag_passage_duration(
        ctx_before_targeting.config,
        ctx_before_targeting.p,
        pos_outbound,
        vel_outbound,
        461.0,
    )
    expected_duration = (
        hooks._edg_mean_anomaly_from_true(nu_exit, e) -
        hooks._edg_mean_anomaly_from_true(nu_outbound, e)
    ) / sqrt(planet.μ / a^3)
    @test outbound_duration ≈ expected_duration rtol=1e-10
    @test outbound_duration < 1_000.0

    ctx_targeting_once = _edg_test_context(
        guidance_modes=(:targeting,),
        max_energy_submodes=(:heat_rate, :structural_load),
    )
    state = ctx_targeting_once.state
    state.selected_mode[1] = :targeting
    state.targeting_active[1] = true
    state.target_energy_jkg[1] = -9.5e6
    _EDG_SM.calcControlEffect!(ctx_targeting_once.control, ctx_targeting_once.u, ctx_targeting_once.p, 0.0, 1)
    first_switch = state.targeting_switch_s[1]
    first_solve_time = state.last_switch_solve_t[1]
    @test isfinite(first_switch)
    state.target_energy_jkg[1] = -9.0e6
    _EDG_SM.calcControlEffect!(ctx_targeting_once.control, ctx_targeting_once.u, ctx_targeting_once.p, 20.0, 1)
    @test state.targeting_switch_s[1] == first_switch
    @test state.last_switch_solve_t[1] == first_solve_time
end
