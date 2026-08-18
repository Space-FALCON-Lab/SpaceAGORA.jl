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
using StaticArrays

const _EDG_SM = SimulationModel
const _ODYSSEY_LIMIT_QDOT_W_CM2 = 0.15
const _ODYSSEY_LIMIT_Q_J_CM2 = 30.0
const _ODYSSEY_LIMIT_DYN_PRESS_PA = 0.5
const _ODYSSEY_RA_INITIAL_M = 28_559.615e3
const _ODYSSEY_HP_INITIAL_M = 77e3
const _ODYSSEY_EPOCH = (year=2001, month=11, day=6, hour=19, minute=0, second=0.0)

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
    )
    left = _EDG_SM.Link{0}(
        root=false,
        m=10.0,
        dims=MVector{3, Float64}(0.01, 3.76 / 2.0, 1.93),
        ref_area=(3.76 / 2.0) * 1.93,
        r=MVector{3, Float64}(0.0, -2.0, 0.0),
    )
    right = _EDG_SM.Link{0}(
        root=false,
        m=10.0,
        dims=MVector{3, Float64}(0.01, 3.76 / 2.0, 1.93),
        ref_area=(3.76 / 2.0) * 1.93,
        r=MVector{3, Float64}(0.0, 2.0, 0.0),
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
    p = _EDG_SM.ODEParams(n_sats=1, args=args)
    p.shared_buffers.et_start[] = _EDG_SM.ephemerides_time_seconds(args.initial_time, args.environment_model.ephemerides_model)
    u = ComponentVector(sc=[(pos=collect(pos), vel=collect(vel), mass=461.0, heat_loads=zeros(3))])
    return (; config, state, guidance, control, args, p, u, spacecraft)
end

@testset "Energy depletion GNC config and guidance" begin
    @test _ODYSSEY_LIMIT_QDOT_W_CM2 == 0.15
    @test _ODYSSEY_LIMIT_Q_J_CM2 == 30.0
    @test _ODYSSEY_LIMIT_DYN_PRESS_PA == 0.5
    @test _ODYSSEY_RA_INITIAL_M == 28_559.615e3
    @test _ODYSSEY_HP_INITIAL_M == 77e3

    cfg = _EDG_SM.AerobrakingEnergyDepletionConfig(
        guidance_modes=(:max_energy_depletion,),
        max_energy_submodes=(:heat_rate,),
    )
    @test cfg.guidance_modes == (:max_energy_depletion,)
    @test cfg.max_energy_submodes == (:heat_rate,)
    @test cfg.targeting_certification_samples == 9
    @test cfg.targeting_energy_order_tolerance_jkg == 1e-3
    @test cfg.targeting_heat_load_tolerance_j_cm2 == 1e-6
    @test_throws ArgumentError _EDG_SM.AerobrakingEnergyDepletionConfig(guidance_modes=(:not_a_mode,))
    @test_throws ArgumentError _EDG_SM.AerobrakingEnergyDepletionConfig(targeting_certification_samples=1)
    @test_throws ArgumentError _EDG_SM.AerobrakingEnergyDepletionConfig(targeting_energy_order_tolerance_jkg=-1.0)
    @test_throws ArgumentError _EDG_SM.AerobrakingEnergyDepletionConfig(targeting_heat_load_tolerance_j_cm2=-1.0)

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
    @test !isfinite(ctx_targeting.state.targeting_switch_s[1])
    @test !ctx_targeting.state.targeting_active[1]
    @test ctx_targeting.state.selected_mode[1] == :safe_low_drag
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

@testset "EDG-T candidate-map certification" begin
    candidate_times = [0.0, 1.0, 2.0, 3.0]

    ordered_energy = [10.0, 9.0, 8.0, 7.0]
    feasible_heat = [1.0, 2.0, 3.0, 4.0]
    ordered_evaluator(t) = begin
        idx = round(Int, t) + 1
        (energy_jkg=ordered_energy[idx], heat_load_j_cm2=feasible_heat[idx])
    end
    ordered = _EDG_SM.ControlHooks._edg_certify_targeting_candidates(
        candidate_times,
        ordered_evaluator;
        heat_load_limit_j_cm2=5.0,
        energy_order_tolerance_jkg=1e-6,
        heat_load_tolerance_j_cm2=1e-6,
    )
    @test ordered.times == candidate_times
    @test ordered.failure == :none

    nonmonotone_energy = [10.0, 9.0, 9.5, 8.0]
    nonmonotone_evaluator(t) = begin
        idx = round(Int, t) + 1
        (energy_jkg=nonmonotone_energy[idx], heat_load_j_cm2=feasible_heat[idx])
    end
    nonmonotone = _EDG_SM.ControlHooks._edg_certify_targeting_candidates(
        candidate_times,
        nonmonotone_evaluator;
        heat_load_limit_j_cm2=5.0,
        energy_order_tolerance_jkg=1e-6,
        heat_load_tolerance_j_cm2=1e-6,
    )
    @test nonmonotone.times == candidate_times[1:2]
    @test nonmonotone.failure == :energy_order
    @test nonmonotone.failure_time_s == candidate_times[3]

    infeasible_heat = [1.0, 2.0, 6.0, 4.0]
    heat_evaluator(t) = begin
        idx = round(Int, t) + 1
        (energy_jkg=ordered_energy[idx], heat_load_j_cm2=infeasible_heat[idx])
    end
    heat_limited = _EDG_SM.ControlHooks._edg_certify_targeting_candidates(
        candidate_times,
        heat_evaluator;
        heat_load_limit_j_cm2=5.0,
        energy_order_tolerance_jkg=1e-6,
        heat_load_tolerance_j_cm2=1e-6,
    )
    @test heat_limited.times == candidate_times[1:2]
    @test heat_limited.failure == :heat_load
    @test heat_limited.failure_time_s == candidate_times[3]
end
