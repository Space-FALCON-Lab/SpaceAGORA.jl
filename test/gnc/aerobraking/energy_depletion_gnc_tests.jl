if !isdefined(@__MODULE__, :REPO_ROOT)
    const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
end

if !isdefined(@__MODULE__, :SimulationModel)
    using SpaceAGORA
    const SimulationModel = SpaceAGORA.SimulationModel
end

using Test
using ComponentArrays
using LinearAlgebra
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
const _VEX_EPOCH = (year=2014, month=5, day=19, hour=4, minute=0, second=0.0)

function _edg_test_context(;
    guidance_modes=(:max_energy_depletion,),
    max_energy_submodes=(:heat_rate, :structural_load, :heat_load),
    heat_load_switch_solver=:closed_form,
    target_apoapsis_radius_m=4.84e6,
    heat_rate_limit_w_cm2=_ODYSSEY_LIMIT_QDOT_W_CM2,
    heat_load_limit_j_cm2=1e-4,
    structural_load_limit_pa=_ODYSSEY_LIMIT_DYN_PRESS_PA,
    density_reference_kg_m3=1e-4,
    gravity_models=(_EDG_SM.InverseSquaredGravityModel(),),
    density_model=_EDG_SM.ExponentialAtmosphereModel(density_reference_kg_m3, 100e3, 20e3; temperature_k=150.0),
    planning_horizon_s=5_000.0,
    heat_load_security_mode=true,
)
    planet = _EDG_SM.Mars()
    planet.L_PI .= @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    rp = planet.Rp_e + _ODYSSEY_HP_INITIAL_M
    semi_major = 0.5 * (_ODYSSEY_RA_INITIAL_M + rp)
    periapsis_speed = sqrt(planet.μ * (2.0 / rp - 1.0 / semi_major))
    pos = SVector{3, Float64}(rp, 0.0, 0.0)
    vel = SVector{3, Float64}(0.0, periapsis_speed, 0.0)
    ic = _EDG_SM.CartesianInitialCondition(pos, vel)

    root = _EDG_SM.Link(
        root=true,
        m=391.0,
        dims=MVector{3, Float64}(2.2, 2.6, 1.7),
        ref_area=2.2 * 1.7,
        reflection_coefficient=0.9,
    )
    left = _EDG_SM.Link(
        root=false,
        m=10.0,
        dims=MVector{3, Float64}(0.01, 3.76 / 2.0, 1.93),
        ref_area=(3.76 / 2.0) * 1.93,
        r=MVector{3, Float64}(0.0, -2.0, 0.0),
        reflection_coefficient=0.9,
    )
    right = _EDG_SM.Link(
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
    # A run reports the EDG control clock when it starts; direct calls here use the
    # configured 1 s control period from the origin.
    state.control_clock_period_s = 1.0
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
        planning_horizon_s=planning_horizon_s,
        heat_load_security_mode=heat_load_security_mode,
        switch_recompute_interval_s=1.0,
    )
    guidance = _EDG_SM.AerobrakingEnergyDepletionGuidanceModel(config, state)
    control = _EDG_SM.AerobrakingEnergyDepletionControlModel(config, state)
    args = _EDG_SM.SimulationConfiguration(
        environment_model=_EDG_SM.EnvironmentModel(
            planet=planet,
            EI=160.0,
            density_model=density_model,
            ephemerides_model=_EDG_SM.SimpleEphemeridesModel(),
            topography=false,
            topo_degree=0,
            topo_order=0,
            wind=false,
            thermal_model=_EDG_SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        ),
        dynamics_model=_EDG_SM.DynamicsModel([spacecraft], gravity_models),
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
    p.shared_buffers.et_start[] = _EDG_SM.ephemerides_epoch_seconds(args.environment_model.ephemerides_model, args.initial_time)
    u = ComponentVector(sc=[(pos=collect(pos), vel=collect(vel), mass=461.0, heat_loads=zeros(3))])
    return (; config, state, guidance, control, args, p, u, spacecraft)
end

# A prediction context with the control clock of the context's EDG state, for direct
# calls of prediction kernels that take no EDG state.
_edg_clocked(ctx) = _EDG_SM.EDGServices._borrow_context(ctx.p, _EDG_SM.EDGServices._edg_state_clock(ctx.state))

@testset "VEX N-body target-energy correction" begin
    # The prediction's gravity query adds the run's N-body acceleration; the target
    # energy correction itself is tested in SpaceAGORAEDG.
    gravity = _EDG_SM.EDGServices._edg_prediction_gravity_query
    spice_path = joinpath(REPO_ROOT, "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")
    planet = _EDG_SM.Venus("", spice_path)

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
    p_without_nbody = _EDG_SM.ODEParams(n_sats=1, args=args_without_nbody)
    p_with_nbody = _EDG_SM.ODEParams(n_sats=1, args=args_with_nbody)
    et = _EDG_SM.ephemerides_epoch_seconds(ephemerides_model, epoch)
    p_without_nbody.shared_buffers.et_start[] = et
    p_with_nbody.shared_buffers.et_start[] = et

    position = SVector{3, Float64}(_VEX_RA_INITIAL_M, 0.0, 0.0)
    velocity = SVector{3, Float64}(0.0, 1_000.0, 0.0)
    mass = 650.0
    acceleration_without_nbody = gravity(
        p_without_nbody,
        position,
        velocity,
        mass,
        0.0,
    )
    acceleration_with_nbody = gravity(
        p_with_nbody,
        position,
        velocity,
        mass,
        0.0,
    )
    expected_nbody = _EDG_SM.DynamicEffectors.PerturbationEffectors.nbody_acceleration_ii_at_epoch(
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
    _EDG_SM.update_guidance!(ctx.guidance, ctx.u, ctx.p, 0.0, 1)
    @test ctx.state.selected_mode[1] == :max_energy_depletion
    @test ctx.state.energy_bracketing_count[1] == 0
    @test !ctx.state.energy_bracketing_evaluated[1]

    ctx_targeting = _edg_test_context(density_reference_kg_m3=1e-7, guidance_modes=(:targeting, :max_energy_depletion))
    _EDG_SM.update_guidance!(ctx_targeting.guidance, ctx_targeting.u, ctx_targeting.p, 0.0, 1)
    @test ctx_targeting.state.energy_bracketing_count[1] == 1
    @test ctx_targeting.state.energy_bracketing_evaluated[1]
    # The tiny fixture heat budget cannot certify a targeting family.
    @test isnan(ctx_targeting.state.bracket_min_energy_jkg[1])
    @test ctx_targeting.state.targeting_diagnostics[1].status == :no_admissible_family

    ctx_fallback = _edg_test_context(density_reference_kg_m3=1e-7, guidance_modes=(:targeting, :max_energy_depletion), target_apoapsis_radius_m=1.0e12)
    _EDG_SM.update_guidance!(ctx_fallback.guidance, ctx_fallback.u, ctx_fallback.p, 0.0, 1)
    @test ctx_fallback.state.selected_mode[1] == :safe_low_drag
    @test !ctx_fallback.state.targeting_active[1]
end

@testset "Energy depletion control switches and constraints" begin
    for solver in (:closed_form, :tpbvp_integration)
        ctx = _edg_test_context(heat_load_switch_solver=solver)
        _EDG_SM.update_guidance!(ctx.guidance, ctx.u, ctx.p, 0.0, 1)
        _EDG_SM.update_control!(ctx.control, ctx.u, ctx.p, 0.0, 1)
        switches = ctx.state.heat_load_switches_s[1]
        @test length(switches) == 2
        @test all(isfinite, switches)
        @test switches[1] < switches[2]
    end

    ctx_before_drag = _edg_test_context(heat_load_switch_solver=:closed_form)
    planet = ctx_before_drag.args.environment_model.planet
    ctx_before_drag.u.sc[1].pos .= (planet.Rp_e + 300e3, 0.0, 0.0)
    _EDG_SM.update_guidance!(ctx_before_drag.guidance, ctx_before_drag.u, ctx_before_drag.p, 0.0, 1)
    _EDG_SM.update_control!(ctx_before_drag.control, ctx_before_drag.u, ctx_before_drag.p, 0.0, 1)
    @test ctx_before_drag.state.heat_load_switches_s[1] == (Inf, Inf)
    @test !ctx_before_drag.state.heat_load_switch_solved[1]

    ctx_once = _edg_test_context(heat_load_switch_solver=:closed_form)
    _EDG_SM.update_guidance!(ctx_once.guidance, ctx_once.u, ctx_once.p, 0.0, 1)
    _EDG_SM.update_control!(ctx_once.control, ctx_once.u, ctx_once.p, 0.0, 1)
    first_switches = ctx_once.state.heat_load_switches_s[1]
    ctx_once.u.sc[1].heat_loads .= 100.0
    _EDG_SM.update_control!(ctx_once.control, ctx_once.u, ctx_once.p, 10.0, 1)
    @test ctx_once.state.heat_load_switches_s[1] == first_switches
    @test ctx_once.state.heat_load_switch_solved[1]
    @test ctx_once.state.last_switch_solve_t[1] == 0.0

    ctx_targeting = _edg_test_context(density_reference_kg_m3=1e-7, guidance_modes=(:targeting,), max_energy_submodes=(:heat_rate, :structural_load))
    ctx_targeting.state.selected_mode[1] = :targeting
    ctx_targeting.state.targeting_active[1] = true
    ctx_targeting.state.bracket_min_energy_jkg[1] = -1.0e7
    ctx_targeting.state.bracket_max_energy_jkg[1] = -9.0e6
    ctx_targeting.state.target_energy_jkg[1] = -9.5e6
    ctx_targeting.state.targeting_switch_s[1] = 10.0 # Plan supplied by guidance.
    _EDG_SM.update_control!(ctx_targeting.control, ctx_targeting.u, ctx_targeting.p, 0.0, 1)
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
    _EDG_SM.update_control!(ctx_constraints.control, ctx_constraints.u, ctx_constraints.p, 0.0, 1)
    @test ctx_constraints.state.last_alpha_rad[1] ==
        min(ctx_constraints.state.last_alpha_heat_rate_rad[1], ctx_constraints.state.last_alpha_structural_rad[1])

    @test !isdefined(_EDG_SM.ControlHooks, :_EDG_TPBVP_HEAT_LOAD_PREDICTION_FACTOR)

    root_q = copy(ctx_constraints.spacecraft.links[1].q)
    @test ctx_constraints.spacecraft.links[2].α == ctx_constraints.state.last_alpha_rad[1]
    @test ctx_constraints.spacecraft.links[3].α == ctx_constraints.state.last_alpha_rad[1]
    @test ctx_constraints.spacecraft.links[1].q == root_q

    sc1 = ctx_constraints.u.sc[1]
    x_ctx = _EDG_SM.StateSample(SVector{3, Float64}(sc1.pos), SVector{3, Float64}(sc1.vel), sc1.mass; spacecraft_index=1)
    force, torque = _EDG_SM.control_wrench(ctx_constraints.control, x_ctx, ctx_constraints.p, 0.0)
    @test force == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test torque == SVector{3, Float64}(0.0, 0.0, 0.0)
end

@testset "EDG Phase 2 and 3 scheduling parity" begin
    ctx_before_targeting = _edg_test_context(
        guidance_modes=(:targeting, :max_energy_depletion),
    )
    planet = ctx_before_targeting.args.environment_model.planet
    ctx_before_targeting.u.sc[1].pos .= (planet.Rp_e + 300e3, 0.0, 0.0)
    _EDG_SM.update_guidance!(ctx_before_targeting.guidance, ctx_before_targeting.u, ctx_before_targeting.p, 0.0, 1)
    @test !ctx_before_targeting.state.energy_bracketing_evaluated[1]
    @test ctx_before_targeting.state.energy_bracketing_count[1] == 0
    @test ctx_before_targeting.state.selected_mode[1] == :inactive

    ctx_targeting_once = _edg_test_context(
        density_reference_kg_m3=1e-7,
        guidance_modes=(:targeting,),
        max_energy_submodes=(:heat_rate, :structural_load),
    )
    state = ctx_targeting_once.state
    state.selected_mode[1] = :targeting
    state.targeting_active[1] = true
    state.target_energy_jkg[1] = -9.5e6
    state.targeting_switch_s[1] = 21.0
    state.last_switch_solve_t[1] = 0.0
    for t in (0.0,20.0,22.0)
        # Even a supplied energy change cannot make control replace the
        # guidance-owned plan. The public control hook only applies it.
        state.target_energy_jkg[1] = -9.0e6
        _EDG_SM.update_control!(ctx_targeting_once.control,ctx_targeting_once.u,ctx_targeting_once.p,t,1)
        @test state.targeting_switch_s[1] == 21.0
        @test state.last_switch_solve_t[1] == 0.0
    end
end

@testset "EDG no-switch plans are cached without a low-drag pulse" begin
    ctx = _edg_test_context(
        heat_load_limit_j_cm2=1e12,
        heat_rate_limit_w_cm2=Inf,
        structural_load_limit_pa=Inf,
    )
    _EDG_SM.update_guidance!(ctx.guidance, ctx.u, ctx.p, 0.0, 1)
    _EDG_SM.update_control!(ctx.control, ctx.u, ctx.p, 30.0, 1)
    @test ctx.state.heat_load_switches_s[1] == (30.0, 30.0)
    @test ctx.state.heat_load_switch_solved[1]
    @test ctx.state.last_alpha_rad[1] == ctx.config.max_alpha_rad
    for t in (30.0 + 1/3, 31.0, 35.0)
        _EDG_SM.update_control!(ctx.control, ctx.u, ctx.p, t, 1)
        @test ctx.state.last_alpha_rad[1] == ctx.config.max_alpha_rad
        @test ctx.state.heat_load_switches_s[1] == (30.0, 30.0)
        @test ctx.state.last_switch_solve_t[1] == 30.0
    end
    # A cached no-switch decision belongs only to this atmospheric passage.
    original_position = copy(ctx.u.sc[1].pos)
    ctx.u.sc[1].pos .= (ctx.args.environment_model.planet.Rp_e + 300e3, 0.0, 0.0)
    _EDG_SM.update_control!(ctx.control, ctx.u, ctx.p, 40.0, 1)
    @test !ctx.state.heat_load_switch_solved[1]
    @test !ctx.state.heat_load_drag_passage_active[1]
    @test ctx.state.heat_load_switches_s[1] == (Inf, Inf)
    ctx.u.sc[1].pos .= original_position
    _EDG_SM.update_control!(ctx.control, ctx.u, ctx.p, 41.0, 1)
    @test ctx.state.heat_load_switch_solved[1]
    @test ctx.state.heat_load_switches_s[1] == (41.0, 41.0)
    @test ctx.state.last_switch_solve_t[1] == 41.0
    @test ctx.state.last_alpha_rad[1] == ctx.config.max_alpha_rad
end

@testset "EDG callback clock follows a nonzero solve start" begin
    ctx=_edg_test_context()
    ctx.state.control_clock_period_s=NaN   # set again by the solve's start below
    cb=only(_EDG_SM.SimulationCallbacks.get_control_callbacks(1,ctx.args))
    de=_EDG_SM.SimulationCallbacks
    # End before the first periodic affect so this isolates initialization,
    # including its default u_modified! semantics.
    prob=de.ODEProblem((du,u,p,t)->fill!(du,0.0),[1.0],(12.345,12.445),ctx.p)
    sol=de.solve(prob,de.Tsit5();callback=cb,save_everystep=false)
    @test ctx.state.control_clock_origin_s==12.345
    @test ctx.state.control_clock_period_s==ctx.args.control_model.control_rates[1]
    @test sol.u[end]==[1.0]
    @test _EDG_SM.EDGServices._edg_targeting_control_clock(_edg_clocked(ctx)).origin==12.345
end

@testset "Retained EDG compatibility configuration" begin
 cfg = _EDG_SM.AerobrakingEnergyDepletionConfig()
 @test cfg.targeting_certification_samples == 9
 @test cfg.targeting_energy_order_tolerance_jkg == 1e-3
 @test cfg.targeting_heat_load_tolerance_j_cm2 == 1e-6
 @test_throws ArgumentError _EDG_SM.AerobrakingEnergyDepletionConfig(targeting_certification_samples=1)
 @test_throws ArgumentError _EDG_SM.AerobrakingEnergyDepletionConfig(targeting_energy_order_tolerance_jkg=-1)
 @test_throws ArgumentError _EDG_SM.AerobrakingEnergyDepletionConfig(targeting_heat_load_tolerance_j_cm2=-1)
end
