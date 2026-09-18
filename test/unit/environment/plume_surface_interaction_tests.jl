module PlumeInteractionTests
using Test, SpaceAGORA, StaticArrays, LinearAlgebra, ComponentArrays
using DiffEqCallbacks: SavedValues
const SM = SpaceAGORA.SimulationModel
const SC = SM.SimulationCallbacks
const PSI = SM.DynamicEffectors.PlumeSurfaceInteraction
const Z3 = SVector(0.0, 0.0, 0.0)
const EYE = SMatrix{3,3,Float64}(I)
const PLANET = SM.make_no_gram_planet(:earth)
control(values) = (; actuators=(; thrust_n=Float64.(values)))

function state_environment(height; q=nothing, position_axis=SVector(0.0,0.0,1.0), terrain_lat=0.0)
    pos = (PLANET.Rp_e + height) * position_axis
    x = SM.StateSample(pos, Z3, 500.0; q_ib=q)
    # A deliberately unrelated geodetic latitude exposes accidental use of it
    # for terrain, whose contract requires planetocentric coordinates from pos_pp.
    pf = SM.PlanetFrameSample(EYE, pos, Z3, height, terrain_lat, 0.0)
    return x, SM.EnvironmentSample(PLANET; planet_frame=pf)
end

@testset "validated plume construction and pure geometry" begin
    cfg = PlumeSurfaceConfig()
    values = Tuple(getfield(cfg, name) for name in fieldnames(typeof(cfg)))
    @test PlumeSurfaceConfig(values...) == cfg
    for (i, name) in enumerate(fieldnames(typeof(cfg)))
        @test_throws ArgumentError PlumeSurfaceConfig(; NamedTuple{(name,)}((NaN,))...)
        @test_throws ArgumentError PlumeSurfaceConfig(Base.setindex(values, Inf, i)...)
    end
    for kwargs in ((; nozzle_exit_radius_m=0.0), (; nozzle_offset_m=-1.0),
            (; plume_half_angle_deg=0.0), (; plume_half_angle_deg=90.0),
            (; friction_coefficient=-1.0), (; threshold_shear_pa=-1.0),
            (; erosion_efficiency=-1.0), (; ground_effect_max_fraction=1.1),
            (; ejecta_speed_min_mps=20.0, ejecta_speed_max_mps=10.0))
        @test_throws ArgumentError PlumeSurfaceConfig(;kwargs...)
    end
    @test_throws ArgumentError PlumeSurfaceState(0)
    @test_throws ArgumentError PlumeSurfaceInteractionModel(control([1,2]); num_sats=1)
    @test_throws ArgumentError PlumeSurfaceInteractionModel(control([-1]))
    @test_throws ArgumentError PlumeSurfaceInteractionModel(control([NaN]))
    @test_throws ArgumentError PlumeSurfaceInteractionModel((;))
    for f in (plume_surface_footprint, plume_quantities, plume_ground_effect_force)
        @test_throws ArgumentError f(cfg, -1.0, 1.0)
        @test_throws ArgumentError f(cfg, 1.0, -1.0)
        @test_throws ArgumentError f(cfg, Inf, 1.0)
    end
    @test plume_quantities(cfg, 0.0, 2.0).erosion_kg_s == 0.0
    @test plume_quantities(cfg, 1000.0, cfg.max_height_m + 1).ground_effect_n == 0.0
    @test plume_ground_effect_force(cfg, 1000.0, 0.0) ≈ 1000cfg.ground_effect_max_fraction
    @test plume_ground_effect_force(cfg, 1000.0, 2cfg.nozzle_exit_radius_m * cfg.ground_effect_cutoff) == 0.0

    cfg0 = PlumeSurfaceConfig(nozzle_offset_m=0.0)
    cfg1 = PlumeSurfaceConfig(nozzle_offset_m=1.0)
    model0 = PlumeSurfaceInteractionModel(control([1000.0]); config=cfg0)
    model1 = PlumeSurfaceInteractionModel(control([1000.0]); config=cfg1)
    x, env = state_environment(2.0)
    f0, torque0 = wrench(model0, x, env, 0.0)
    f1, torque1 = wrench(model1, x, env, 0.0)
    @test norm(f1) > norm(f0) > 0
    @test f1[3] > 0
    @test torque0 == torque1 == Z3
    @test norm(cross(-SVector(0.0,0.0,1.0) * cfg1.nozzle_offset_m, f1)) == 0.0
    for (q, axis) in ((SVector(0.0,0.0,0.0,1.0), SVector(0.0,0.0,1.0)),
                      (SVector(0.0,0.0,0.0,1.0), SVector(1.0,0.0,0.0)))
        up, upe = state_environment(2.0; q=q, position_axis=axis)
        @test wrench(model1, up, upe, 0.0) == (Z3, Z3)
        @test PSI._plume_sample(model1, up, upe, 0.0, 1).height_m == Inf
    end
    down, downe = state_environment(2.0; q=SVector(1.0,0.0,0.0,0.0))
    @test wrench(model1, down, downe, 0.0) == (f1, torque1)
    model2 = PlumeSurfaceInteractionModel(control([1000.0, 2000.0]); config=cfg1)
    @test_throws ArgumentError wrench(model2, x, env, 0.0)
    @test wrench(model2, x, env, 0.0, 2)[1] ≈ 2wrench(model2, x, env, 0.0, 1)[1]
    before = deepcopy(model2.state)
    for t in (0.0, 100.0, -1.0, 0.01)
        wrench(model2, x, env, t, 2)
        SM.DynamicEffectors.wrench_caching!(model2, x, env, t, nothing, 1)
    end
    @test all(name -> isequal(getfield(before,name), getfield(model2.state,name)), fieldnames(PlumeSurfaceState))

    # Pole sample must use latitude +90 from its Cartesian position, not the
    # deliberately supplied geodetic latitude zero.
    grid = DEMGrid(fill(1.0f0,2,2), 80.0,90.0,-10.0,10.0; reference_radius_m=PLANET.Rp_e)
    terrain = DEMTerrainModel([grid]; reference_radius_m=PLANET.Rp_e)
    terrain_model = PlumeSurfaceInteractionModel(control([1000]), terrain)
    sample = PSI._plume_sample(terrain_model, x, env, 0.0, 1)
    @test sample.height_m ≈ 1.0
end

mutable struct HeldControl <: SpaceAGORA.AbstractControlEffectorModel
    actuators::NamedTuple{(:thrust_n,),Tuple{Vector{Float64}}}
    initial::Vector{Float64}
    jump::Bool
end
function SM.ControlHooks.calcControlEffect!(model::HeldControl, u::ComponentVector,
        p::SM.ODEParams, t::Float64, i::Int64)
    model.actuators.thrust_n[i] = model.initial[i] * (model.jump && t >= 1.0 ? 2.0 : 1.0)
    return nothing
end
SM.ControlHooks.calcControlForceTorque(::HeldControl, u::AbstractVector,
    p::SM.ODEParams, i::Int64, t::Float64) = (Z3,Z3)
SM.ControlHooks.calcControlMassFlowRate(::HeldControl, u::AbstractVector,
    p::SM.ODEParams, i::Int64, t::Float64) = 0.0
struct OscillatingForce <: SpaceAGORA.AbstractForceTorqueModel end
SM.wrench(::OscillatingForce, x::SM.StateSample, env::SM.EnvironmentSample, t::Float64) =
    (SVector(500sin(200t), 0.0, 0.0), Z3)

function make_config(; dt=0.4, duration=2.0, jump=false, oscillate=false, cfg=nothing, height=20.0)
    cfg === nothing && (cfg=PlumeSurfaceConfig(nozzle_offset_m=0.0,
        threshold_shear_pa=0.0, ejecta_speed_min_mps=50.0, ejecta_speed_max_mps=50.0,
        ground_effect_max_fraction=0.0))
    values = [1000.0,2000.0]
    held = HeldControl((;thrust_n=copy(values)), values, jump)
    plume = PlumeSurfaceInteractionModel(held; config=cfg)
    spacecraft = SM.SpacecraftModel[]
    for i in 1:2
        root = SM.Link(root=true, m=500.0, ref_area=1.0)
        ic = SM.CartesianInitialCondition([PLANET.Rp_e+height,0.0,0.0], [-1.0,0.0,0.0])
        push!(spacecraft, SM.SpacecraftModel(SM.Joint[],[root],root,true,500.0,0.0,root.inertia,0,0,ic,100+i))
    end
    args = SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_time=duration,
            data_rate=0.25, num_steps_to_save=9),
        environment_model=SM.make_no_gram_environment(planet=PLANET),
        dynamics_model=SM.DynamicsModel(spacecraft, oscillate ? (plume,OscillatingForce()) : (plume,)),
        guidance_model=SM.GuidanceModel(guidance_effectors=(),guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(),navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(held,),control_rates=[1.0]),
        initial_time=SM.InitialTime(year=2020),
        integration_tolerances=SM.IntegrationTolerances(dt_max_orbit=dt,
            reltol_orbit=1e-10,abstol_orbit=1e-11),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))
    return args,plume
end
function run_saved(args; isolate=true)
    saved = SavedValues(Float64,SM.SaveData)
    fields = SC.default_save_fields(args)
    callback = SC.get_data_saving_callback(2,args,fields,saved)
    sol = run_simulation(args; isolate_state=isolate,return_solution=true,
        visualization=false,extra_callbacks=(callback,))
    return sol,saved
end

@testset "accepted-step plume integral, saves, reset and two-spacecraft isolation" begin
    args,original = make_config(jump=true)
    rates = [plume_quantities(original.config,F,20.0).erosion_kg_s for F in [1000.0,2000.0]]
    @test rates[2] ≈ 2rates[1]
    sol,saved = run_saved(args)
    @test string(sol.retcode) == "Success"
    used = sol.prob.p.args.dynamics_model.dynamic_effectors[1]
    @test used !== original
    @test used.control === only(sol.prob.p.args.control_model.control_effectors)
    @test all(iszero,original.state.eroded_kg)
    @test used.state.eroded_kg ≈ 3rates rtol=1e-10
    @test length(saved.t) > 2
    for (t,row) in zip(saved.t,saved.saveval)
        factor = min(t,1.0) + 2max(t-1.0,0.0)
        @test row[:plume_eroded_kg] ≈ factor*rates rtol=1e-10 atol=1e-10
        @test row[:plume_erosion_kg_s] ≈ (t < 1 ? rates : 2rates) rtol=1e-10
        @test row[:plume_height_m] ≈ fill(20.0-t,2) atol=1e-7
    end
    # Reuse the caller-owned effector: the new run must reset its diagnostics.
    args2,shared = make_config()
    run_saved(args2;isolate=false)
    first_mass = copy(shared.state.eroded_kg)
    run_saved(args2;isolate=false)
    @test shared.state.eroded_kg ≈ first_mass
    @test first_mass ≈ 2rates rtol=1e-10

    # The actual engine reuses one callback set across checkpoint segments.
    # Initialization at each segment must retain, not reset, this run's mass.
    mktempdir() do directory
        segment_args,_ = make_config()
        settings=SM.SimulationSettings(results=false,verbose=false,generate_plots=false,
            normalize=false,save_csv=false,checkpoint_enabled=true,
            checkpoint_interval_s=0.5,checkpoint_directory=directory)
        segment_args=SM.SimConfig._with_configuration(segment_args;simulation_settings=settings)
        segmented=run_simulation(segment_args;return_solution=true,
            return_solver_metadata=true,visualization=false)
        @test length(segmented.solver_trace) == 4
        segment_model=segmented.solution.prob.p.args.dynamics_model.dynamic_effectors[1]
        @test segment_model.state.eroded_kg ≈ 2rates rtol=1e-10
    end

    # Positive ground effect must reach the real translational RHS, not only
    # the diagnostic callback. Primary thrust remains zero in this fixture.
    positive_args,_ = make_config(duration=0.2,height=2.0,cfg=PlumeSurfaceConfig())
    positive,_ = run_saved(positive_args)
    @test positive.u[end].sc[2].vel[1] > positive.u[end].sc[1].vel[1] > -1.0
    @test positive.prob.p.args.dynamics_model.dynamic_effectors[1].state.ground_effect_n[1] > 0

    # Its force is already rejected by the existing backbone dispatch policy,
    # before any second-order callback sampling is attempted.
    backbone = SM.SimConfig._with_configuration(args2;
        control_model=SM.ControlModel(control_effectors=(),control_rates=Float64[]),
        solver_config=SM.SolverConfig(solver_mode=:gravity_backbone_split))
    reason = SpaceAGORA.SimulationEngine._gravity_backbone_reject_reason(backbone)
    @test occursin("does not support PlumeSurfaceInteractionModel", reason)
    @test_throws ArgumentError run_simulation(backbone;visualization=false)

    # Oscillatory forcing creates rejected trial stages. Constant-rate erosion
    # still has an exact time integral because speed is clamped and threshold=0.
    reject_args,reject_model = make_config(dt=2.0,oscillate=true)
    rejected,_ = run_saved(reject_args)
    @test rejected.stats.nreject > 0
    actual = rejected.prob.p.args.dynamics_model.dynamic_effectors[1]
    @test actual.state.eroded_kg ≈ 2rates rtol=1e-10
end

@testset "diagnostic mass converges with accepted-step refinement" begin
    cfg=PlumeSurfaceConfig(nozzle_offset_m=0.0, threshold_shear_pa=0.5,
        ejecta_speed_min_mps=50.0,ejecta_speed_max_mps=50.0,
        ground_effect_max_fraction=0.0)
    duration=2.0; height=6.0; count=20000
    # Independent fine midpoint time quadrature along the known straight-line
    # trajectory; no numerical ODE path is used to construct this reference.
    reference=sum(plume_quantities(cfg,1000.0,height-(k-0.5)*duration/count).erosion_kg_s
        for k in 1:count)*duration/count
    @test reference > 0
    errors=Float64[]
    for dt in (0.5,0.1,0.02)
        args,_=make_config(dt=dt,cfg=cfg,height=height)
        sol,_=run_saved(args)
        mass=sol.prob.p.args.dynamics_model.dynamic_effectors[1].state.eroded_kg[1]
        push!(errors,abs(mass-reference))
    end
    @test errors[3] < errors[1]
    @test errors[3] < 1e-3*reference
end
end # module PlumeInteractionTests
