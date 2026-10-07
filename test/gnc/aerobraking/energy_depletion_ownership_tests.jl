
module EDGOwnershipTests
using Test, SpaceAGORA, StaticArrays, LinearAlgebra
import .._edg_test_context
const S = SpaceAGORA.SimulationModel
const A = S.EDGAlgorithms
const P = S.EDGServices
const C = S.ControlHooks
const G = S.GuidanceHooks
const root = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(root, "test", "contracts", "edg_ownership_checks.jl"))
state_values(s) = Tuple(copy(getfield(s, f)) for f in fieldnames(typeof(s)))
geometry(sc) = [(r=Tuple(l.r), q=Tuple(l.q), alpha=l.α, a=Tuple(l.aᵇ), b=Tuple(l.bᵇ),
                 inertia=Tuple(l.inertia), force=Tuple(l.net_force), torque=Tuple(l.net_torque)) for l in sc.links]
decision(x, t=0.0, i=1; links=x.control.aoa_effector.controlled_panel_links) =
    A.control_decision!(x.config, x.state, links, x.u, x.p, t, i)

# Retained pre-extraction public hook sequence. Its numerical callees are also
# checked against the matched-parent fingerprint and verbatim move inventory.
function reference_control!(model, u, p, t, i)
    state = model.state
    1 <= i <= length(state.selected_mode) || return nothing
    config = model.config
    if state.selected_mode[i] == :inactive
        state.selected_mode[i] = (:max_energy_depletion in config.guidance_modes) ? :max_energy_depletion : :safe_low_drag
    end
    sc = C._edg_control_sat_state(u, i)
    env = C._edg_environment_state(u, p, Float64(t), i)
    spacecraft = p.args.dynamics_model.spacecraft[i]
    pos, vel, mass = C._edg_control_pos_vel_mass(sc)
    heat_load = C._edg_max_heat_load_for_links(sc, model.aoa_effector.controlled_panel_links)
    C._edg_recompute_switches!(model, p, env, spacecraft, pos, vel, mass, heat_load, Float64(t), i)
    heat_load_low_drag_active = C._edg_heat_load_low_drag_active(model, Float64(t), i)
    base_alpha = C._edg_base_alpha(model, Float64(t), i)
    alpha = C._edg_command_alpha!(model, p, u, env, spacecraft, base_alpha, heat_load, heat_load_low_drag_active, i)
    C._apply_solar_panel_aoa!(model.aoa_effector, spacecraft, alpha)
    return nothing
end

@testset "Typed EDG owners and compatibility" begin
    @test S.AerobrakingEnergyDepletionConfig === G.AerobrakingEnergyDepletionConfig
    @test S.AerobrakingEnergyDepletionState === G.AerobrakingEnergyDepletionState
    @test S.AerobrakingEnergyDepletionGuidanceModel === G.AerobrakingEnergyDepletionGuidanceModel
    @test S.AerobrakingEnergyDepletionControlModel === C.AerobrakingEnergyDepletionControlModel
    @test S.calcGuidanceEffect! === G.calcGuidanceEffect!
    @test S.calcControlEffect! === C.calcControlEffect!
    @test G._edg_sat_state === C._edg_control_sat_state === P._edg_control_sat_state
    @test G._edg_pos_vel_mass === C._edg_control_pos_vel_mass === P._edg_control_pos_vel_mass
    for (name, path) in EDGOwnershipChecks.EXPECTED
        symbol = Symbol(name)
        owner = endswith(path, "/services.jl") ? P : A
        @test isdefined(owner, symbol)
        if isdefined(C, symbol) && name ∉ ("_edg_recompute_switches!", "_edg_base_alpha", "_edg_command_alpha!", "_edg_heat_load_low_drag_active")
            @test getfield(C, symbol) === getfield(owner, symbol)
        end
    end
    gnc_modules = [G, C, S.NavigationHooks]
    @test isempty(EDGOwnershipChecks.runtime_violations(A, P, gnc_modules))
    # A same-name function in another loaded module is not an alias.
    foreign = Module(gensym(:EDGForeignWitness))
    Core.eval(foreign, :(_edg_prediction_time_grid(x) = x))
    @test any(message -> occursin("separate EDG binding _edg_prediction_time_grid", message),
        EDGOwnershipChecks.runtime_violations(A, P, [gnc_modules; foreign]))
    # Qualified extensions are legal Julia, but do not belong to the EDG owner.
    Core.eval(foreign, :(const Algorithms = $A))
    before_methods = collect(methods(A._edg_prediction_time_grid))
    try
        Core.eval(foreign, :(Algorithms._edg_prediction_time_grid(::Val{:edg_foreign_witness}) = nothing))
        @test any(message -> occursin("foreign method on _edg_prediction_time_grid", message),
            Base.invokelatest(EDGOwnershipChecks.runtime_violations, A, P, gnc_modules))
    finally
        for method in setdiff(Base.invokelatest(() -> collect(methods(A._edg_prediction_time_grid))), before_methods)
            Base.delete_method(method)
        end
    end
    @test isempty(EDGOwnershipChecks.runtime_violations(A, P, gnc_modules))
    @test !isdefined(A, :ControlHooks)
    @test !isdefined(A, :SolarPanelAngleOfAttackControlModel)
    @test !isdefined(A, :_apply_solar_panel_aoa!)
end

@testset "Decision, telemetry and panel application match retained hook sequence" begin
    for mode in (:inactive, :safe_low_drag, :targeting, :max_energy_depletion),
        t in (prevfloat(10.0), 10.0, nextfloat(10.0), 20.0, nextfloat(20.0))
        x = _edg_test_context(heat_load_limit_j_cm2=30.0)
        x.state.selected_mode[1] = mode
        x.state.targeting_active[1] = mode == :targeting
        x.state.targeting_switch_s[1] = 10.0
        x.state.heat_load_switches_s[1] = (10.0,20.0)
        x.state.heat_load_switch_solved[1] = true
        old = deepcopy(x); hook = deepcopy(x)
        before = geometry(x.spacecraft); loads = copy(x.u.sc[1].heat_loads)
        result = decision(x, t)
        @test result.apply
        @test result.command isa S.CommandTypes.AerobrakingControlCommand
        @test geometry(x.spacecraft) == before
        @test x.u.sc[1].heat_loads == loads
        @test result.diagnostics.heat_load_j_cm2 == x.state.last_heat_load_j_cm2[1]
        @test isequal(result.diagnostics.targeting_switch_s, x.state.targeting_switch_s[1])
        C._apply_solar_panel_aoa!(x.control.aoa_effector, x.spacecraft, result.command.alpha_command)
        reference_control!(old.control, old.u, old.p, t, 1)
        @test S.calcControlEffect!(hook.control, hook.u, hook.p, t, 1) === nothing
        @test isequal(state_values(x.state), state_values(old.state))
        @test isequal(state_values(hook.state), state_values(old.state))
        @test geometry(x.spacecraft) == geometry(old.spacecraft) == geometry(hook.spacecraft)
    end
end

@testset "Public hooks preserve measured links, fallback and cached precedence" begin
    for (links, measured, expected_alpha) in (((2,), 3.0, pi/2), ((3,), 12.0, 1e-4))
        x = _edg_test_context(max_energy_submodes=(:heat_load,), heat_load_limit_j_cm2=10.0)
        x.u.sc[1].heat_loads .= [0.0, 3.0, 12.0]
        x.state.selected_mode[1] = :max_energy_depletion
        x.state.heat_load_switch_solved[1] = true
        x.state.heat_load_drag_passage_active[1] = true
        x.state.heat_load_switches_s[1] = (Inf, Inf)
        model = S.AerobrakingEnergyDepletionControlModel(x.config, x.state;
            aoa_effector=S.SolarPanelAngleOfAttackControlModel(links))
        before = geometry(x.spacecraft)
        @test x.config.controlled_panel_links == (2, 3)
        @test S.calcControlEffect!(model, x.u, x.p, 1.0, 1) === nothing
        @test x.state.last_heat_load_j_cm2[1] == measured
        @test x.state.last_alpha_rad[1] == expected_alpha
        @test x.spacecraft.links[only(links)].α == expected_alpha
        @test geometry(x.spacecraft)[links == (2,) ? 3 : 2] ==
              before[links == (2,) ? 3 : 2]
        @test x.u.sc[1].heat_loads == [0.0, 3.0, 12.0]
    end
    for modes in ((:targeting, :max_energy_depletion), (:targeting,))
        x = _edg_test_context(guidance_modes=modes)
        x.u.sc[1].pos .= (x.args.environment_model.planet.Rp_e + 300e3, 0.0, 0.0)
        before = geometry(x.spacecraft)
        @test S.calcGuidanceEffect!(x.guidance, x.u, x.p, 0.0, 1) === nothing
        @test x.state.selected_mode[1] == (length(modes) == 2 ? :max_energy_depletion : :safe_low_drag)
        @test x.state.safe_low_drag[1] == (length(modes) == 1)
        @test !x.state.targeting_active[1]
        @test !x.state.energy_bracketing_evaluated[1]
        @test x.state.energy_bracketing_count[1] == 0
        @test geometry(x.spacecraft) == before
    end
    x = _edg_test_context()
    x.u.sc[1].pos .= (x.args.environment_model.planet.Rp_e + 300e3, 0.0, 0.0)
    x.state.selected_mode[1] = :targeting
    x.state.targeting_active[1] = true
    x.state.targeting_switch_s[1] = 10.0
    @test decision(x, 8.0).switch_action == :cached
    @test x.state.targeting_switch_s[1] == 10.0
end

@testset "Targeting solve retains heat accumulated after entry bracketing" begin
    # Bounded scenario from the independent S10c differential. Start with zero
    # heat, establish a reachable target, then accumulate heat before control.
    options = (guidance_modes=(:targeting, :max_energy_depletion),
        max_energy_submodes=(:heat_rate, :structural_load, :heat_load),
        heat_rate_limit_w_cm2=Inf, structural_load_limit_pa=Inf,
        heat_load_limit_j_cm2=30.0,
        density_model=S.ExponentialAtmosphereModel(1e-8, 100e3, 20e3; temperature_k=150.0),
        planning_horizon_s=400.0)
    base = _edg_test_context(; options...)
    sc = base.u.sc[1]
    r0, v0 = SVector{3,Float64}(sc.pos), SVector{3,Float64}(sc.vel)
    duration = A._edg_drag_passage_duration(base.config, base.p, r0, v0, Float64(sc.mass))
    target = A._edg_targeting_outcome_with_heat_load(base.config, base.p, base.spacecraft,
        r0, v0, Float64(sc.mass), 0.0, 0.373 * (duration + 1.0), 0.0;
        heat_rate_control=true, structural_control=true)
    @test 10.0 < target.heat_load_j_cm2 < 20.0
    for accumulated in (0.0, 20.0)
        x = _edg_test_context(; options..., target_apoapsis_radius_m=target.apoapsis_radius_m)
        S.calcControlEffect!(x.control, x.u, x.p, 0.0, 1)
        S.calcGuidanceEffect!(x.guidance, x.u, x.p, 0.0, 1)
        @test x.state.targeting_active[1]
        @test x.state.selected_mode[1] == :targeting
        @test x.state.energy_bracketing_evaluated[1]
        @test !isfinite(x.state.targeting_switch_s[1])
        x.u.sc[1].heat_loads .= [0.0, accumulated, 0.0]
        @test S.calcControlEffect!(x.control, x.u, x.p, 1.0, 1) === nothing
        @test x.state.last_heat_load_j_cm2[1] == accumulated
        @test x.state.last_switch_solve_t[1] == 1.0
        @test x.state.targeting_active[1] == (accumulated == 0.0)
        @test x.state.selected_mode[1] == (accumulated == 0.0 ? :targeting : :max_energy_depletion)
        @test isfinite(x.state.targeting_switch_s[1]) == (accumulated == 0.0)
        @test x.u.sc[1].heat_loads == [0.0, accumulated, 0.0]
    end
end

@testset "Switch endpoints, sentinels and panel selection" begin
    x = _edg_test_context(max_energy_submodes=(:heat_load,), heat_load_limit_j_cm2=30.0)
    x.state.selected_mode[1] = :max_energy_depletion
    x.state.heat_load_switch_solved[1] = true
    x.state.heat_load_switches_s[1] = (10.,20.)
    for (t, active) in ((prevfloat(10.),false),(10.,true),(nextfloat(10.),true),(20.,true),(nextfloat(20.),false))
        @test A._edg_heat_load_low_drag_active(x.config,x.state,t,1) == active
        @test A._edg_base_alpha(x.config,x.state,t,1) == (active ? x.config.min_alpha_rad : x.config.max_alpha_rad)
    end
    x.state.selected_mode[1] = :targeting; x.state.targeting_active[1] = true; x.state.targeting_switch_s[1] = 10.
    @test A._edg_base_alpha(x.config,x.state,prevfloat(10.),1) == x.config.max_alpha_rad
    @test A._edg_base_alpha(x.config,x.state,10.,1) == x.config.min_alpha_rad
    x.state.selected_mode[1] = :max_energy_depletion
    x.state.heat_load_switches_s[1] = (Inf,Inf)
    x.u.sc[1].heat_loads .= [0.,1.,90.]
    before = geometry(x.spacecraft)
    result = decision(x,0.;links=(2,))
    @test x.config.controlled_panel_links == (2,3)
    @test result.switch_action == :cached
    @test result.diagnostics.heat_load_j_cm2 == 1.
    @test result.command.alpha_command == x.config.max_alpha_rad
    C._apply_solar_panel_aoa!(S.SolarPanelAngleOfAttackControlModel((2,)),x.spacecraft,result.command.alpha_command)
    @test geometry(x.spacecraft)[3] == before[3]
    @test A._edg_max_heat_load_for_links(x.u.sc[1],x.config.controlled_panel_links) == 90.
    @test decision(x,0.;links=(2,3)).command.alpha_command == x.config.min_alpha_rad
    @test A._edg_max_heat_load_for_links((heat_loads=[NaN,-1.,2.],),(1,2,3,8)) == 2.
    @test A._edg_max_heat_load_for_links((pos=zeros(3),),(2,)) == 0.
    x.state.selected_mode[1]=:targeting; x.state.targeting_active[1]=true; x.state.targeting_switch_s[1]=Inf
    x.state.target_energy_jkg[1]=NaN
    @test decision(x,5.).switch_action == :solved
    @test x.state.targeting_switch_s[1] == 5. + 0.5*x.config.planning_horizon_s
end

mutable struct QueryDensity <: S.AbstractDensityModel
    events::Vector{Any}
    fail::Bool
end
function S.EnvironmentModels.getDensity(m::QueryDensity,h::Float64,lat::Float64,lon::Float64,t::Float64,wind::Bool,p)
    push!(m.events,(:density,h,lat,lon,t,wind))
    m.fail && throw(ArgumentError("density witness"))
    return 1e-9,150.,SVector(60.,-35.,4.)
end
struct QueryFrames <: S.AbstractEphemeridesModel
    events::Vector{Any}
end
S.EphemeridesModels.ephemerides_requires_spice(::QueryFrames) = false
function S.EphemeridesModels.planet_frame_lpi(planet,t::Float64,m::QueryFrames)
    push!(m.events,(:frame,t))
    return SMatrix{3,3,Float64}(I)
end
function spy_context(; kwargs...)
    events=Any[]; density=QueryDensity(events,false)
    x=_edg_test_context(; density_model=density,kwargs...)
    old=x.args.environment_model
    fields=(; (f=>getfield(old,f) for f in fieldnames(typeof(old)))...)
    env=S.EnvironmentModel(; merge(fields,(ephemerides_model=QueryFrames(events),wind=true))...)
    args=S.SimConfig._with_configuration(x.args;environment_model=env)
    p=S.ODEParams(n_sats=1,args=args);p.shared_buffers.et_start[]=12345.
    return (;x...,args,p,events,density)
end
@testset "Lazy projection, query coordinates, time, wind and failure order" begin
    x=spy_context(max_energy_submodes=(:heat_rate,),heat_rate_limit_w_cm2=Inf)
    for i in (0,2)
        prior=state_values(x.state)
        r=decision(x,0.,i)
        @test !r.apply && r.switch_action == :invalid_index && r.diagnostics === nothing
        @test S.calcControlEffect!(x.control,x.u,x.p,0.,i) === nothing
        @test S.calcGuidanceEffect!(x.guidance,x.u,x.p,0.,i) === nothing
        @test isequal(state_values(x.state),prior)
    end
    @test isempty(x.events)
    x.state.energy_bracketing_evaluated[1]=true
    A._edg_run_target_energy_bracketing!(x.config,x.state,x.u,x.p,7.,1)
    @test isempty(x.events)
    x.u.sc[1].pos .= [3.0e6,1.0e6,1.2e6]
    r=SVector{3,Float64}(x.u.sc[1].pos);v=SVector{3,Float64}(x.u.sc[1].vel)
    env=P._edg_targeting_prediction_environment(x.p,r,v,7.)
    @test length(x.events)==3
    @test x.events[1:2]==[(:frame,12352.),(:frame,12352.)]
    lla=S.FrameTransforms.rtolatlong(r,x.args.environment_model.planet)
    @test x.events[3]==(:density,lla[1],lla[2],lla[3],7.,true)
    lat,lon=lla[2:3]
    east=SVector(-sin(lon),cos(lon),0.)
    north=SVector(-sin(lat)*cos(lon),-sin(lat)*sin(lon),cos(lat))
    up=SVector(cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat))
    expected=v-cross(x.args.environment_model.planet.ω,r)-(60.0*east-35.0*north+4.0*up)
    @test env.vel_pp_rw ≈ expected rtol=1e-13
    @test !isapprox(env.vel_pp_rw,v-cross(x.args.environment_model.planet.ω,r)+(60.0*east-35.0*north+4.0*up))
    empty!(x.events)
    @test P._edg_sample_prediction_atmosphere(x.p,-2.,7.) == (1e-9,150.)
    @test x.events == [(:density,0.,0.,0.,7.,true)]
    empty!(x.events)
    current_env=P._edg_environment_state(x.u,x.p,7.,1)
    @test current_env.speed ≈ norm(expected) rtol=1e-13
    @test x.events == [(:frame,12352.),(:density,lla[1],lla[2],lla[3],7.,true)]
    empty!(x.events)
    x.state.selected_mode[1]=:targeting;x.state.targeting_active[1]=true;x.state.targeting_switch_s[1]=2.
    @test decision(x,7.).switch_action==:cached
    @test length(x.events)==2 # only current control environment, no predictor
    empty!(x.events);x.density.fail=true;x.state.selected_mode[1]=:inactive
    before=geometry(x.spacecraft)
    @test_throws ArgumentError decision(x,7.)
    @test x.state.selected_mode[1]==:max_energy_depletion
    @test geometry(x.spacecraft)==before
    @test length(x.events)==2
    empty!(x.events);x.density.fail=false
    x.state.selected_mode[1]=:targeting;x.state.targeting_active[1]=true;x.state.targeting_switch_s[1]=Inf
    x.u.sc[1].pos .= (x.args.environment_model.planet.Rp_e+300e3,0.,0.)
    @test decision(x,8.).switch_action==:outside_pass
    @test x.state.targeting_switch_s[1]==Inf
    @test length(x.events)==2
end

@testset "Structured/flat inputs, missing mass and spacecraft index" begin
    x=_edg_test_context(max_energy_submodes=(:heat_rate,),heat_rate_limit_w_cm2=Inf)
    sc=x.u.sc[1];r,v,m=P._edg_control_pos_vel_mass(sc)
    @test P._edg_control_pos_vel_mass(collect(sc)[1:7]) == (r,v,m)
    @test isequal(P._edg_control_pos_vel_mass(collect(sc)[1:6]),(r,v,NaN))
    @test isequal(P._edg_control_pos_vel_mass((pos=collect(r),vel=collect(v))),(r,v,NaN))
    @test P._edg_control_sat_state((sc=[sc,sc],),2)===sc
    @test P._edg_control_sat_state(collect(sc),1)==collect(sc)
    @test A._edg_predict_mass(x.spacecraft,NaN) == A._edg_predict_mass(x.spacecraft,Inf)
    @test A._edg_predict_mass(x.spacecraft,321.)==321.
    structured=decision(x)
    y=_edg_test_context(max_energy_submodes=(:heat_rate,),heat_rate_limit_w_cm2=Inf)
    flat=A.control_decision!(y.config,y.state,(2,3),collect(y.u.sc[1])[1:7],y.p,0.,1)
    @test structured.command.alpha_command==flat.command.alpha_command
    @test isequal(state_values(x.state),state_values(y.state))
end

@testset "Architecture negative controls" begin
    sources=EDGOwnershipChecks.source_map(root)
    @test isempty(EDGOwnershipChecks.violations(sources))
    function refused(path,extra)
        bad=copy(sources);bad[path]=get(bad,path,"")*"\n"*extra
        return !isempty(EDGOwnershipChecks.violations(bad))
    end
    @test refused("src/gnc/control/other.jl","function _edg_prediction_time_grid(t)\n    return [t]\nend")
    @test refused(EDGOwnershipChecks.OWNER*"targeting.jl","const reverse = ControlHooks")
    @test refused(EDGOwnershipChecks.OWNER*"targeting.jl","rotate_link(link, axis, angle)")
    @test refused(EDGOwnershipChecks.OWNER*"targeting.jl","link.α = 0.0")
    @test refused("src/gnc/guidance/target_energy_bracketing.jl","_control_module()._edg_prediction_time_grid(1.0)")
    bad=copy(sources)
    path="src/gnc/control/targeting_control.jl"
    bad[path]=replace(bad[path],"    return EDGAlgorithms._edg_base_alpha"=>"    t += 1.0\n    return EDGAlgorithms._edg_base_alpha")
    @test !isempty(EDGOwnershipChecks.violations(bad))
    for definition in (
        "_edg_prediction_time_grid(t) = [t]",
        "    function _edg_prediction_time_grid(t)\n        [t]\n    end",
        "@noinline function _edg_prediction_time_grid(t)\n    [t]\nend",
        "function EDGAlgorithms._edg_prediction_time_grid(t)\n    [t]\nend",
        "control_decision!(args...) = nothing",
    )
        @test refused("src/gnc/control/other.jl", definition)
        @test refused("src/core/other.jl", definition)
    end
    @test refused(path, "_edg_base_alpha(model, t, i) = t")
    bad = copy(sources)
    bad[path] = replace(bad[path], "model.config, model.state, t, i)" =>
        "model.config, model.state, i, t)")
    @test !isempty(EDGOwnershipChecks.violations(bad))
    for include_call in (
        "include(joinpath(dirname(dirname(dirname(@__DIR__))), \"control\", \"struct_load_control.jl\"))",
        "Base.include(@__MODULE__, joinpath(dirname(@__DIR__), \"control\", \"other.jl\"))",
    )
        @test EDGOwnershipChecks.has_control_include(include_call)
        @test refused(EDGOwnershipChecks.OWNER * "algorithms.jl", include_call)
    end
    @test !EDGOwnershipChecks.has_control_include("# include(\"control.jl\")")
    @test !EDGOwnershipChecks.has_control_include("message = \"include(control) is forbidden\"")
    aggregator = sources[EDGOwnershipChecks.OWNER * "algorithms.jl"]
    @test !isempty(EDGOwnershipChecks.aggregator_violations(
        replace(aggregator, "heat_rate.jl" => "other.jl")))
    @test !isempty(EDGOwnershipChecks.aggregator_violations(
        aggregator * "\ninclude(joinpath(@__DIR__, \"heat_rate.jl\"))"))
    missing=copy(sources);delete!(missing,EDGOwnershipChecks.OWNER*"targeting.jl")
    @test !isempty(EDGOwnershipChecks.violations(missing))
end

@testset "Actuator errors preserve preceding decision and partial application" begin
    x=_edg_test_context(max_energy_submodes=(:heat_rate,),heat_rate_limit_w_cm2=Inf)
    # The old actuator updates link 2 before rejecting 99; decision state is
    # already committed. Extraction must not prevalidate, reorder or roll back.
    effector=S.SolarPanelAngleOfAttackControlModel((2,99))
    model=S.AerobrakingEnergyDepletionControlModel(x.config,x.state;aoa_effector=effector)
    x.state.selected_mode[1]=:safe_low_drag
    old=deepcopy((x=x,model=model))
    @test_throws ArgumentError S.calcControlEffect!(model,x.u,x.p,1.,1)
    @test_throws ArgumentError reference_control!(old.model,old.x.u,old.x.p,1.,1)
    @test isequal(state_values(x.state),state_values(old.x.state))
    @test geometry(x.spacecraft)==geometry(old.x.spacecraft)
    @test x.spacecraft.links[2].α==x.config.min_alpha_rad
    @test x.state.last_alpha_rad[1]==x.config.min_alpha_rad
end

end # module EDGOwnershipTests
