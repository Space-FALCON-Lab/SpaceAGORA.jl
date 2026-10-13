using Test
@testset "Journal targeting family and candidate mapping" begin
    CH=SimulationModel.ControlHooks
    candidate(k; admissible=true, heat=true, energy=-10.0-k) =
        (switch_s=Float64(k),energy_jkg=energy,
         certificate=(instantaneous_admissible=admissible,heat_feasible=heat))
    family=[candidate(k) for k in 0:6]
    mapped(target)=c->(residual_jkg=c.energy_jkg-target,)
    for (target,mode,status) in ((-9.0,:safe_low_drag,:above_upper_energy),
            (-17.0,:max_energy_depletion,:below_lower_energy),
            (-12.5,:targeting,:quantized_bracket),(-13.0,:targeting,:exact_sample),
            (-10.0,:targeting,:exact_sample),(-16.0,:targeting,:exact_sample))
        x=CH._edg_select_mapped_candidate(family,mapped(target))
        @test x.mode==mode
        @test x.status==status
        if mode==:targeting
            @test abs(x.mapping.residual_jkg)<=0.5
            @test x.left_residual_jkg>=0>=x.right_residual_jkg
        end
    end
    @test_throws ArgumentError CH._edg_select_mapped_candidate(NamedTuple[],mapped(-11.0))
    # A hole must stop the connected prefix, never be skipped to reach a target.
    for (bad,expected) in ((:heat,:heat_boundary),(:instantaneous,:instantaneous_infeasibility))
        seen=Int[]
        evaluate(k)=(push!(seen,k);candidate(k;heat=!(bad==:heat && k==3),admissible=!(bad==:instantaneous && k==3)))
        x=CH._edg_ordered_feasible_prefix(0:6,evaluate,-100.0)
        @test length(x.accepted)==3
        @test x.status==expected
        @test seen==[0,1,2,3]
    end
    x=CH._edg_ordered_feasible_prefix(0:6,k->candidate(k;energy=k==3 ? -11.0 : -10.0-k),-100.0)
    @test x.status==:ordering_not_resolved
    @test length(x.accepted)==3
    x=CH._edg_ordered_feasible_prefix(0:6,k->candidate(k),-12.5)
    @test x.status==:depletion_reference_boundary
    @test last(x.accepted).energy_jkg==-12.0
    @test x.rejected.energy_jkg==-13.0
    x=CH._edg_ordered_feasible_prefix(0:6,k->candidate(k),-13.0)
    @test x.status==:depletion_reference_reached
    @test length(x.accepted)==4
    x=CH._edg_ordered_feasible_prefix(0:6,k->candidate(k;admissible=false),-13.0)
    @test isempty(x.accepted)
    @test x.status==:instantaneous_infeasibility
    # Candidate-specific mappings can expose a nonmonotone outer residual.
    x=CH._edg_select_mapped_candidate(family,c->(residual_jkg=c.switch_s==3 ? 100.0 : 3.5-c.switch_s,))
    @test x.status==:mapping_ordering_not_resolved
    @test x.mode==:safe_low_drag
    @test SimulationModel.AerobrakingEnergyDepletionConfig().switch_recompute_interval_s==2.0
    @test CH._edg_replan_due(0.0,-Inf,2.0)
    for start in (0.0,10.0,123.0)
        @test all(!CH._edg_replan_due(start+k/3,start,2.0) for k in 0:5)
        @test CH._edg_replan_due(start+2,start,2.0)
    end
end

@testset "Journal new-passage planning reset" begin
    c=_edg_test_context(guidance_modes=(:targeting,))
    c.state.energy_bracketing_evaluated[1]=true
    c.state.targeting_active[1]=true
    c.state.targeting_switch_s[1]=123.0
    c.state.last_switch_solve_t[1]=100.0
    c.state.targeting_diagnostics[1]=(status=:quantized_bracket,)
    SimulationModel.ControlHooks._edg_capture_entry_heat!(c.args,c.u,1)
    @test !c.state.energy_bracketing_evaluated[1]
    @test !c.state.targeting_active[1]
    @test c.state.targeting_switch_s[1]==Inf
    @test c.state.last_switch_solve_t[1]==-Inf
    @test c.state.targeting_diagnostics[1].status==:unplanned
end

module EDGJournalSecurityFixture
using StaticArrays
struct ODEParams end
const RATE_SCALE=Ref(1.0)
const CALLS=Ref(0)
_edg_predict_mass(sc,mass)=mass
_edg_remaining_closed_form_heat_load_trajectory(args...)=(CALLS[]+=1;(time=[0.0,1.0,2.0],))
_edg_constrained_heat_load_alpha_profile(config,p,sc,links,track,profile;kwargs...)=profile
_edg_profile_heat_rates(config,p,track,profile;kwargs...)=RATE_SCALE[].*profile
_edg_in_drag_passage(p,env)=true
end
for name in (:_edg_profile_heat_load,:_edg_heat_load_security_required,:_edg_update_heat_load_security!)
    EDGE2cRepairTests.load_definitions!(EDGJournalSecurityFixture,"src/gnc/control/heat_load_control.jl",[name])
end
@testset "Journal planned-profile heat safeguard and release" begin
    F=EDGJournalSecurityFixture
    config=(min_alpha_rad=0.0,max_alpha_rad=2.0,controlled_panel_links=(2,3),
        max_energy_submodes=(:heat_load,),heat_load_limit_j_cm2=4.0,heat_load_security_mode=true)
    p=F.ODEParams();r=SVector(1.0,0.0,0.0);v=SVector(0.0,1.0,0.0)
    for t in (0.0,5000.0)
        F.RATE_SCALE[]=1.0
        @test first(F._edg_heat_load_security_required(config,p,nothing,r,v,1.0,nothing,0.0,t))
        @test !first(F._edg_heat_load_security_required(config,p,nothing,r,v,1.0,nothing,0.0,t;planned_switches=(t,t+2)))
        @test !first(F._edg_heat_load_security_required(config,p,nothing,r,v,1.0,nothing,2.0,t;planned_switches=(t+1,t+2)))
        st=(selected_mode=[:max_energy_depletion],heat_load_security_active=[false],heat_load_switches_s=[(Inf,Inf)])
        model=(;config,state=st);F.CALLS[]=0
        F._edg_update_heat_load_security!(model,p,nothing,nothing,r,v,1.0,0.0,t,1)
        @test st.heat_load_security_active[1]
        @test st.heat_load_switches_s[1]==(Inf,Inf)
        F.RATE_SCALE[]=0.5
        F._edg_update_heat_load_security!(model,p,nothing,nothing,r,v,1.0,0.0,t+1/3,1)
        @test !st.heat_load_security_active[1]
        @test F.CALLS[]==4 # Both updates evaluate the nominal profile twice.
        @test st.heat_load_switches_s[1]==(Inf,Inf)
        st.selected_mode[1]=:targeting
        F._edg_update_heat_load_security!(model,p,nothing,nothing,r,v,1.0,100.0,t+2/3,1)
        @test !st.heat_load_security_active[1]
        @test F.CALLS[]==4 # Depletion-mode safeguard does not rewrite EDG-T.
    end
end

@testset "Journal physical feasible family, J2 mapping and live replanning" begin
    SM=SimulationModel; CH=SM.ControlHooks
    # A separate low-density outbound regression fixture, not a relaxation of the
    # frozen Mars comparison. All three finite physical limits retain their values.
    c=_edg_test_context(guidance_modes=(:targeting,:max_energy_depletion),
        density_reference_kg_m3=1e-7,heat_load_limit_j_cm2=30.0,
        gravity_models=(SM.InverseSquaredJ2GravityModel(),))
    planet=c.args.environment_model.planet
    rp=planet.Rp_e+100e3; ra=_ODYSSEY_RA_INITIAL_M; a=(rp+ra)/2
    e=(ra-rp)/(ra+rp); sl=a*(1-e^2); radius=planet.Rp_e+150e3
    nu=acos((sl/radius-1)/e)
    r=SVector(radius*cos(nu),radius*sin(nu),0.0)
    v=sqrt(planet.μ/sl)*SVector(-sin(nu),e+cos(nu),0.0)
    known=CH._edg_predict_targeting_outcome(c.config,c.p,c.spacecraft,r,v,461.0,0.0,10.0;
        heat_rate_control=true,structural_control=true)
    coast=CH._edg_vacuum_apoapsis_correction(c.p,known.exit_position,known.exit_velocity,
        461.0,known.duration_s)
    @test coast.event_reached
    target=norm(coast.position)
    fields=(; (k=>getfield(c.config,k) for k in fieldnames(typeof(c.config)))...)
    config=SM.AerobrakingEnergyDepletionConfig(;merge(fields,
        (target_apoapsis_radius_m=target,switch_recompute_interval_s=2.0))...)
    guide=SM.AerobrakingEnergyDepletionGuidanceModel(config,c.state)
    # The predictor deliberately requires the matching configured effector.
    afields=(; (k=>getfield(c.args,k) for k in fieldnames(typeof(c.args)))...)
    control=SM.AerobrakingEnergyDepletionControlModel(config,c.state)
    args=SM.SimulationConfiguration(;merge(afields,
        (guidance_model=SM.GuidanceModel((guide,),[1.0]),
         control_model=SM.ControlModel((control,),[1.0])))...)
    params=SM.ODEParams{1}(args=args)
    params.shared_buffers.et_start[]=c.p.shared_buffers.et_start[]

    u=deepcopy(c.u);u.sc[1].pos .= r;u.sc[1].vel .= v
    SM.calcGuidanceEffect!(guide,u,params,0.0,1)
    plan=c.state.targeting_diagnostics[1]
    @test plan.mode==:targeting
    @test plan.certificate.instantaneous_admissible
    @test plan.certificate.heat_feasible
    @test plan.candidate_count > 10
    @test abs(plan.predicted_apoapsis_error_m)<0.01
    @test plan.switch_s==10.0
    @test abs(plan.coast_energy_change_jkg)>1.0 # Nonzero J2 work exercises correction.
    @test plan.target_energy_jkg ≈ known.energy_jkg atol=0.001 rtol=0.0
    @test plan.left_residual_jkg>=0>=plan.right_residual_jkg
    SM.calcGuidanceEffect!(guide,u,params,1.0,1)
    @test c.state.energy_bracketing_count[1]==1
    # Advance to the state reached under the accepted high-angle command.
    idx=findfirst(==(2.0),known.track.time)
    @test idx!==nothing
    u.sc[1].pos .= known.track.positions[idx];u.sc[1].vel .= known.track.velocities[idx]
    SM.calcGuidanceEffect!(guide,u,params,2.0,1)
    @test c.state.energy_bracketing_count[1]==2
    @test c.state.last_switch_solve_t[1]==2.0
    @test c.state.targeting_active[1]
    @test c.state.targeting_switch_s[1]==10.0
    # Passed-switch cache is independent of later target residuals.
    SM.calcGuidanceEffect!(guide,u,params,11.0,1)
    @test c.state.energy_bracketing_count[1]==2
end
@testset "Closed-form remaining passage keeps its entry reference" begin
    CH=SimulationModel.ControlHooks
    c=_edg_test_context(heat_load_limit_j_cm2=30.0)
    planet=c.args.environment_model.planet
    interface=1000.0*c.args.environment_model.EI
    rp=planet.Rp_e+100e3; ra=20_000e3; a=(rp+ra)/2; e=(ra-rp)/(ra+rp)
    sl=a*(1-e^2); nu_exit=acos((sl/(planet.Rp_e+interface)-1)/e)
    function state_at(nu)
        radius=sl/(1+e*cos(nu))
        pos=radius*SVector(cos(nu),sin(nu),0.0)
        vel=sqrt(planet.μ/sl)*SVector(-sin(nu),e+cos(nu),0.0)
        env=(altitude_m=radius-planet.Rp_e,)
        return pos,vel,env
    end
    for fraction in (-1.0,-0.8,-0.2,-1e-8,0.0,1e-8,0.2,0.8,0.999)
        pos,vel,env=state_at(fraction*nu_exit)
        track=CH._edg_remaining_closed_form_heat_load_trajectory(c.config,c.p,c.spacecraft,pos,vel,461.0,0.0,env;prediction_step_s=1.0)
        shifted=CH._edg_remaining_closed_form_heat_load_trajectory(c.config,c.p,c.spacecraft,pos,vel,461.0,50_000.0,env;prediction_step_s=1.0)
        @test shifted==track # Absolute mission time cannot restart the passage.
        @test all(isfinite,track.h) && all(isfinite,track.speed) && all(isfinite,track.rho)
        @test first(track.h)==env.altitude_m
        @test first(track.speed)≈norm(vel) atol=1e-10
        @test first(track.gamma)≈asin(dot(pos,vel)/(norm(pos)*norm(vel))) atol=1e-12
        @test last(track.h)≈interface atol=1e-6
        @test maximum(track.h)<=interface+1e-6
        @test 0.0<last(track.time)<1000.0
        if fraction>=0.0
            @test all(diff(track.h).>=0.0)
            @test all(track.gamma.>=-1e-12)
        end
        if fraction==-1.0
            entry=CH._edg_closed_form_heat_load_trajectory(c.config,c.p,c.spacecraft,pos,vel,461.0,0.0,env;prediction_step_s=1.0)
            @test track.time≈entry.time rtol=1e-10
            @test track.h≈entry.h rtol=1e-10
            @test track.speed≈entry.speed rtol=1e-10
            @test track.gamma≈entry.gamma atol=1e-10
        end
    end
    # At periapsis, restarting the symmetric formula holds altitude constant.
    # The remaining passage must instead climb to the outgoing interface.
    pos,vel,env=state_at(0.0)
    parent=CH._edg_closed_form_heat_load_trajectory(c.config,c.p,c.spacecraft,pos,vel,461.0,0.0,env)
    repaired=CH._edg_remaining_closed_form_heat_load_trajectory(c.config,c.p,c.spacecraft,pos,vel,461.0,0.0,env)
    @test maximum(parent.h)-minimum(parent.h)==0.0
    @test last(repaired.h)-first(repaired.h)>50_000.0
    # Continuity through periapsis: no sign-dependent horizon jump or singular fit.
    tracks=map((-1e-8,0.0,1e-8)) do fraction
        pos,vel,env=state_at(fraction*nu_exit)
        CH._edg_remaining_closed_form_heat_load_trajectory(c.config,c.p,c.spacecraft,pos,vel,461.0,0.0,env)
    end
    @test abs(last(tracks[1].time)-last(tracks[3].time))<1e-4
    @test abs(last(tracks[1].speed)-last(tracks[3].speed))<1e-3
    # A too-short horizon cannot masquerade as a complete heat forecast.
    short=(planning_horizon_s=1.0,)
    @test_throws ArgumentError CH._edg_remaining_closed_form_reference(short,c.p,pos,vel,461.0,env)
end
