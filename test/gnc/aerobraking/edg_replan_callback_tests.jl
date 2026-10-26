# SR-1 uses the assembled production callbacks and physical planner. The fixed
# orbit and prescribed heat ramp isolate scheduling, not trajectory accuracy.
struct EDGReplanCallbackProbe <: SimulationModel.AbstractTypes.AbstractControlEffectorModel
    state::SimulationModel.AerobrakingEnergyDepletionState
    records::Vector{NamedTuple}
end
function SimulationModel.update_control!(probe::EDGReplanCallbackProbe,u,p,t::Float64,i::Int64)
    st=probe.state
    push!(probe.records,(time=t,count=st.energy_bracketing_count[i],
        last=st.last_switch_solve_t[i],switch=st.targeting_switch_s[i],
        diagnostic=st.targeting_diagnostics[i],alpha=p.args.dynamics_model.spacecraft[i].links[2].α,
        decision_alpha=st.last_alpha_rad[i]))
    return nothing
end

@testset "EDG assembled callbacks have one targeting replan owner" begin
    SM=SimulationModel; CB=SM.SimulationCallbacks
    c=_edg_test_context(guidance_modes=(:targeting,:max_energy_depletion),
        density_reference_kg_m3=1e-7,heat_load_limit_j_cm2=30.0,
        gravity_models=(SM.InverseSquaredJ2GravityModel(),))
    c.args.control_model.control_rates[1]=1/3
    c.state.control_clock_period_s=1/3   # as the solve reports it at its start
    c.args.guidance_model.guidance_rates[1]=1/3
    planet=c.args.environment_model.planet
    rp=planet.Rp_e+100e3;ra=_ODYSSEY_RA_INITIAL_M;a=(rp+ra)/2
    e=(ra-rp)/(ra+rp);sl=a*(1-e^2);radius=planet.Rp_e+150e3
    nu=acos((sl/radius-1)/e)
    r=SVector(radius*cos(nu),radius*sin(nu),0.0)
    v=sqrt(planet.μ/sl)*SVector(-sin(nu),e+cos(nu),0.0)
    # The reachable target for a switch at 10 s, as SpaceAGORAEDG's kernels give it
    # (recorded on v0.6.0; SpaceAGORAEDG's tests check the kernels still do).
    target=9.750451931596396e6
    fields=(;(k=>getfield(c.config,k) for k in fieldnames(typeof(c.config)))...)
    config=SM.AerobrakingEnergyDepletionConfig(;merge(fields,
        (target_apoapsis_radius_m=target,switch_recompute_interval_s=2.0))...)
    guide=SM.AerobrakingEnergyDepletionGuidanceModel(config,c.state)
    control=SM.AerobrakingEnergyDepletionControlModel(config,c.state)
    records=NamedTuple[];probe=EDGReplanCallbackProbe(c.state,records)
    afields=(;(k=>getfield(c.args,k) for k in fieldnames(typeof(c.args)))...)
    args=SM.SimulationConfiguration(;merge(afields,(
        simulation_settings=SM.SimulationSettings(results=false,verbose=false,generate_plots=false,save_csv=false),
        guidance_model=SM.GuidanceModel((guide,),[1/3]),
        control_model=SM.ControlModel((control,probe),[1/3,1/3])))...)
    params=SM.ODEParams(n_sats=1, args=args)
    params.shared_buffers.et_start[]=c.p.shared_buffers.et_start[]
    u=deepcopy(c.u);u.sc[1].pos .= r;u.sc[1].vel .= v
    u.sc[1].heat_loads .= [100.0,8.0,9.0]
    c.state.heat_load_entry_j_cm2[1]=[100.0,6.0,6.0]
    SM.update_guidance!(guide,u,params,0.0,1)
    initial=c.state.targeting_diagnostics[1]
    @test initial.mode==:targeting
    @test initial.certificate.heat_feasible
    @test initial.certificate.heat_load_j_cm2>=3.0
    @test c.state.targeting_switch_s[1]==10.0
    callbacks=CB.get_callbacks(1,args.dynamics_model.dynamic_effectors,args)
    function fixture_rhs!(du,u,p,t)
        fill!(du,0.0)
        du.sc[1].heat_loads[2]=0.05
        du.sc[1].heat_loads[3]=0.05
    end
    prob=CB.ODEProblem(fixture_rhs!,u,(0.0,4.1),params)
    sol=CB.solve(prob,CB.Tsit5();callback=callbacks,dtmax=0.1,
        abstol=1e-10,reltol=1e-10,save_everystep=false)
    @test length(records)==12
    for (k,rec) in enumerate(records)
        @test rec.time ≈ k/3 atol=2e-14
        expected=k<6 ? 1 : k<12 ? 2 : 3
        @test rec.count==expected # Observed after EDG control, before ordinary guidance.
        @test rec.last ≈ (expected-1)*2.0 atol=2e-14
        @test rec.diagnostic.mode==:targeting
        @test rec.diagnostic.certificate.heat_feasible
        @test rec.diagnostic.certificate.instantaneous_admissible
        @test rec.diagnostic.certificate.heat_load_j_cm2>=3.0+0.05*rec.last
        @test rec.alpha==rec.decision_alpha==config.max_alpha_rad # The panel follows the high-drag command before the switch.
    end
    @test c.state.energy_bracketing_count[1]==3
    # Accumulated heat must be rechecked by guidance at the next due tick.
    # Control alone consumes the accepted plan and cannot consume that tick.
    hot=deepcopy((guide=guide,control=control,params=params,u=sol.u[end]))
    hot.u.sc[1].heat_loads .= [100.0,40.0,41.0]
    hotstate=hot.guide.state
    @test hotstate === hot.control.state
    @test hotstate.heat_load_entry_j_cm2[1]==[100.0,6.0,6.0] && isempty(hotstate.heat_load_exit_j_cm2[1])
    SM.update_control!(hot.control,hot.u,hot.params,6.0,1)
    @test hotstate.energy_bracketing_count[1]==3
    @test hotstate.last_switch_solve_t[1]==4.0
    SM.update_guidance!(hot.guide,hot.u,hot.params,6.0,1)
    @test hotstate.energy_bracketing_count[1]==4
    @test hotstate.last_switch_solve_t[1]==6.0
    @test !hotstate.targeting_active[1]
    @test hotstate.targeting_diagnostics[1].status==:no_admissible_family
    @test !hotstate.targeting_diagnostics[1].certificate.heat_feasible
    @test hotstate.targeting_diagnostics[1].certificate.heat_load_j_cm2>=35.0
    @test hot.u.sc[1].heat_loads==[100.0,40.0,41.0]
    # Standalone control consumes a guidance-owned plan without solving or
    # advancing its clock, including an otherwise due replan tick.
    snapshot=(c.state.targeting_switch_s[1],c.state.last_switch_solve_t[1],c.state.targeting_diagnostics[1])
    SM.update_control!(control,sol.u[end],params,6.0,1)
    @test (c.state.targeting_switch_s[1],c.state.last_switch_solve_t[1],c.state.targeting_diagnostics[1])==snapshot
    # After the accepted switch, the guidance cache remains fixed.
    switched=c.state.targeting_switch_s[1]+1/3
    SM.update_guidance!(guide,sol.u[end],params,switched,1)
    SM.update_control!(control,sol.u[end],params,switched,1)
    @test c.state.energy_bracketing_count[1]==3
    @test c.spacecraft.links[2].α==config.min_alpha_rad
end
