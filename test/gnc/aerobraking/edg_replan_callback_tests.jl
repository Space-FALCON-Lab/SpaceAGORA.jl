# SR-1 uses the assembled production callbacks and physical planner. The fixed
# orbit and prescribed heat ramp isolate scheduling, not trajectory accuracy.
struct EDGReplanCallbackProbe <: SimulationModel.AbstractTypes.AbstractControlEffectorModel
    state::SimulationModel.AerobrakingEnergyDepletionState
    records::Vector{NamedTuple}
end
function SimulationModel.calcControlEffect!(probe::EDGReplanCallbackProbe,u,p,t::Float64,i::Int64)
    st=probe.state
    push!(probe.records,(time=t,count=st.energy_bracketing_count[i],
        last=st.last_switch_solve_t[i],switch=st.targeting_switch_s[i],
        diagnostic=st.targeting_diagnostics[i],alpha=p.args.dynamics_model.spacecraft[i].links[2].α))
    return nothing
end

@testset "EDG assembled callbacks have one targeting replan owner" begin
    SM=SimulationModel; CH=SM.ControlHooks; CB=SM.SimulationCallbacks
    c=_edg_test_context(guidance_modes=(:targeting,:max_energy_depletion),
        density_reference_kg_m3=1e-7,heat_load_limit_j_cm2=30.0,
        gravity_models=(SM.InverseSquaredJ2GravityModel(),))
    c.args.control_model.control_rates[1]=1/3
    c.args.guidance_model.guidance_rates[1]=1/3
    planet=c.args.environment_model.planet
    rp=planet.Rp_e+100e3;ra=_ODYSSEY_RA_INITIAL_M;a=(rp+ra)/2
    e=(ra-rp)/(ra+rp);sl=a*(1-e^2);radius=planet.Rp_e+150e3
    nu=acos((sl/radius-1)/e)
    r=SVector(radius*cos(nu),radius*sin(nu),0.0)
    v=sqrt(planet.μ/sl)*SVector(-sin(nu),e+cos(nu),0.0)
    known=CH._edg_predict_targeting_outcome(c.config,c.p,c.spacecraft,r,v,461.0,0.0,10.0;
        heat_rate_control=true,structural_control=true)
    coast=CH._edg_vacuum_apoapsis_correction(c.p,known.exit_position,known.exit_velocity,461.0,known.duration_s)
    @test coast.event_reached
    fields=(;(k=>getfield(c.config,k) for k in fieldnames(typeof(c.config)))...)
    config=SM.AerobrakingEnergyDepletionConfig(;merge(fields,
        (target_apoapsis_radius_m=norm(coast.position),switch_recompute_interval_s=2.0))...)
    guide=SM.AerobrakingEnergyDepletionGuidanceModel(config,c.state)
    control=SM.AerobrakingEnergyDepletionControlModel(config,c.state)
    records=NamedTuple[];probe=EDGReplanCallbackProbe(c.state,records)
    afields=(;(k=>getfield(c.args,k) for k in fieldnames(typeof(c.args)))...)
    args=SM.SimulationConfiguration(;merge(afields,(
        simulation_settings=SM.SimulationSettings(results=false,verbose=false,generate_plots=false,save_csv=false),
        guidance_model=SM.GuidanceModel((guide,),[1/3]),
        control_model=SM.ControlModel((control,probe),[1/3,1/3])))...)
    params=SM.ODEParams{1}(args=args)
    params.shared_buffers.et_start[]=c.p.shared_buffers.et_start[]
    u=deepcopy(c.u);u.sc[1].pos .= r;u.sc[1].vel .= v
    u.sc[1].heat_loads .= [100.0,8.0,9.0]
    c.state.heat_load_entry_j_cm2[1]=[100.0,6.0,6.0]
    @test CH._edg_pass_heat_load_for_links(u.sc[1],(2,3),c.state,1)==3.0
    SM.calcGuidanceEffect!(guide,u,params,0.0,1)
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
        @test rec.alpha==CH._edg_base_alpha(control,rec.time,1)
    end
    @test c.state.energy_bracketing_count[1]==3
    # Standalone control consumes a guidance-owned plan without solving or
    # advancing its clock, including an otherwise due replan tick.
    snapshot=(c.state.targeting_switch_s[1],c.state.last_switch_solve_t[1],c.state.targeting_diagnostics[1])
    SM.calcControlEffect!(control,sol.u[end],params,6.0,1)
    @test (c.state.targeting_switch_s[1],c.state.last_switch_solve_t[1],c.state.targeting_diagnostics[1])==snapshot
    # After the accepted switch, the guidance cache remains fixed.
    switched=c.state.targeting_switch_s[1]+1/3
    SM.calcGuidanceEffect!(guide,sol.u[end],params,switched,1)
    SM.calcControlEffect!(control,sol.u[end],params,switched,1)
    @test c.state.energy_bracketing_count[1]==3
    @test c.spacecraft.links[2].α==config.min_alpha_rad
end
