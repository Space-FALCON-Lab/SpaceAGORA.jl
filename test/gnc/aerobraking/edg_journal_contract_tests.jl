using Test
@testset "Journal targeting family and candidate mapping" begin
    # SpaceAGORAEDG's tests check the candidate family and its mapping.
    @test SimulationModel.AerobrakingEnergyDepletionConfig().switch_recompute_interval_s==2.0
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

@testset "Journal physical feasible family, J2 mapping and live replanning" begin
    SM=SimulationModel
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
    # The reachable target for a switch at 10 s, its predicted energy and the
    # predicted state 2 s into the passage, as SpaceAGORAEDG's kernels give them
    # (recorded on v0.6.0; SpaceAGORAEDG's tests check the kernels still do).
    target=9.750451931596259e6
    known_energy_jkg=-3.2213806175077297e6
    position_at_2s=SVector(3.3890335691519254e6, 1.0467160974865406e6, 0.0)
    velocity_at_2s=SVector(-850.7081980339041, 4121.059596430641, 0.0)
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
    params=SM.ODEParams(n_sats=1, args=args)
    params.shared_buffers.et_start[]=c.p.shared_buffers.et_start[]

    u=deepcopy(c.u);u.sc[1].pos .= r;u.sc[1].vel .= v
    SM.update_guidance!(guide,u,params,0.0,1)
    plan=c.state.targeting_diagnostics[1]
    @test plan.mode==:targeting
    @test plan.certificate.instantaneous_admissible
    @test plan.certificate.heat_feasible
    @test plan.candidate_count > 10
    @test abs(plan.predicted_apoapsis_error_m)<0.01
    @test plan.switch_s==10.0
    @test abs(plan.coast_energy_change_jkg)>1.0 # Nonzero J2 work exercises correction.
    @test plan.target_energy_jkg ≈ known_energy_jkg atol=0.001 rtol=0.0
    @test plan.left_residual_jkg>=0>=plan.right_residual_jkg
    SM.update_guidance!(guide,u,params,1.0,1)
    @test c.state.energy_bracketing_count[1]==1
    # Advance to the state reached under the accepted high-angle command.
    u.sc[1].pos .= position_at_2s;u.sc[1].vel .= velocity_at_2s
    SM.update_guidance!(guide,u,params,2.0,1)
    @test c.state.energy_bracketing_count[1]==2
    @test c.state.last_switch_solve_t[1]==2.0
    @test c.state.targeting_active[1]
    @test c.state.targeting_switch_s[1]==10.0
    # Passed-switch cache is independent of later target residuals.
    SM.update_guidance!(guide,u,params,11.0,1)
    @test c.state.energy_bracketing_count[1]==2
end
