# These checks exercise production control and physical telemetry without changing limits.
if !isdefined(@__MODULE__, :REPO_ROOT)
    const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
end
if !isdefined(@__MODULE__, :_edg_test_context)
    source = read(joinpath(REPO_ROOT, "test/gnc/aerobraking/energy_depletion_gnc_tests.jl"), String)
    include_string(@__MODULE__, first(split(source, "\n@testset")), "edg_context_definitions.jl")
end
using Test
let SM=SimulationModel, CH=SimulationModel.ControlHooks
@testset "E2c accumulated heat and physical-limit diagnostics" begin
    for modes in ((:max_energy_depletion,), (:targeting,:max_energy_depletion)),
            security in (false,true), base in (0.4,pi/2), load in (29.0,30.0,34.56)
        c=_edg_test_context(guidance_modes=modes,max_energy_submodes=(:heat_load,),heat_load_limit_j_cm2=30.0)
        fields=(; (key=>getfield(c.config,key) for key in fieldnames(typeof(c.config)))...)
        config=SM.AerobrakingEnergyDepletionConfig(;merge(fields,(heat_load_security_mode=security,))...)
        model=SM.AerobrakingEnergyDepletionControlModel(config,c.state)
        env=CH._edg_environment_state(c.u,c.p,0.0,1)
        angle=CH._edg_command_alpha!(model,c.p,c.u,env,c.spacecraft,base,load,false,1)
        @test angle == (load >= 30.0 ? config.min_alpha_rad : base)
        @test c.state.last_heat_budget_status[1] == (load >= 30.0 ? :exhausted : :available)
        @test c.state.last_heat_load_j_cm2[1] == load
        @test c.state.last_heat_rate_status[1] == :disabled
        @test c.state.last_structural_load_status[1] == :disabled
    end
    for (enabled,limit,value,minimum,status) in (
            (false,0.5,1.0,0.6,:disabled), (true,Inf,1.0,0.6,:unbounded),
            (true,NaN,1.0,0.6,:invalid_limit), (true,-1.0,1.0,0.6,:invalid_limit),
            (true,0.5,NaN,0.4,:unobserved), (true,0.5,0.5,0.4,:within_limit),
            (true,0.5,0.6,0.4,:command_above_limit),
            (true,0.5,0.6,0.55,:above_limit_at_minimum_angle),
            (true,0.5,0.4,0.55,:within_limit))
        @test CH._edg_constraint_status(enabled,limit,value,minimum) == status
    end
    for (submodes,limit,load,expected) in (
            ((:heat_rate,),30.0,34.56,:disabled), ((:heat_load,),Inf,34.56,:unbounded),
            ((:heat_load,),0.0,34.56,:invalid_limit), ((:heat_load,),NaN,34.56,:invalid_limit),
            ((:heat_load,),30.0,NaN,:unobserved))
        config=SM.AerobrakingEnergyDepletionConfig(max_energy_submodes=submodes,heat_load_limit_j_cm2=limit)
        @test CH._edg_heat_budget_status(config,load)==expected
    end
    c=_edg_test_context(heat_load_limit_j_cm2=30.0)
    env=CH._edg_environment_state(c.u,c.p,0.0,1)
    for scale in (1.0,1e-6)
        scaled=merge(env,(rho=env.rho*scale,dynamic_pressure=env.dynamic_pressure*scale))
        angle=CH._edg_command_alpha!(c.control,c.p,c.u,scaled,c.spacecraft,pi/2,34.56,false,1)
        @test angle==c.config.min_alpha_rad
        @test c.state.last_heat_budget_status[1]==:exhausted
        @test c.state.last_heat_rate_w_cm2[1]==c.state.last_minimum_heat_rate_w_cm2[1]
        @test c.state.last_structural_load_pa[1]==c.state.last_minimum_structural_load_pa[1]
        @test c.state.last_structural_load_status[1]==(scale==1.0 ? :above_limit_at_minimum_angle : :within_limit)
        @test c.state.last_heat_rate_w_cm2[1]>0.0 # Minimum angle cannot erase incoming heat.
    end
    # Disabling only accumulated-heat control preserves the existing angle choice.
    c=_edg_test_context(max_energy_submodes=(:heat_rate,),heat_load_limit_j_cm2=30.0)
    env=CH._edg_environment_state(c.u,c.p,0.0,1)
    a=CH._edg_command_alpha!(c.control,c.p,c.u,env,c.spacecraft,pi/2,29.0,false,1)
    b=CH._edg_command_alpha!(c.control,c.p,c.u,env,c.spacecraft,pi/2,34.56,false,1)
    @test a==b
    @test c.state.last_heat_budget_status[1]==:disabled
end

end
