using Test, ComponentArrays, StaticArrays
if !isdefined(@__MODULE__, :_edg_test_context)
    source = read(joinpath(REPO_ROOT, "test/gnc/aerobraking/energy_depletion_gnc_tests.jl"), String)
    include_string(@__MODULE__, first(split(source, "\n@testset")), "edg_context_definitions.jl")
end

@testset "EDG passage heat bookkeeping preserves the control law" begin
    CH = SimulationModel.ControlHooks
    @testset "Panel identity, exact baselines, and fresh-run reset" begin
        c = _edg_test_context(heat_load_limit_j_cm2=30.0)
        CH._edg_initialize_heat_accounting!(c.args, c.u, 0.0)
        @test length(CH._edg_heat_states(c.args)) == 1 # Shared control/guidance state.
        c.u.sc[1].heat_loads .= [999.0, 100.0, 20.0]
        original = copy(c.u)
        CH._edg_capture_entry_heat!(c.args, c.u, 1)
        @test c.u == original
        c.u.sc[1].heat_loads .+= [100.0, 5.0, 8.0]
        @test CH._edg_max_heat_load_for_links(c.u.sc[1], (2,3)) == 105.0
        @test CH._edg_pass_heat_load_for_links(c.u.sc[1], (2,3), c.state, 1) == 8.0
        @test CH._edg_pass_heat_load_for_links(c.u.sc[1], (2,), c.state, 1) == 5.0
        @test CH._edg_pass_heat_load_for_links(c.u.sc[1], (3,), c.state, 1) == 8.0
        CH._edg_capture_exit_heat!(c.args, c.u, 1)
        c.u.sc[1].heat_loads .+= [0.0,50.0,50.0] # Coast telemetry stays cumulative.
        @test CH._edg_pass_heat_load_for_links(c.u.sc[1], (2,3), c.state, 1) == 8.0
        CH._edg_capture_entry_heat!(c.args, c.u, 1)
        @test CH._edg_pass_heat_load_for_links(c.u.sc[1], (2,3), c.state, 1) == 0.0
        # A fresh simulation has zero integrated heat and must discard old baselines.
        fill!(c.u.sc[1].heat_loads, 0.0)
        CH._edg_initialize_heat_accounting!(c.args, c.u, 0.0)
        c.u.sc[1].heat_loads .= [0.0, 2.0, 3.0]
        @test CH._edg_pass_heat_load_for_links(c.u.sc[1], (2,3), c.state, 1) == 3.0
        @test_throws ArgumentError CH._edg_initialize_heat_accounting!(c.args,c.u,5000.0)
    end
    @testset "Spacecraft accounting is independent" begin
        state = SimulationModel.AerobrakingEnergyDepletionState(num_sats=2)
        config = SimulationModel.AerobrakingEnergyDepletionConfig()
        model = SimulationModel.AerobrakingEnergyDepletionControlModel(config,state)
        args = (control_model=(control_effectors=(model,),),guidance_model=(guidance_effectors=(),))
        u = ComponentVector(sc=[(heat_loads=[0.0,26.0,20.0],),(heat_loads=[0.0,80.0,90.0],)])
        CH._edg_initialize_heat_accounting!(args,u,0.0)
        CH._edg_capture_entry_heat!(args,u,1)
        @test state.heat_load_entry_j_cm2[2] == zeros(3)
        CH._edg_capture_entry_heat!(args,u,2)
        u.sc[1].heat_loads .+= [0,4,8]
        u.sc[2].heat_loads .+= [0,9,2]
        @test CH._edg_pass_heat_load_for_links(u.sc[1],(2,3),state,1)==8.0
        @test CH._edg_pass_heat_load_for_links(u.sc[2],(2,3),state,2)==9.0
        CH._edg_capture_entry_heat!(args,u,1)
        @test CH._edg_pass_heat_load_for_links(u.sc[2],(2,3),state,2)==9.0
        unrelated=(control_model=(control_effectors=(nothing,),),guidance_model=(guidance_effectors=(),))
        @test isempty(CH._edg_heat_states(unrelated))
        @test isnothing(CH._edg_initialize_heat_accounting!(unrelated,u,5000.0))
    end
    @testset "Real controller and planner ignore prior-passage offsets" begin
        for modes in ((:max_energy_depletion,), (:targeting,:max_energy_depletion))
            zero = _edg_test_context(guidance_modes=modes, heat_load_limit_j_cm2=30.0)
            offset = _edg_test_context(guidance_modes=modes, heat_load_limit_j_cm2=30.0)
            for (c, prior) in ((zero, [0.0,0.0,0.0]), (offset,[1000.0,26.0,40.0]))
                CH._edg_initialize_heat_accounting!(c.args,c.u,0.0)
                c.u.sc[1].heat_loads .= prior
                CH._edg_capture_entry_heat!(c.args,c.u,1)
                c.u.sc[1].heat_loads .+= [0.0,8.0,5.0]
                SimulationModel.calcGuidanceEffect!(c.guidance,c.u,c.p,0.0,Int64(1))
                SimulationModel.calcControlEffect!(c.control,c.u,c.p,0.0,Int64(1))
            end
            for field in (:selected_mode,:targeting_active,:target_energy_jkg,
                    :bracket_min_energy_jkg,:bracket_max_energy_jkg,:targeting_switch_s,
                    :heat_load_switches_s,:last_alpha_rad,:last_heat_rate_w_cm2,
                    :last_structural_load_pa,:last_pass_heat_load_j_cm2,:last_heat_budget_status)
                @test isequal(getfield(zero.state,field),getfield(offset.state,field))
            end
            @test offset.state.last_heat_load_j_cm2[1] == 45.0
            @test offset.state.last_pass_heat_load_j_cm2[1] == 8.0
            @test offset.state.last_heat_budget_status[1] == :available
            @test offset.u.sc[1].heat_loads == [1000.0,34.0,45.0]
        end
    end
end

# Exercise the existing continuous entry event, with an analytic trajectory and
# heat integral. This tests accounting, not a new aerobraking campaign.
if !isdefined(@__MODULE__, :SimulationEngine)
    include(joinpath(REPO_ROOT, "src/simulation/engine/simulation_engine.jl"))
end
using DifferentialEquations
@testset "EDG baseline is captured at the actual entry root" begin
    CH = SimulationModel.ControlHooks
    c = _edg_test_context(heat_load_limit_j_cm2=30.0)
    fields = (; (key=>getfield(c.args,key) for key in fieldnames(typeof(c.args)))...)
    args = SimulationModel.SimulationConfiguration(;merge(fields,(
        integration_tolerances=SimulationModel.IntegrationTolerances(
            dt_max_atmosphere=0.1,dt_max_orbit=0.1,reltol_atmosphere=1e-12,reltol_orbit=1e-12,
            abstol_atmosphere=1e-12,abstol_orbit=1e-12),))...)
    c=merge(c,(args=args,p=SimulationModel.ODEParams{1}(args=args)))
    radius = c.p.args.environment_model.planet.Rp_e
    c.u.sc[1].pos .= [radius+170000.0,0.0,0.0]
    c.u.sc[1].heat_loads .= [0.0,100.0,200.0]
    CH._edg_initialize_heat_accounting!(c.args,c.u,0.0)
    function analytic_passes!(du,u,p,t)
        fill!(du,0.0)
        du.sc[1].pos[1] = -4000pi*sin(pi*t/10)
        du.sc[1].heat_loads .= [0.0,1.0,2.0]
    end
    callback = SimulationModel.SimulationCallbacks.get_drag_state_callback(1)
    problem = ODEProblem(analytic_passes!,copy(c.u),(0.0,40.0),c.p)
    sol = solve(problem,Tsit5();callback=callback,dtmax=0.1,reltol=1e-10,abstol=1e-10)
    entry2 = 20.0 + 10acos(0.75)/pi
    @test string(sol.retcode)=="Success"
    @test c.state.heat_load_entry_j_cm2[1] ≈ [0.0,100+entry2,200+2entry2] atol=1e-5
    @test sol.u[end].sc[1].heat_loads ≈ [0.0,140.0,280.0] atol=1e-8
    exit2 = 40.0 - 10acos(0.75)/pi
    @test c.state.heat_load_exit_j_cm2[1] ≈ [0.0,100+exit2,200+2exit2] atol=1e-5
    @test CH._edg_pass_heat_load_for_links(sol.u[end].sc[1],(2,3),c.state,1) ≈ 2*(exit2-entry2) atol=1e-5
    @test !c.p.shared_buffers.in_atmosphere[1]
    # EDG still receives the entry root when integration tolerances match.
    equal_args=SimulationModel.SimulationConfiguration(;merge(fields,(
        integration_tolerances=SimulationModel.IntegrationTolerances(
            dt_max_atmosphere=1.0,dt_max_orbit=1.0,reltol_atmosphere=1e-8,reltol_orbit=1e-8,
            abstol_atmosphere=1e-8,abstol_orbit=1e-8),))...)
    @test SimulationModel.SimulationCallbacks._requires_drag_state_callback((),equal_args)
end
