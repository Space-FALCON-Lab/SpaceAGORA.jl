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
function _edg_counted_drag_callback(num_sats)
    callback = SimulationModel.SimulationCallbacks.get_drag_state_callback(num_sats)
    entries, exits = Float64[], Float64[]
    up!(integrator, idx) = (push!(exits, Float64(integrator.t)); callback.affect!(integrator, idx))
    down!(integrator, idx) = (push!(entries, Float64(integrator.t)); callback.affect_neg!(integrator, idx))
    counted = NamedTuple{(:condition, :affect!, :affect_neg!, :save_positions)}(
        (callback.condition, up!, down!, callback.save_positions))
    return counted, entries, exits
end

@testset "EDG baseline is captured at the actual entry root" begin
    for stationary in (false, true)
        CH = SimulationModel.ControlHooks
        c = _edg_test_context(heat_load_limit_j_cm2=30.0)
        fields = (; (key=>getfield(c.args,key) for key in fieldnames(typeof(c.args)))...)
        environment = c.args.environment_model
        if stationary
            planet = environment.planet
            planet_fields = (; (key=>getfield(planet,key) for key in fieldnames(typeof(planet)))...)
            stationary_planet = typeof(planet)(;merge(planet_fields, (ω=zero(planet.ω),))...)
            environment_fields = (; (key=>getfield(environment,key) for key in fieldnames(typeof(environment)))...)
            environment = SimulationModel.EnvironmentModel(;merge(environment_fields,(planet=stationary_planet,))...)
        end
        args = SimulationModel.SimulationConfiguration(;merge(fields,(
            environment_model=environment,
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
        counted_drag, drag_entries, drag_exits = _edg_counted_drag_callback(1)
        callback = SimulationModel.SimulationCallbacks.get_edg_heat_callback(1;
            drag_callback=counted_drag)
        problem = ODEProblem(analytic_passes!,copy(c.u),(0.0,40.0),c.p)
        sol = solve(problem,Tsit5();callback=callback,dtmax=0.1,reltol=1e-10,abstol=1e-10)
        @test length(drag_entries) == 2
        @test length(drag_exits) == 2
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
        @test !SimulationModel.SimulationCallbacks._requires_drag_state_callback((),equal_args)
        callbacks = SimulationModel.SimulationCallbacks.get_callbacks(1, (), equal_args)
        heat_type = typeof(SimulationModel.SimulationCallbacks.get_edg_heat_callback(1))
        @test any(cb -> cb isa heat_type, callbacks.continuous_callbacks)
    end
end

@testset "EDG heat uses its geodetic passage at high latitude" begin
    CH = SimulationModel.ControlHooks
    CB = SimulationModel.SimulationCallbacks
    for latitude in (60.0, 90.0, -90.0)
        c = _edg_test_context(heat_load_limit_j_cm2=30.0)
        planet = c.args.environment_model.planet
        direction = SVector(cosd(latitude), 0.0, sind(latitude))
        # Solve only the static geometric equation to give the analytic radial
        # oscillator a known EDG entry radius. No atmospheric model is fitted.
        lower, upper = planet.Rp_p + 150000.0, planet.Rp_e + 170000.0
        for _ in 1:60
            mid = (lower + upper) / 2
            if CH.rtolatlong(mid * direction, planet)[1] < 160000.0
                lower = mid
            else
                upper = mid
            end
        end
        entry_radius = (lower + upper) / 2
        fields = (; (key=>getfield(c.args,key) for key in fieldnames(typeof(c.args)))...)
        args = SimulationModel.SimulationConfiguration(;merge(fields, (
            integration_tolerances=SimulationModel.IntegrationTolerances(
                dt_max_atmosphere=0.1, dt_max_orbit=0.1,
                reltol_atmosphere=1e-10, reltol_orbit=1e-10,
                abstol_atmosphere=1e-10, abstol_orbit=1e-10),))...)
        p = SimulationModel.ODEParams{1}(args=args)
        p.shared_buffers.et_start[] = c.p.shared_buffers.et_start[]
        u = copy(c.u)
        u.sc[1].pos .= (entry_radius + 10000.0) * direction
        u.sc[1].heat_loads .= [0.0, 100.0, 200.0]
        CH._edg_initialize_heat_accounting!(args, u, 0.0)
        # The gate and bookkeeping must agree even in the shell between the
        # old spherical callback and EDG's geodetic entry.
        @test norm(u.sc[1].pos) - planet.Rp_e < 160000.0
        @test CH._edg_heat_boundary_distance(u, p, 0.0, 1) > 0.0
        env = CH._edg_environment_state(u, p, 0.0, 1)
        @test !CH._edg_in_drag_passage(p, env)
        @test CH._edg_heat_boundary_distance(u, p, 0.0, 1) ≈ env.altitude_m - 160000.0
        # Start above both surfaces so both spherical and geodetic callbacks
        # cross on every passage, with a clear shell between their events.
        u.sc[1].pos .= (entry_radius + 30000.0) * direction
        function analytic_high_latitude!(du, u, p, t)
            fill!(du, 0.0)
            du.sc[1].pos .= (-6000pi*sin(pi*t/10)) * direction
            du.sc[1].heat_loads .= [0.0, 1.0, 2.0]
        end
        # Both callbacks run together: the old spherical events must neither
        # replace nor freeze the geodetic snapshots.
        counted_drag, drag_entries, drag_exits = _edg_counted_drag_callback(1)
        callback = CB.get_edg_heat_callback(1; drag_callback=counted_drag)
        sol = solve(ODEProblem(analytic_high_latitude!, u, (0.0,40.0), p), Tsit5();
            callback=callback, dtmax=0.1, reltol=1e-10, abstol=1e-10)
        @test length(drag_entries) == 2
        @test length(drag_exits) == 2
        entry2 = 20.0 + 10acos(0.5)/pi
        exit2 = 40.0 - 10acos(0.5)/pi
        @test string(sol.retcode) == "Success"
        @test c.state.heat_load_entry_j_cm2[1] ≈ [0.0,100+entry2,200+2entry2] atol=1e-5
        @test c.state.heat_load_exit_j_cm2[1] ≈ [0.0,100+exit2,200+2exit2] atol=1e-5
        @test CH._edg_pass_heat_load_for_links(sol.u[end].sc[1],(2,3),c.state,1) ≈ 2*(exit2-entry2) atol=1e-5
        @test sol.u[end].sc[1].heat_loads ≈ [0.0,140.0,280.0] atol=1e-8
        @test isempty(filter(x -> !isfinite(x), c.state.heat_load_entry_j_cm2[1]))
        # A second solve genuinely starts within EDG's passage; its first exit
        # freezes heat against the initial zero integral, without inventing an entry.
        inside = copy(u)
        inside.sc[1].pos .= (entry_radius - 10000.0) * direction
        fill!(inside.sc[1].heat_loads, 0.0)
        CH._edg_initialize_heat_accounting!(args, inside, 0.0)
        function analytic_exit!(du, u, p, t)
            fill!(du, 0.0)
            du.sc[1].pos .= 1000.0 * direction
            du.sc[1].heat_loads .= [0.0, 1.0, 2.0]
        end
        exit_sol = solve(ODEProblem(analytic_exit!, inside, (0.0,35.0), p), Tsit5();
            callback=callback, dtmax=0.1, reltol=1e-10, abstol=1e-10)
        @test length(drag_entries) == 2
        @test length(drag_exits) == 3
        @test c.state.heat_load_entry_j_cm2[1] == zeros(3)
        @test c.state.heat_load_exit_j_cm2[1] ≈ [0.0,10.0,20.0] atol=1e-5
        @test CH._edg_pass_heat_load_for_links(exit_sol.u[end].sc[1],(2,3),c.state,1) ≈ 20.0 atol=1e-5
        @test exit_sol.u[end].sc[1].heat_loads ≈ [0.0,35.0,70.0] atol=1e-8
    end
end
