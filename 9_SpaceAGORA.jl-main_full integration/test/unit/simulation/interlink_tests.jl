module InterlinkTests

using Test
using ComponentArrays
using OrdinaryDiffEq
using Random
using StaticArrays
using PythonCall
include(joinpath(@__DIR__, "..", "..", "..", "examples", "gve_sma_interlinks.jl"))
using SpaceAGORA.SimulationModel.SchedulingPolicies: matching_workspace!

const SE = SpaceAGORA.SimulationEngine
const EP = SM.ExternalPropagation

edge_key(graph, edge) = map(v -> graph.vertices[v], graph.edges[edge])
semimajor_axis(state, mu) = inv(2 / norm(state.pos) - dot(state.vel, state.vel) / mu)
eccentricity(state, mu) = norm(cross(state.vel, cross(state.pos, state.vel)) / mu - state.pos / norm(state.pos))

function endpoint_fixture(; range=200e3)
    spacecraft = [SpacecraftModel(id=100 + index) for index in 1:2]
    state = ComponentVector(sc=[
        (pos=[7e6, 1000.0, 0.0], vel=[0.0, 7500.0, 0.0], mass=100.0,
            laser_dv=zeros(3), laser_delta_sma=0.0),
        (pos=[7e6, 0.0, 0.0], vel=[0.0, 7500.0, 0.0], mass=200.0,
            laser_dv=zeros(3), laser_delta_sma=0.0)])
    model = ConstellationInterLinkModel()
    key = register_allowed_link!(model, spacecraft, (1, 1), (2, 1);
        parameters=InterLinkParameters(; range))
    runtime = (args=(interlink_model=model, scheduling_policy_model=SchedulingPolicyModel(:gve_sma),
        dynamics_model=(spacecraft=spacecraft,), environment_model=(planet=make_no_gram_planet(:earth),)),
        is_active=trues(2))
    return spacecraft, state, model, key, runtime
end

function continuous_laser_fixture!(derivative, state, parameters, time)
    derivative .= 0.0
    for satellite in eachindex(state.sc)
        derivative.sc[satellite].pos .= state.sc[satellite].vel
    end
    SE._apply_interlink_rhs!(derivative, state, parameters)
    return nothing
end

function with_settings(args; overrides...)
    settings = args.simulation_settings
    names = fieldnames(typeof(settings))
    fields = NamedTuple{names}(map(name -> getfield(settings, name), names))
    return SM.SimConfig._with_configuration(args;
        simulation_settings=SimulationSettings(; merge(fields, (; overrides...))...))
end

function exhaustive_matching_score(model, policy)
    best = zero(Rational{BigInt})
    scores = matching_workspace!(model, policy).score
    for mask in 0:((1 << length(scores)) - 1)
        used = Set{Tuple{Int, Int}}()
        total = zero(Rational{BigInt})
        valid = true
        for edge in eachindex(scores)
            mask & (1 << (edge - 1)) == 0 && continue
            key = edge_key(model.graph, edge)
            if !model.graph.available[edge] || scores[edge] <= 0 || any(in(used), key)
                valid = false
                break
            end
            union!(used, key)
            total += Rational{BigInt}(scores[edge])
        end
        valid && (best = max(best, total))
    end
    return best
end

struct InterlinkExternalOwner <: EP.AbstractExternalPropagator end
EP.external_spacecraft(::InterlinkExternalOwner) = [2]

withenv("SPACEAGORA_RHS_CALIBRATE" => "off", "SPACEAGORA_RHS_IDENTIFY" => "0",
    "SPACEAGORA_RHS_EXECUTION_MODE" => "serial", "SPACEAGORA_VISUALIZATION" => "0") do
@testset "Constellation graph and laser interlink integration" begin
    @testset "Public API and configuration compatibility" begin
        @test SpaceAGORA.InterLinkParameters === SM.InterLinkParameters
        @test SpaceAGORA.schedule_interlinks! === SM.schedule_interlinks!
        spacecraft = SpacecraftModel(n_terminal=2, battery_energy_index=75, tempurature_index=40)
        copied = SE._without_gnc(spacecraft)
        @test copied.n_terminal == 2
        @test copied.battery_energy_index == 75.0
        @test copied.tempurature_index == 40.0
        @test copied.attachments === spacecraft.attachments
        old_fields = fieldnames(SpacecraftModel)[1:15]
        positional = SpacecraftModel((getfield(spacecraft, field) for field in old_fields)...)
        @test positional.n_terminal == 1
        @test positional.battery_energy_index == positional.tempurature_index == 100.0
        @test SpacecraftModel(n_terminal=0).n_terminal == 0
        @test_throws ArgumentError SpacecraftModel(n_terminal=-1)
        @test_throws ArgumentError SpacecraftModel(battery_energy_index=101)
        @test_throws ArgumentError SpacecraftModel(tempurature_index=NaN)

        args = gve_sma_configuration(duration=20.0)
        cloned = SM.SimConfig._with_configuration(args; solver_config=SM.SolverConfig(solver_mode=:dp8))
        @test cloned.interlink_model === args.interlink_model
        @test cloned.scheduling_policy_model === args.scheduling_policy_model
        legacy = SimulationConfiguration((getfield(args, field) for field in fieldnames(typeof(args))[1:12])...)
        @test legacy.interlink_model === nothing
        @test !hasproperty(SE.build_initial_conditions(legacy).sc[1], :laser_dv)
        @test legacy.scheduling_policy_model !== SimulationConfiguration(
            environment_model=args.environment_model, dynamics_model=args.dynamics_model,
            initial_time=args.initial_time).scheduling_policy_model
        @test SE._config_has_cross_satellite_coupling(args)
        @test !SE._config_has_cross_satellite_coupling(legacy)
    end

    @testset "Indexed registration and availability" begin
        spacecraft = [SpacecraftModel(id=100 + i, n_terminal=2) for i in 1:4]
        model = ConstellationInterLinkModel(forbidden_pairs=[(4, 2)])
        graph = model.graph
        first_key = @inferred register_allowed_link!(model, spacecraft, (4, 2), (2, 1);
            parameters=InterLinkParameters(P=12.0, B=3.0, range=10.0))
        second_key = register_allowed_link!(model, spacecraft, (1, 2), (2, 1);
            parameters=InterLinkParameters(P=24.0, B=5.0, range=20.0))
        @test graph.vertices == [(2, 1), (4, 2), (1, 2)]
        @test graph.edges == [(1, 2), (3, 1)]
        @test graph.vertex_neighbors == [[2, 1], [1], [2]]
        @test graph.satellite_edges == [[2], [1, 2], Int[], [1]]
        @test graph.edge_index_lookup == Dict(first_key => 1, second_key => 2)
        @test graph.available == graph.active == [false, false]
        @test [p.P for p in graph.edge_properties] == [12.0, 24.0]
        @test !hasproperty(graph, :score)
        snapshot = deepcopy(graph)
        @test_throws ArgumentError register_allowed_link!(model, spacecraft, reverse(first_key)...)
        @test_throws ArgumentError register_allowed_link!(model, spacecraft, (1, 3), (4, 1))
        @test_throws ArgumentError register_allowed_link!(model, spacecraft, (1, 1), (1, 2))
        @test_throws ArgumentError register_allowed_link!(model, spacecraft, (0, 1), (4, 1))
        @test all(getfield(graph, field) == getfield(snapshot, field) for field in fieldnames(InterLinkGraph))
        for parameters in ((P=-1,), (B=0,), (range=0,), (P=Inf,), (range=NaN,))
            @test_throws ArgumentError InterLinkParameters(; parameters...)
        end
        @test_throws ArgumentError ConstellationInterLinkModel(constellation_link_type=:invalid)
        state = (sc=[(pos=[Float64(i), 0.0, 0.0],) for i in 1:4],)
        @test update_availability!(model, spacecraft, state) == [false, true]
        available = update_availability!(model, spacecraft, state)
        available[2] = false
        @test graph.available[2]
        empty!(model.forbidden_pairs)
        @test update_availability!(model, spacecraft, state) == [true, true]
        @test update_availability!(model, spacecraft, state; is_active=[true, true, true, false]) == [false, true]
        for (position, expected) in ((2.0, false), (12.0, true), (nextfloat(12.0), false))
            state.sc[4].pos[1] = position
            @test update_availability!(model, spacecraft, state)[1] == expected
        end
        spacecraft[1].battery_energy_index = 0.0
        spacecraft[2].tempurature_index = 0.0
        @test update_availability!(model, spacecraft, state)[2]
        @test isempty(update_availability!(ConstellationInterLinkModel(), spacecraft, state))
    end

    @testset "Exact tree and cyclic terminal matching" begin
        spacecraft = [SpacecraftModel(n_terminal=2) for _ in 1:5]
        model = ConstellationInterLinkModel()
        policy = SchedulingPolicyModel(:gve_sma)
        central = register_allowed_link!(model, spacecraft, (1, 1), (2, 1))
        left = register_allowed_link!(model, spacecraft, (1, 1), (3, 1))
        right = register_allowed_link!(model, spacecraft, (2, 1), (4, 1))
        fill!(model.graph.available, true)
        work = matching_workspace!(model, policy)
        work.score .= [10.0, 6.0, 6.0]
        @test select_interlinks!(model, policy) == [left, right]
        @test model.graph.active == [false, true, true]
        for scale in (1.0, 1e-20, floatmin(Float64))
            work.score .= [2 * scale, scale, nextfloat(scale)]
            @test select_interlinks!(model, policy) == [left, right]
            work.score[3] = prevfloat(scale)
            @test select_interlinks!(model, policy) == [central]
        end
        work.score .= [10.0, 6.0, 6.0]
        copied_model, copied_policy = deepcopy((model, policy))
        @test copied_policy.matching.graph === copied_model.graph
        @test copied_policy.matching.graph !== model.graph
        @test select_interlinks!(copied_model, copied_policy) == [left, right]
        extra = register_allowed_link!(model, spacecraft, (1, 2), (5, 1))
        @test matching_workspace!(model, policy).score == [10.0, 6.0, 6.0, 0.0]
        model.graph.available[end] = true
        matching_workspace!(model, policy).score[end] = 1.0
        @test select_interlinks!(model, policy) == sort([left, right, extra])
        @test isempty(select_interlinks!(model, SchedulingPolicyModel(:gve_sma)))

        triangle = ConstellationInterLinkModel()
        for (i, j) in ((1, 2), (2, 3), (1, 3))
            register_allowed_link!(triangle, spacecraft, (i, 1), (j, 1))
        end
        fill!(triangle.graph.available, true)
        triangle_policy = SchedulingPolicyModel(:gve_sma)
        matching_workspace!(triangle, triangle_policy).score .= [3.0, 2.0, 1.0]
        @test select_interlinks!(triangle, triangle_policy) == [((1, 1), (2, 1))]
        @test any(triangle_policy.matching.cyclic)
        PythonCall.GIL.@unlock begin
            @test_throws ArgumentError select_interlinks!(triangle, triangle_policy)
        end
        @test select_interlinks!(triangle, triangle_policy) == [((1, 1), (2, 1))]
        triangle.graph.available[1] = false
        @test select_interlinks!(triangle, triangle_policy) == [((2, 1), (3, 1))]
        @test !any(triangle_policy.matching.cyclic)

        rng = MersenneTwister(718)
        pairs = [(i, j) for i in 1:5 for j in (i + 1):6]
        six = [SpacecraftModel() for _ in 1:6]
        for _ in 1:40
            random_model = ConstellationInterLinkModel()
            for (i, j) in shuffle(rng, pairs)[1:8]
                register_allowed_link!(random_model, six, (i, 1), (j, 1))
            end
            random_model.graph.available .= rand(rng, Bool, 8)
            random_policy = SchedulingPolicyModel(:gve_sma)
            scores = matching_workspace!(random_model, random_policy).score
            scores .= rand(rng, [-1.0, 0.0, 1e-20, 1.0, nextfloat(1.0), 3.0], 8)
            expected = exhaustive_matching_score(random_model, random_policy)
            selected = select_interlinks!(random_model, random_policy)
            actual = sum((Rational{BigInt}(scores[random_model.graph.edge_index_lookup[key]])
                for key in selected); init=zero(Rational{BigInt}))
            @test actual == expected
            @test length(Set(endpoint for key in selected for endpoint in key)) == 2 * length(selected)
            @test selected == select_interlinks!(random_model, random_policy)
        end
    end

    @testset "Forces and GVE policies" begin
        spacecraft, state, model, key, runtime = endpoint_fixture()
        mu = runtime.args.environment_model.planet.μ
        parameters = model.graph.edge_properties[1]
        first_force = force_on_endpoint(state.sc[1], state.sc[2], parameters)
        second_force = force_on_endpoint(state.sc[2], state.sc[1], parameters)
        @test first_force == -second_force
        @test norm(first_force) ≈ parameters.B * parameters.P / 299_792_458.0
        @test iszero(force_on_endpoint(state.sc[1], state.sc[1], parameters))
        @test iszero(laser_force_on_spacecraft(model, state, 1))
        update_availability!(model, spacecraft, state)
        rates = [2 * semimajor_axis(state.sc[i], mu)^2 / mu *
            dot(state.sc[i].vel, i == 1 ? first_force : second_force) / state.sc[i].mass for i in 1:2]
        for (policy, expected) in (
            (SchedulingPolicyModel(:gve_sma), rates[1]),
            (SchedulingPolicyModel(:gve_sma; target_idx=2), rates[2]),
            (SchedulingPolicyModel(:gve_sma; target_idx=nothing), sum(rates)),
            (SchedulingPolicyModel(:gve_sma; target_idx=99, target_idxs=[2, 1]), sum(rates)),
            (SchedulingPolicyModel(:gve_sma_first_index; target_idx=2), rates[1]),
        )
            score_candidates!(model, policy, state, mu; candidate_force=force_on_endpoint)
            @test only(matching_workspace!(model, policy).score) ≈ expected
        end
        eccentric_policy = SchedulingPolicyModel(:gve_eccentricity)
        score_candidates!(model, eccentric_policy, state, mu; candidate_force=force_on_endpoint)
        plus, minus = copy(state.sc[1]), copy(state.sc[1])
        plus.vel .+= first_force / plus.mass
        minus.vel .-= first_force / minus.mass
        @test isapprox(only(eccentric_policy.matching.score),
            (eccentricity(plus, mu) - eccentricity(minus, mu)) / 2; rtol=1e-5, atol=1e-14)
        for targets in (Int[], [0], [1, 1])
            @test_throws ArgumentError SchedulingPolicyModel(:gve_sma; target_idxs=targets)
        end
        @test_throws ArgumentError SchedulingPolicyModel(:invalid)
        @test_throws ArgumentError SchedulingPolicyModel(:gve_eccentricity; target_idxs=[1])
        @test_throws ArgumentError score_candidates!(model, SchedulingPolicyModel(:gve_sma; target_idxs=[3]),
            state, mu; candidate_force=force_on_endpoint)
        @test_throws ArgumentError score_candidates!(model, SchedulingPolicyModel(:gve_sma),
            state, 0.0; candidate_force=force_on_endpoint)
        @test schedule_interlinks!(model, runtime.args.scheduling_policy_model, spacecraft, state, mu;
            candidate_force=force_on_endpoint) == [key]
        derivative = zero(state)
        SE._apply_interlink_rhs!(derivative, state, runtime)
        @test derivative.sc[1].laser_dv ≈ first_force / state.sc[1].mass
        @test derivative.sc[1].vel == derivative.sc[1].laser_dv
        @test derivative.sc[1].laser_delta_sma ≈ rates[1]
        @test state.sc[1].mass * derivative.sc[1].vel + state.sc[2].mass * derivative.sc[2].vel ≈ zeros(3)
        state.sc[2].pos[3] = 1000.0
        @test laser_force_on_spacecraft(model, state, 1) != first_force
        runtime.is_active[2] = false
        @test isempty(schedule_interlinks!(model, runtime.args.scheduling_policy_model, spacecraft, state, mu;
            candidate_force=force_on_endpoint, is_active=runtime.is_active))
        @test iszero(laser_force_on_spacecraft(model, state, 1))
    end

    @testset "Continuous positions and switching cache refresh" begin
        spacecraft, state, model, key, runtime = endpoint_fixture(range=1000.5)
        state.sc[1].vel[2] = 1.0
        state.sc[2].vel[2] = 0.0
        problem = ODEProblem(continuous_laser_fixture!, state, (0.0, 2.0), runtime;
            callback=interlink_scheduler_callback(candidate_force=force_on_endpoint))
        solution = solve(problem, Tsit5(); adaptive=false, dt=1.0)
        @test [sample.time for sample in model.history] == [0.0, 1.0, 2.0]
        @test model.history[1].active == [key]
        @test all(isempty(sample.active) for sample in model.history[2:end])
        for index in 1:2
            acceleration = force_on_endpoint(state.sc[index], state.sc[3 - index],
                model.graph.edge_properties[1]) / state.sc[index].mass
            @test isapprox(solution.u[end].sc[index].vel, state.sc[index].vel + acceleration; atol=1e-10, rtol=0)
            @test isapprox(solution.u[end].sc[index].pos, state.sc[index].pos + 2 * state.sc[index].vel + 1.5 * acceleration; atol=1e-9, rtol=0)
            @test isapprox(solution.u[end].sc[index].laser_dv, acceleration; atol=1e-12, rtol=0)
        end
    end

    @testset "Engine, RHS routes, isolation and rejected trials" begin
        args = gve_sma_configuration(duration=20.0)
        solution = run_simulation(args; return_solution=true)
        model = solution.prob.p.args.interlink_model
        policy = solution.prob.p.args.scheduling_policy_model
        @test solution.t[end] == 20.0
        @test all(isfinite, solution.u[end])
        @test !isempty(model.history[1].active)
        @test length(model.history) == solution.stats.naccept + 1
        @test solution.prob.f.jac_prototype === nothing
        @test isempty(args.interlink_model.history)
        @test args.scheduling_policy_model.matching === nothing
        @test policy.matching.graph === model.graph
        @test !any(args.interlink_model.graph.active)
        history_before = [(sample.time, copy(sample.active)) for sample in model.history]
        flags_before = (copy(model.graph.available), copy(model.graph.active), copy(policy.matching.score))
        state_before = copy(solution.u[end])
        for time in (10.0, 2.0, 19.0)
            SE.spacecraft_dynamics!(similar(state_before), state_before, solution.prob.p, time)
        end
        @test state_before == solution.u[end]
        @test [(sample.time, sample.active) for sample in model.history] == history_before
        @test (model.graph.available, model.graph.active, policy.matching.score) == flags_before
        for (mode, expected) in (("satellite", :satellite_batch), ("flat", :flat_constellation_effector_queue))
            routed = withenv("SPACEAGORA_RHS_EXECUTION_MODE" => mode) do
                run_simulation(args; return_solution=true)
            end
            p = routed.prob.p
            @test SE._rhs_execution_plan(p.args, p, p.args.dynamics_model.dynamic_effectors, 3).mode == expected
            @test isapprox(routed.u[end], solution.u[end]; rtol=1e-12, atol=1e-8)
            @test any(!iszero(sc.laser_delta_sma) for sc in routed.u[end].sc)
        end
        stiff_args = SM.SimConfig._with_configuration(args; solver_config=SM.SolverConfig(solver_mode=:rodas5p))
        stiff = run_simulation(stiff_args; return_solution=true)
        @test stiff.prob.f.jac_prototype === nothing
        @test isapprox(stiff.u[end], solution.u[end]; rtol=1e-8, atol=1e-5)
        multi_target = gve_sma_configuration(duration=20.0)
        register_allowed_link!(multi_target.interlink_model, multi_target.dynamics_model.spacecraft,
            (2, 1), (3, 1); parameters=InterLinkParameters(range=500e3))
        multi_target = SM.SimConfig._with_configuration(multi_target;
            scheduling_policy_model=SchedulingPolicyModel(:gve_sma; target_idxs=[3, 1, 2]))
        multi_solution = run_simulation(multi_target; return_solution=true)
        multi_ran = multi_solution.prob.p.args
        @test multi_solution.t[end] == 20.0
        @test all(isfinite, multi_solution.u[end])
        @test any(multi_ran.scheduling_policy_model.matching.cyclic)
        @test length(multi_ran.interlink_model.history) == multi_solution.stats.naccept + 1
        @test all(length(sample.active) == 1 for sample in multi_ran.interlink_model.history)
        orientation_args = SM.SimConfig._with_configuration(args;
            mission_configuration=MissionConfiguration(mission_type=MissionTime, mission_time=20.0,
                orientation_sim=true, data_rate=10.0))
        oriented = run_simulation(orientation_args; return_solution=true)
        for i in 1:3
            @test norm(oriented.u[end].sc[i].q) ≈ 1.0
            @test isapprox(oriented.u[end].sc[i].pos, solution.u[end].sc[i].pos; atol=1e-6, rtol=1e-12)
            @test isapprox(oriented.u[end].sc[i].laser_dv, solution.u[end].sc[i].laser_dv; atol=1e-10)
        end
        problem = deepcopy(solution.prob)
        empty!(problem.p.args.interlink_model.history)
        rejected = solve(problem, Tsit5(); dt=20.0, dtmax=20.0, reltol=1e-13, abstol=1e-13)
        rejected_history = rejected.prob.p.args.interlink_model.history
        @test rejected.stats.nreject > 0
        @test length(rejected_history) == rejected.stats.naccept + 1
        @test all(diff([sample.time for sample in rejected_history]) .> 0)
        reference_problem = deepcopy(solution.prob)
        empty!(reference_problem.p.args.interlink_model.history)
        reference = solve(reference_problem, Tsit5(); dt=0.1, dtmax=0.25, reltol=1e-13, abstol=1e-13)
        for index in 1:3
            @test isapprox(rejected.u[end].sc[index].laser_dv, reference.u[end].sc[index].laser_dv; atol=1e-9)
            @test isapprox(rejected.u[end].sc[index].laser_delta_sma, reference.u[end].sc[index].laser_delta_sma; atol=1e-6)
        end
        baseline_args = gve_sma_configuration(duration=20.0, with_interlinks=false)
        baseline = run_simulation(baseline_args; return_solution=true)
        zero_power = run_simulation(gve_sma_configuration(duration=20.0, power=0.0); return_solution=true)
        @test !hasproperty(baseline.u[end].sc[1], :laser_dv)
        @test !any(field -> field.name in (:laser_dv, :laser_delta_sma), default_save_fields(baseline_args))
        for i in 1:3
            @test isapprox(baseline.u[end].sc[i].pos, zero_power.u[end].sc[i].pos; atol=1e-6, rtol=1e-12)
            @test iszero(zero_power.u[end].sc[i].laser_dv)
            @test iszero(zero_power.u[end].sc[i].laser_delta_sma)
        end
        run_simulation(args; isolate_state=false)
        first_history = [(sample.time, copy(sample.active)) for sample in args.interlink_model.history]
        run_simulation(args; isolate_state=false)
        @test first_history == [(sample.time, sample.active) for sample in args.interlink_model.history]
    end

    @testset "Results, checkpoint continuation and invalid combinations" begin
        mktempdir() do directory
            args = with_settings(gve_sma_configuration(duration=20.0);
                results=true, results_directory=directory)
            result = run_simulation(args; return_results=true, return_solution=true)
            csv = CSV.read(joinpath(directory, "simulation_results.csv"), DataFrame)
            @test names(result.table) == names(csv)
            sma_columns = ["sc$(i)_laser_delta_sma" for i in 1:3]
            dv_columns = ["sc$(i)_laser_dv_$(axis)" for i in 1:3 for axis in 1:3]
            @test filter(name -> occursin("laser_delta_sma", name), names(csv)) == sma_columns
            @test filter(name -> occursin("laser_dv", name), names(csv)) == dv_columns
            @test [last(csv[!, name]) for name in sma_columns] ≈ [sc.laser_delta_sma for sc in result.solution.u[end].sc]
            @test [last(csv[!, name]) for name in dv_columns] ≈ vcat([collect(sc.laser_dv) for sc in result.solution.u[end].sc]...)
            @test !isempty(result.configuration.interlink_model.history)
            checkpointed = with_settings(args; results=false, checkpoint_enabled=true, checkpoint_interval_s=10.0)
            segments = run_simulation(checkpointed; return_results=true, return_solution=true)
            times = [sample.time for sample in segments.configuration.interlink_model.history]
            @test all(diff(times) .> 0)
            @test first(times) == 0.0 && last(times) == 20.0
            resumed_args = with_settings(SM.SimConfig._with_configuration(checkpointed;
                mission_configuration=gve_sma_configuration(duration=30.0).mission_configuration);
                resume_from_checkpoint=true)
            resumed = run_simulation(resumed_args; return_results=true, return_solution=true)
            @test first(resumed.configuration.interlink_model.history).time == 20.0
            @test last(resumed.configuration.interlink_model.history).time == 30.0
            @test resumed.solution.prob.u0.sc[1].laser_delta_sma ≈ segments.solution.u[end].sc[1].laser_delta_sma
            full = run_simulation(gve_sma_configuration(duration=30.0); return_solution=true)
            @test isapprox(resumed.solution.u[end], full.u[end]; atol=1e-5, rtol=1e-9)
            mismatch = SM.SimConfig._with_configuration(resumed_args; interlink_model=nothing)
            @test_throws ArgumentError run_simulation(mismatch)
        end
        args = gve_sma_configuration(duration=20.0)
        for mode in (:split_imex, :multirate, :gravity_backbone_split)
            @test_throws ArgumentError run_simulation(SM.SimConfig._with_configuration(args;
                solver_config=SM.SolverConfig(solver_mode=mode)))
        end
        @test_throws ArgumentError run_simulation(SM.SimConfig._with_configuration(args;
            interlink_model=ConstellationInterLinkModel(constellation_link_type=:radio)))
        @test_throws ArgumentError run_constellation_ensemble(args)
        @test_throws ArgumentError run_constellation_ensemble(args; allow_gnc_effectors=true)
        @test_throws ArgumentError SE._validate_interlinks!(SM.SimConfig._with_configuration(args;
            external_propagators=(InterlinkExternalOwner(),)), :tsit5)
        invalid = deepcopy(args)
        invalid.dynamics_model.spacecraft[1].n_terminal = 0
        @test_throws ArgumentError run_simulation(invalid)
        articulated = deepcopy(args)
        sc = articulated.dynamics_model.spacecraft[1]
        panel = Link(m=1.0)
        push!(sc.links, panel)
        push!(sc.joints, Joint(sc.root, panel; joint_type=:hinge, axis=[0.0, 0.0, 1.0]))
        @test_throws ArgumentError SE._validate_interlinks!(articulated, :tsit5)
    end

    @testset "One-hour step-cap convergence" begin
        coarse = run_gve_sma_case(dt_max_orbit=10.0, save_results=false)
        fine = run_gve_sma_case(dt_max_orbit=5.0, save_results=false)
        baseline = run_simulation(gve_sma_configuration(with_interlinks=false); return_solution=true)
        for (solution, step_cap) in ((coarse, 10.0), (fine, 5.0))
            history = solution.prob.p.args.interlink_model.history
            intervals = diff([sample.time for sample in history])
            @test solution.t[end] == 3600.0
            @test all(isfinite, solution.u[end])
            @test length(history) == solution.stats.naccept + 1
            @test all(intervals .> 0)
            @test 0.9 * step_cap < maximum(intervals) <= step_cap + 1e-9
            @test !isempty(interlink_switches(history))
            @test sum(sc.laser_delta_sma for sc in solution.u[end].sc) > 0
            for index in 1:3
                mu = solution.prob.p.args.environment_model.planet.μ
                measured = semimajor_axis(solution.u[end].sc[index], mu) - semimajor_axis(baseline.u[end].sc[index], mu)
                @test isapprox(measured, solution.u[end].sc[index].laser_delta_sma; atol=0.01, rtol=0)
            end
        end
        @test maximum(abs(coarse.u[end].sc[i].laser_delta_sma - fine.u[end].sc[i].laser_delta_sma) for i in 1:3) <= 1.0
        @test maximum(norm(coarse.u[end].sc[i].laser_dv - fine.u[end].sc[i].laser_dv) for i in 1:3) <= 1e-3
        @test maximum(norm(coarse.u[end].sc[i].pos - fine.u[end].sc[i].pos) for i in 1:3) <= 5.0
        coarse_switches = interlink_switches(coarse.prob.p.args.interlink_model.history)
        fine_switches = interlink_switches(fine.prob.p.args.interlink_model.history)
        @test [sample.active for sample in coarse_switches] == [sample.active for sample in fine_switches]
        @test maximum(abs(a.time - b.time) for (a, b) in zip(coarse_switches, fine_switches)) <= 10.0
    end
end
end

end # module InterlinkTests
