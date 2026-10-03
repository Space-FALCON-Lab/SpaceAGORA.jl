using Test
using ComponentArrays
using OrdinaryDiffEq
using Random
include(joinpath(@__DIR__, "..", "II_examples", "gve_sma_interlinks.jl"))
include(joinpath(@__DIR__, "scalability.jl"))
include(joinpath(@__DIR__, "walker.jl"))

function endpoint_fixture(; eligibility=nothing)
    spacecraft = [SpacecraftModel(id=index) for index in 1:2]
    state = ComponentVector(sc=[
        (pos=[7e6, 1000.0, 0.0], vel=[0.0, 7500.0, 0.0], mass=100.0,
            laser_dv=zeros(3), laser_delta_sma=0.0),
        (pos=[7e6, 0.0, 0.0], vel=[0.0, 7500.0, 0.0], mass=200.0,
            laser_dv=zeros(3), laser_delta_sma=0.0)])
    model = InterLinkModel(; eligibility)
    key = register_candidate!(model, spacecraft, (1, 1), (2, 1))
    runtime = (args=(interlink_model=model, scheduling_policy_model=SchedulingPolicyModel(),
        dynamics_model=(spacecraft=spacecraft,), environment_model=(planet=make_no_gram_planet(:earth),)),
        is_active=trues(2))
    return spacecraft, state, model, key, runtime
end

function continuous_laser_fixture!(derivative, state, parameters, time)
    derivative .= 0.0
    for satellite in eachindex(state.sc)
        derivative.sc[satellite].pos .= state.sc[satellite].vel
    end
    SimulationEngine._apply_interlink_rhs!(derivative, state, parameters)
    return nothing
end

semimajor_axis(state, mu) = inv(2 / norm(state.pos) - dot(state.vel, state.vel) / mu)

@testset "V5 interlink reconstruction" begin
    @testset "Registration and availability" begin
        spacecraft = [SpacecraftModel(id=100 + index) for index in 1:4]
        @test spacecraft[1].n_terminal == 1
        @test spacecraft[1].battery_energy_index == 100.0
        @test spacecraft[1].tempurature_index == 100.0
        @test SpacecraftModel(n_terminal=0).n_terminal == 0
        @test_throws ArgumentError SpacecraftModel(n_terminal=-1)
        @test_throws ArgumentError SpacecraftModel(battery_energy_index=101)
        @test_throws ArgumentError SpacecraftModel(tempurature_index=NaN)
        @test_throws ArgumentError InterLinkParameters(P=-1)
        @test_throws ArgumentError InterLinkParameters(B=0)
        @test_throws ArgumentError InterLinkParameters(range=0)
        @test_throws ArgumentError InterLinkModel(active_link_penalty=-1)
        @test_throws ArgumentError InterLinkModel(battery_energy_threshold=101)
        model = InterLinkModel(tempurature_threshold=25.0)
        key = register_candidate!(model, spacecraft, (2, 1), (1, 1);
            parameters=InterLinkParameters(range=10.0))
        @test key == ((1, 1), (2, 1))
        @test spacecraft[1].id == 101
        connection = model.linkgraph[key]
        @test !connection.state.available && !connection.state.active && connection.state.score == 0.0
        @test_throws ArgumentError register_candidate!(model, spacecraft, (1, 1), (2, 1))
        @test_throws ArgumentError register_candidate!(model, spacecraft, (1, 2), (3, 1))
        @test_throws ArgumentError register_candidate!(model, spacecraft, (1, 1), (5, 1))
        @test_throws ArgumentError register_candidate!(model, spacecraft, (1, 1), (1, 1))
        state = (sc=[(pos=[Float64(index), 0.0, 0.0],) for index in 1:4],)
        @test update_availability!(model, key, spacecraft, state)
        connection.state.active = true
        @test update_availability!(model, key, spacecraft, state)
        state.sc[2].pos[1] = 11.0
        @test update_availability!(model, key, spacecraft, state)
        state.sc[2].pos[1] = 11.1
        @test !update_availability!(model, key, spacecraft, state)
        state.sc[2].pos[1] = 2.0
        spacecraft[1].battery_energy_index = 49.0
        @test !update_availability!(model, key, spacecraft, state)
        spacecraft[1].battery_energy_index = 50.0
        @test update_availability!(model, key, spacecraft, state)
        spacecraft[2].tempurature_index = 24.0
        @test !update_availability!(model, key, spacecraft, state)
        spacecraft[2].tempurature_index = 100.0
        push!(model.forbidden_pairs, (1, 2))
        @test !update_availability!(model, key, spacecraft, state)
        empty!(model.forbidden_pairs)
        @test !update_availability!(model, key, spacecraft, state; is_active=[true, false, true, true])
        model.eligibility = (key, spacecraft, state) -> false
        @test !update_availability!(model, key, spacecraft, state)
    end

    @testset "Exact penalized terminal matching" begin
        spacecraft = [SpacecraftModel(n_terminal=2) for _ in 1:5]
        model = InterLinkModel(active_link_penalty=1.0)
        central = register_candidate!(model, spacecraft, (1, 1), (2, 1))
        left = register_candidate!(model, spacecraft, (1, 1), (3, 1))
        right = register_candidate!(model, spacecraft, (2, 1), (4, 1))
        for (key, score) in ((central, 10.0), (left, 6.0), (right, 6.0))
            model.linkgraph[key].state.available = true
            model.linkgraph[key].state.score = score
        end
        @test select_interlinks!(model) == [left, right]
        @test select_interlinks!(model) == [left, right]
        @test !model.linkgraph[central].state.active
        model.active_link_penalty = 10.0
        @test isempty(select_interlinks!(model))
        @test all(!connection.state.active for connection in values(model.linkgraph))
        @test isempty(select_interlinks!(InterLinkModel()))
        model = InterLinkModel()
        first_key = register_candidate!(model, spacecraft, (1, 1), (2, 1))
        second_key = register_candidate!(model, spacecraft, (1, 2), (3, 1))
        for connection in values(model.linkgraph)
            connection.state.available = true
            connection.state.score = 1.0
        end
        @test select_interlinks!(model) == [first_key, second_key]

        generator = MersenneTwister(52)
        for trial in 1:8
            model = InterLinkModel(active_link_penalty=1.0)
            candidates = [register_candidate!(model, spacecraft, (first_satellite, 1), (second_satellite, 1))
                for first_satellite in 1:4 for second_satellite in (first_satellite + 1):5]
            for connection in values(model.linkgraph)
                connection.state.available = rand(generator, Bool)
                connection.state.score = rand(generator, -2:10)
            end
            best_score = 0.0
            for mask in 0:((1 << length(candidates)) - 1)
                subset = [key for (index, key) in enumerate(candidates) if !iszero(mask & (1 << (index - 1)))]
                endpoints = [endpoint for key in subset for endpoint in key]
                length(unique(endpoints)) == length(endpoints) || continue
                all(model.linkgraph[key].state.available for key in subset) || continue
                best_score = max(best_score, sum((model.linkgraph[key].state.score - 1.0 for key in subset); init=0.0))
            end
            selected = select_interlinks!(model)
            @test sum((model.linkgraph[key].state.score - 1.0 for key in selected); init=0.0) == best_score
            @test selected == select_interlinks!(model)
        end
    end

    @testset "Endpoint forces and physical scores" begin
        spacecraft, state, model, key, runtime = endpoint_fixture()
        mu = runtime.args.environment_model.planet.μ
        @test schedule_interlinks!(model, SchedulingPolicyModel(), spacecraft, state, mu) == [key]
        forces = [laser_force_on_spacecraft(model, state, index) for index in 1:2]
        @test norm(forces[1]) ≈ 1e6 / 299_792_458.0
        @test forces[1] == -forces[2]
        @test forces[1][2] > 0.0
        @test forces[1] / state.sc[1].mass ≈ -2 * forces[2] / state.sc[2].mass
        @test reverse([laser_force_on_spacecraft(model, state, index) for index in 2:-1:1]) == forces
        @test force_on_endpoint(state.sc[1], state.sc[2], InterLinkParameters(P=20_000.0)) ≈ 2 * forces[1]
        @test force_on_endpoint(state.sc[1], state.sc[2], InterLinkParameters(B=200.0)) ≈ 2 * forces[1]
        original = copy(state)
        derivative = zero(state)
        SimulationEngine._apply_interlink_rhs!(derivative, state, runtime)
        @test state == original
        @test derivative.sc[1].vel == forces[1] / 100.0
        @test derivative.sc[2].vel == forces[2] / 200.0
        @test derivative.sc[1].laser_dv == derivative.sc[1].vel
        @test derivative.sc[1].laser_delta_sma ≈ model.linkgraph[key].state.score
        @test isempty(model.history)
        for index in 1:2
            forward, backward = copy(state), copy(state)
            acceleration = forces[index] / state.sc[index].mass
            forward.sc[index].vel .+= acceleration
            backward.sc[index].vel .-= acceleration
            finite_difference = (semimajor_axis(forward.sc[index], mu) - semimajor_axis(backward.sc[index], mu)) / 2
            @test isapprox(semimajor_axis_rate(state.sc[index], forces[index], mu), finite_difference; atol=1e-7)
        end
        stage = copy(state)
        stage.sc[2].pos[1] += 1000.0
        @test laser_force_on_spacecraft(model, stage, 1) != forces[1]
        @test model.linkgraph[key].state.active
        stage.sc[2].pos[1] += 300e3
        @test norm(laser_force_on_spacecraft(model, stage, 1)) > 0
        @test isempty(schedule_interlinks!(model, SchedulingPolicyModel(), spacecraft, stage, mu))
        @test iszero(laser_force_on_spacecraft(model, stage, 1))
        @test SchedulingPolicyModel().target_idx == 1
        @test_throws ArgumentError SchedulingPolicyModel(:gve_sma; target_idx=0)
        @test isempty(schedule_interlinks!(model, SchedulingPolicyModel(:gve_sma; target_idx=2), spacecraft, state, mu))
        @test model.linkgraph[key].state.score ≈ semimajor_axis_rate(state.sc[2], forces[2], mu)
        @test_throws ArgumentError score_candidates!(model, SchedulingPolicyModel(target_idx=3), state, mu)
    end

    @testset "Continuous positions and switching cache refresh" begin
        spacecraft, state, model, key, runtime = endpoint_fixture(
            eligibility=(key, spacecraft, state) -> state.sc[1].pos[2] < 1001.0)
        for current in state.sc
            current.vel[2] = 1.0
        end
        problem = ODEProblem(continuous_laser_fixture!, state, (0.0, 2.0), runtime;
            callback=interlink_scheduler_callback())
        solution = solve(problem, Tsit5(); adaptive=false, dt=1.0)
        @test [sample.time for sample in model.history] == [0.0, 1.0, 2.0]
        @test model.history[1].active == [key]
        @test all(isempty(sample.active) for sample in model.history[2:end])
        for index in 1:2
            acceleration = force_on_endpoint(state.sc[index], state.sc[3 - index], model.linkgraph[key].parameters) / state.sc[index].mass
            @test isapprox(solution.u[end].sc[index].vel, state.sc[index].vel + acceleration; atol=1e-10, rtol=0)
            @test isapprox(solution.u[end].sc[index].pos, state.sc[index].pos + 2 * state.sc[index].vel + 1.5 * acceleration; atol=1e-9, rtol=0)
            @test isapprox(solution.u[end].sc[index].laser_dv, acceleration; atol=1e-12, rtol=0)
        end
    end

    @testset "Engine integration and rejected trials" begin
        args = gve_sma_configuration(duration=20.0)
        solution = run_simulation(args; return_solution=true)
        model = solution.prob.p.args.interlink_model
        @test solution.t[end] == 20.0
        @test all(isfinite, solution.u[end])
        @test !isempty(model.history[1].active)
        @test length(model.history) == solution.destats.naccept + 1
        @test isempty(args.interlink_model.history)
        @test all(!connection.state.active for connection in values(args.interlink_model.linkgraph))
        history_before = deepcopy(model.history)
        flags_before = [(connection.state.available, connection.state.active, connection.state.score)
            for connection in values(model.linkgraph)]
        state_before = copy(solution.u[end])
        for time in (10.0, 2.0, 19.0)
            SimulationEngine.spacecraft_dynamics!(similar(state_before), state_before, solution.prob.p, time)
        end
        @test state_before == solution.u[end]
        @test [(sample.time, sample.active) for sample in model.history] == [(sample.time, sample.active) for sample in history_before]
        @test [(connection.state.available, connection.state.active, connection.state.score)
            for connection in values(model.linkgraph)] == flags_before
        problem = deepcopy(solution.prob)
        empty!(problem.p.args.interlink_model.history)
        rejected = solve(problem, Tsit5(); dt=20.0, dtmax=20.0, reltol=1e-13, abstol=1e-13)
        rejected_history = rejected.prob.p.args.interlink_model.history
        @test rejected.destats.nreject > 0
        @test length(rejected_history) == rejected.destats.naccept + 1
        @test all(diff([sample.time for sample in rejected_history]) .> 0)
        reference_problem = deepcopy(solution.prob)
        empty!(reference_problem.p.args.interlink_model.history)
        reference = solve(reference_problem, Tsit5(); dt=0.1, dtmax=0.25, reltol=1e-13, abstol=1e-13)
        for index in 1:3
            @test isapprox(rejected.u[end].sc[index].laser_dv, reference.u[end].sc[index].laser_dv; atol=1e-9)
            @test isapprox(rejected.u[end].sc[index].laser_delta_sma, reference.u[end].sc[index].laser_delta_sma; atol=1e-6)
        end
        baseline = run_simulation(gve_sma_configuration(duration=20.0, with_interlinks=false); return_solution=true)
        @test baseline.t[end] == 20.0
        @test !hasproperty(baseline.u[end].sc[1], :laser_dv)
    end

    @testset "One-hour gve_sma step-cap comparison" begin
        coarse = run_gve_sma_case(dt_max_orbit=10.0)
        fine = run_gve_sma_case(dt_max_orbit=5.0)
        baseline = run_simulation(gve_sma_configuration(with_interlinks=false); return_solution=true)
        for (solution, step_cap) in ((coarse, 10.0), (fine, 5.0))
            args = solution.prob.p.args
            history = args.interlink_model.history
            intervals = diff([sample.time for sample in history])
            @test solution.t[end] == 3600.0
            @test all(isfinite, solution.u[end])
            @test length(history) == solution.destats.naccept + 1
            @test all(intervals .> 0)
            @test 0.9 * step_cap < maximum(intervals) <= step_cap + 1e-9
            @test !isempty(interlink_switches(history))
            @test sum(current.laser_delta_sma for current in solution.u[end].sc) > 0
            @test all(sample -> length(unique([endpoint for key in sample.active for endpoint in key])) == 2 * length(sample.active), history)
            for index in 1:3
                measured = semimajor_axis(solution.u[end].sc[index], args.environment_model.planet.μ) -
                    semimajor_axis(baseline.u[end].sc[index], args.environment_model.planet.μ)
                @test isapprox(measured, solution.u[end].sc[index].laser_delta_sma; atol=0.01, rtol=0)
            end
            output = args.simulation_settings.results_directory
            schedule = CSV.read(joinpath(output, "accepted_step_schedule.csv"), DataFrame)
            @test schedule.time_s == [sample.time for sample in history]
            csv_paths = filter(path -> endswith(path, ".csv") && basename(path) != "accepted_step_schedule.csv",
                readdir(output; join=true))
            @test !isempty(csv_paths)
            saved_diagnostics = false
            for path in csv_paths
                table = CSV.read(path, DataFrame)
                columns = filter(name -> occursin("laser_delta_sma", name), names(table))
                isempty(columns) && continue
                saved_diagnostics = true
                @test isapprox(sum(last(table[!, column]) for column in columns),
                    sum(current.laser_delta_sma for current in solution.u[end].sc); atol=1e-8)
            end
            @test saved_diagnostics
        end
        sma_error = maximum(abs(coarse.u[end].sc[index].laser_delta_sma - fine.u[end].sc[index].laser_delta_sma) for index in 1:3)
        dv_error = maximum(norm(coarse.u[end].sc[index].laser_dv - fine.u[end].sc[index].laser_dv) for index in 1:3)
        position_error = maximum(norm(coarse.u[end].sc[index].pos - fine.u[end].sc[index].pos) for index in 1:3)
        @test sma_error <= 1.0
        @test dv_error <= 1e-3
        @test position_error <= 5.0
        coarse_switches = interlink_switches(coarse.prob.p.args.interlink_model.history)
        fine_switches = interlink_switches(fine.prob.p.args.interlink_model.history)
        @test [sample.active for sample in coarse_switches] == [sample.active for sample in fine_switches]
        @test length(coarse_switches) == length(fine_switches)
        switch_error = length(coarse_switches) == length(fine_switches) ?
            maximum(abs(first_sample.time - second_sample.time) for (first_sample, second_sample) in zip(coarse_switches, fine_switches)) : Inf
        @test switch_error <= 10.0
        @printf("10 s versus 5 s: max delta-sma error=%.9f m, delta-v error=%.9g m/s, position error=%.9f m, switch error=%.6f s\n",
            sma_error, dv_error, position_error, switch_error)
        println("Switch times [s], 10 s cap: ", [sample.time for sample in coarse_switches])
        println("Switch times [s], 5 s cap: ", [sample.time for sample in fine_switches])
    end
end