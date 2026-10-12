using Test
using ComponentArrays
using OrdinaryDiffEq
using Random
include(joinpath(@__DIR__, "..", "II_examples", "gve_sma_interlinks.jl"))
include(joinpath(@__DIR__, "..", "II_examples", "walker_interlinks.jl"))

function endpoint_fixture(; range=200e3)
    spacecraft = [SpacecraftModel(id=index) for index in 1:2]
    state = ComponentVector(sc=[
        (pos=[7e6, 1000.0, 0.0], vel=[0.0, 7500.0, 0.0], mass=100.0,
            laser_dv=zeros(3), laser_delta_sma=0.0),
        (pos=[7e6, 0.0, 0.0], vel=[0.0, 7500.0, 0.0], mass=200.0,
            laser_dv=zeros(3), laser_delta_sma=0.0)])
    model = InterLinkModel()
    key = register_candidate!(model, spacecraft, (1, 1), (2, 1);
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
        @test update_availability!(model, key, spacecraft, state)
    end

    @testset "Walker analytic snapshots" begin
        constellation = walker_constellation(satellites=40, planes=4, range_m=20e6)
        repeated = walker_constellation(satellites=40, planes=4, range_m=20e6)
        @test length(constellation.spacecraft) == 40
        @test length(constellation.model.linkgraph) == 39
        @test constellation.metadata == repeated.metadata
        @test constellation.metadata.per_plane == 10
        @test all(vehicle.n_terminal == 1 for vehicle in constellation.spacecraft)
        initial = walker_snapshot(constellation)
        @test initial == walker_snapshot(repeated)
        radius = constellation.planet.Rp_e + constellation.metadata.altitude
        mean_motion = sqrt(constellation.planet.μ / radius^3)
        @test all(isapprox(norm(current.pos), radius; atol=1e-6, rtol=0) for current in initial.sc)
        @test all(isapprox(dot(current.pos, current.vel), 0.0; atol=1e-4) for current in initial.sc)
        state = deepcopy(initial)
        advance_walker_snapshot!(state, initial, mean_motion, 2pi / mean_motion)
        @test all(isapprox(current.pos, origin.pos; atol=1e-6, rtol=0) for (current, origin) in zip(state.sc, initial.sc))
        @test all(isapprox(current.vel, origin.vel; atol=1e-8, rtol=0) for (current, origin) in zip(state.sc, initial.sc))
        @test initial == walker_snapshot(repeated)
        benchmark = benchmark_walker(constellation; duration=120.0)
        @test benchmark.metrics.time_s == [0.0, 60.0, 120.0]
        @test all(benchmark.metrics.selected .== 1)
        @test [sample.time for sample in benchmark.history] == benchmark.metrics.time_s
        @test nrow(schedule_table(benchmark.history)) == 3
        @test benchmark.history[1].active !== benchmark.history[2].active
        @test_throws ArgumentError walker_constellation(satellites=41, planes=4)
        @test_throws ArgumentError walker_constellation(satellites=40, planes=4, target_idx=41)
        @test_throws ArgumentError walker_constellation(altitude=-1.0)
        @test_throws ArgumentError benchmark_walker(constellation; interval=0.0)
    end

    @testset "Walker force-coupled integration" begin
        constellation = walker_constellation(satellites=40, planes=4, range_m=20e6)
        engine = run_walker_engine(constellation; duration=2.0, dt_max=0.5)
        @test engine.summary.satellites == 40
        @test engine.summary.candidates == 39
        @test engine.summary.retcode == "Success"
        @test engine.solution.t[end] == 2.0
        @test engine.summary.history_samples == engine.summary.accepted_steps + 1
        @test engine.summary.nonempty_schedules == engine.summary.history_samples
        @test engine.summary.target_integrated_laser_dv_norm_m_s > 0
        @test isapprox(engine.summary.target_delta_sma_m,
            engine.summary.target_integrated_laser_delta_sma_m; atol=1e-4, rtol=0)
        @test_throws ArgumentError run_walker_engine(constellation; duration=0.0)
        @test_throws ArgumentError run_walker_engine(constellation; dt_max=0.0)
    end

    @testset "All-pairs first-index Walker integration" begin
        constellation = walker_constellation(satellites=40, planes=4, altitude=1000e3,
            range_m=20e6, all_pairs=true, policy_name=:gve_sma_first_index)
        @test length(constellation.model.linkgraph) == 780
        @test haskey(constellation.model.linkgraph, ((2, 1), (3, 1)))
        @test constellation.metadata.all_pairs
        @test constellation.metadata.altitude == 1000e3
        @test constellation.policy.policy === :gve_sma_first_index
        @test all(length(edges) == 39 for edges in constellation.model.satellite_edges)
        state = walker_snapshot(constellation)
        selected = schedule_interlinks!(constellation.model, constellation.policy,
            constellation.spacecraft, state, constellation.planet.μ; candidate_force=force_on_endpoint)
        validate_walker_selection(constellation.model, selected)
        @test 1 < length(selected) <= 20
        @test any(key[1][1] != 1 && key[2][1] != 1 for key in selected)
        original_scores = copy(constellation.model.scores)
        score_candidates!(constellation.model, SchedulingPolicyModel(:gve_sma_first_index; target_idx=40),
            state, constellation.planet.μ; candidate_force=force_on_endpoint)
        @test constellation.model.scores == original_scores
        partial = run_walker_engine(constellation; duration=20.0, wall_limit_s=1e-9)
        @test !partial.summary.completed
        @test partial.summary.termination_reason == "wall_time_limit"
        @test partial.summary.retcode == "Terminated"
        @test 0 < partial.summary.simulated_duration_s < 20
        @test partial.history[end].time == partial.solution.t[end]
        @test all(length(sample.active) > 1 for sample in partial.history)
        complete = run_walker_engine(constellation; duration=1.0, dt_max=0.5)
        @test complete.summary.completed
        @test complete.summary.termination_reason == "completed"
        @test complete.summary.scoring_policy == "gve_sma_first_index"
        @test complete.summary.all_pairs
        @test complete.summary.altitude_m == 1000e3
        @test all(length(sample.active) > 1 for sample in complete.history)
        @test_throws ArgumentError run_walker_engine(constellation; wall_limit_s=0.0)
        @test_throws ArgumentError run_walker_engine(constellation; wall_limit_s=NaN)
    end

    @testset "Indexed interlink storage" begin
        spacecraft = [SpacecraftModel(n_terminal=2) for _ in 1:4]
        state = (sc=[(pos=[7e6, 10.0 * index, 0.0], vel=[0.0, -7500.0, 0.0], mass=200.0)
            for index in 1:4],)
        model = InterLinkModel()
        first_key = register_candidate!(model, spacecraft, (2, 1), (1, 1); parameters=InterLinkParameters(range=10.0))
        second_key = register_candidate!(model, spacecraft, (3, 1), (2, 2); parameters=InterLinkParameters(range=10.0))
        @test model.edge_keys == [first_key, second_key]
        @test model.edge_index[first_key] == 1
        @test model.connections[1] === model.linkgraph[first_key]
        @test model.satellite_edges == [[1], [1, 2], [2], Int[]]
        available = update_availability!(model, spacecraft, state)
        @test available == [true, true]
        @test model.availability_transition == Int8[1, 1]
        @test update_availability!(model, spacecraft, state) === available
        @test model.availability_transition == Int8[0, 0]
        @test @allocated(update_availability!(model, spacecraft, state)) <= 1024
        state.sc[3].pos[2] += 0.1
        @test update_availability!(model, spacecraft, state) == [true, false]
        @test model.availability_transition == Int8[0, -1]
        state.sc[3].pos[2] -= 0.1
        @test update_availability!(model, second_key, spacecraft, state)
        @test model.availability_transition == Int8[0, 1]
        @test update_availability!(model, spacecraft, state; is_active=[true, false, true, true]) == [false, false]
        @test model.availability_transition == Int8[-1, -1]
        spacecraft[1].battery_energy_index = 49.0
        @test update_availability!(model, spacecraft, state) == [false, true]
        spacecraft[1].battery_energy_index = 100.0
        push!(model.forbidden_pairs, (2, 3))
        @test update_availability!(model, spacecraft, state) == [true, false]
        empty!(model.forbidden_pairs)
        update_availability!(model, spacecraft, state)
        for policy in (:gve_sma, :gve_eccentricity, :gve_sma_first_index)
            score_candidates!(model, SchedulingPolicyModel(policy), state, 3.986e14; candidate_force=force_on_endpoint)
            @test model.scores == [connection.state.score for connection in model.connections]
        end
        for connection in model.connections
            connection.state.active = true
        end
        for satellite in eachindex(spacecraft)
            expected = zeros(3)
            for (key, connection) in model.linkgraph
                connection.state.active || continue
                satellite in (key[1][1], key[2][1]) || continue
                partner = key[1][1] == satellite ? key[2][1] : key[1][1]
                expected .+= force_on_endpoint(state.sc[satellite], state.sc[partner], connection.parameters)
            end
            @test laser_force_on_spacecraft(model, state, satellite) ≈ expected
        end
        third_key = register_candidate!(model, spacecraft, (4, 2), (1, 2))
        @test model.satellite_edges[4] == [3]
        @test model.connections[3] === model.linkgraph[third_key]
        @test length(model.available) == length(model.scores) == 3
    end

    @testset "Exact terminal matching" begin
        spacecraft = [SpacecraftModel(n_terminal=2) for _ in 1:5]
        model = InterLinkModel()
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
        model.linkgraph[central].state.score = 13.0
        @test select_interlinks!(model) == [central]
        for scale in (1.0, 1e-20, floatmin(Float64))
            model.linkgraph[central].state.score = 2 * scale
            model.linkgraph[left].state.score = scale
            model.linkgraph[right].state.score = nextfloat(scale)
            @test select_interlinks!(model) == [left, right]
            model.linkgraph[right].state.score = prevfloat(scale)
            @test select_interlinks!(model) == [central]
        end
        for connection in values(model.linkgraph)
            connection.state.score -= 10.0
        end
        @test isempty(select_interlinks!(model))
        @test all(!connection.state.active for connection in values(model.linkgraph))
        @test isempty(select_interlinks!(InterLinkModel()))
        model = InterLinkModel()
        first_key = register_candidate!(model, spacecraft, (1, 1), (2, 1))
        second_key = register_candidate!(model, spacecraft, (1, 2), (3, 1))
        for connection in values(model.linkgraph)
            connection.state.available = true
            connection.state.score = 0.5
        end
        @test select_interlinks!(model) == [first_key, second_key]

        generator = MersenneTwister(52)
        for trial in 1:32
            model = InterLinkModel()
            candidates = [register_candidate!(model, spacecraft, (first_satellite, 1), (second_satellite, 1))
                for first_satellite in 1:4 for second_satellite in (first_satellite + 1):5]
            for connection in values(model.linkgraph)
                connection.state.available = rand(generator, Bool)
                connection.state.score = rand(generator, -8:40) / 4
            end
            best_score = 0.0
            for mask in 0:((1 << length(candidates)) - 1)
                subset = [key for (index, key) in enumerate(candidates) if !iszero(mask & (1 << (index - 1)))]
                endpoints = [endpoint for key in subset for endpoint in key]
                length(unique(endpoints)) == length(endpoints) || continue
                all(model.linkgraph[key].state.available for key in subset) || continue
                best_score = max(best_score, sum((model.linkgraph[key].state.score for key in subset); init=0.0))
            end
            selected = select_interlinks!(model)
            @test sum((model.linkgraph[key].state.score for key in selected); init=0.0) == best_score
            @test selected == select_interlinks!(model)
        end

        model = InterLinkModel()
        for (first_satellite, second_satellite, score) in ((1, 2, 9.0), (2, 3, 10.0),
                (1, 3, 9.0), (1, 4, 8.0), (3, 5, 8.0))
            key = register_candidate!(model, spacecraft, (first_satellite, 1), (second_satellite, 1))
            model.linkgraph[key].state.available = true
            model.linkgraph[key].state.score = score
        end
        selected = select_interlinks!(model)
        @test sum(model.linkgraph[key].state.score for key in selected) == 18.0
        @test length(unique([endpoint for key in selected for endpoint in key])) == 2 * length(selected)
        for connection in values(model.linkgraph)
            connection.state.available = false
        end
        @test isempty(select_interlinks!(model))
        @test all(!connection.state.active for connection in values(model.linkgraph))

        spacecraft = [SpacecraftModel() for _ in 1:100]
        model = InterLinkModel()
        for first_satellite in 1:99, second_satellite in (first_satellite + 1):100
            key = register_candidate!(model, spacecraft, (first_satellite, 1), (second_satellite, 1))
            model.linkgraph[key].state.available = true
            model.linkgraph[key].state.score = isodd(first_satellite) && second_satellite == first_satellite + 1 ? 2.0 : 1.0
        end
        selected = select_interlinks!(model)
        @test selected == [((satellite, 1), (satellite + 1, 1)) for satellite in 1:2:99]
        @test sum(model.linkgraph[key].state.score for key in selected) == 100.0
        @test all(connection.state.active == (key in selected) for (key, connection) in model.linkgraph)
    end

    @testset "Cached hybrid matching" begin
        spacecraft = [SpacecraftModel(n_terminal=2) for _ in 1:7]
        model = InterLinkModel()
        scale = 1e-20
        for (first, second, score) in ((1, 2, 5.0), (2, 3, 7.0), (1, 3, 6.0),
                (4, 5, scale), (5, 6, 2 * scale), (6, 7, nextfloat(scale)))
            key = register_candidate!(model, spacecraft, (first, 1), (second, 1))
            model.linkgraph[key].state.available = true
            model.linkgraph[key].state.score = score
        end
        expected = [((2, 1), (3, 1)), ((4, 1), (5, 1)), ((6, 1), (7, 1))]
        @test select_interlinks!(model) == expected
        cached = model.matching
        @test count(cached.cyclic) == 3
        @test size(cached.endpoints) == (2, 6)
        @test select_interlinks!(model) == expected
        @test model.matching === cached
        model.linkgraph[((1, 1), (3, 1))].state.available = false
        @test select_interlinks!(model) == expected
        @test !any(cached.cyclic)
        model.linkgraph[((6, 1), (7, 1))].state.score = prevfloat(scale)
        @test select_interlinks!(model) == [((2, 1), (3, 1)), ((5, 1), (6, 1))]
        @test model.matching === cached
        new_key = register_candidate!(model, spacecraft, (1, 2), (7, 2))
        @test model.matching === nothing
        model.linkgraph[new_key].state.available = true
        model.linkgraph[new_key].state.score = 1.0
        @test new_key in select_interlinks!(model)
        @test model.matching !== cached
        @test size(model.matching.endpoints) == (2, 7)

        tied = [((1, 1), (2, 1)), ((2, 1), (3, 1)), ((3, 1), (4, 1)), ((4, 1), (5, 1))]
        selected_orders = map((tied, reverse(tied))) do keys
            current = InterLinkModel()
            for key in keys
                register_candidate!(current, spacecraft, key...)
                current.linkgraph[key].state.available = true
                current.linkgraph[key].state.score = 1.0
            end
            select_interlinks!(current)
        end
        @test selected_orders[1] == selected_orders[2]

        spacecraft = [SpacecraftModel() for _ in 1:4000]
        for shape in (:star, :chain, :disjoint)
            model = InterLinkModel()
            pairs = shape === :star ? [(1, partner) for partner in 2:4000] :
                [(first, first + 1) for first in (shape === :chain ? (1:3999) : (1:2:3999))]
            for (first, second) in pairs
                key = register_candidate!(model, spacecraft, (first, 1), (second, 1))
                model.linkgraph[key].state.available = true
                model.linkgraph[key].state.score = shape === :star ? Float64(second) : (isodd(first) ? 2.0 : 1.0)
            end
            expected = shape === :star ? [((1, 1), (4000, 1))] :
                [((first, 1), (first + 1, 1)) for first in 1:2:3999]
            @test select_interlinks!(model) == expected
            @test !any(model.matching.cyclic)
            @test count(connection.state.active for connection in model.connections) == length(expected)
        end
    end

    @testset "Endpoint forces and physical scores" begin
        spacecraft, state, model, key, runtime = endpoint_fixture()
        mu = runtime.args.environment_model.planet.μ
        @test schedule_interlinks!(model, SchedulingPolicyModel(:gve_sma), spacecraft, state, mu;
            candidate_force=force_on_endpoint) == [key]
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
            score_candidates!(model, SchedulingPolicyModel(:gve_sma; target_idx=index), state, mu;
                candidate_force=force_on_endpoint)
            @test isapprox(model.linkgraph[key].state.score, finite_difference; atol=1e-7)
        end
        circular_state = (pos=[7e6, 0.0, 0.0], vel=[0.0, sqrt(mu / 7e6), 0.0], mass=100.0)
        score_candidates!(model, SchedulingPolicyModel(:gve_eccentricity), (sc=[circular_state, state.sc[2]],), mu;
            candidate_force=(target, partner, parameters) -> [0.0, 1.0, 0.0])
        @test model.linkgraph[key].state.score ≈ 2 / (circular_state.mass * circular_state.vel[2])
        stage = copy(state)
        stage.sc[2].pos[1] += 1000.0
        @test laser_force_on_spacecraft(model, stage, 1) != forces[1]
        @test model.linkgraph[key].state.active
        stage.sc[2].pos[1] += 300e3
        @test norm(laser_force_on_spacecraft(model, stage, 1)) > 0
        @test isempty(schedule_interlinks!(model, SchedulingPolicyModel(:gve_sma), spacecraft, stage, mu;
            candidate_force=force_on_endpoint))
        @test iszero(laser_force_on_spacecraft(model, stage, 1))
        @test SchedulingPolicyModel(:gve_sma).target_idx == 1
        @test_throws ArgumentError SchedulingPolicyModel(:gve_sma; target_idx=0)
        @test isempty(schedule_interlinks!(model, SchedulingPolicyModel(:gve_sma; target_idx=2), spacecraft, state, mu;
            candidate_force=force_on_endpoint))
        @test model.linkgraph[key].state.score ≈ derivative.sc[2].laser_delta_sma
        @test_throws ArgumentError score_candidates!(model, SchedulingPolicyModel(:gve_sma; target_idx=3), state, mu;
            candidate_force=force_on_endpoint)

        model.linkgraph[key].state.available = true
        all_satellites_policy = SchedulingPolicyModel(:gve_sma; target_idx=nothing)
        @test all_satellites_policy.target_idx === nothing
        score_candidates!(model, all_satellites_policy, state, mu; candidate_force=force_on_endpoint)
        both_endpoint_rates = derivative.sc[1].laser_delta_sma + derivative.sc[2].laser_delta_sma
        @test model.linkgraph[key].state.score ≈ both_endpoint_rates

        eccentricity_policy = SchedulingPolicyModel(:gve_eccentricity; target_idx=1)
        score_candidates!(model, eccentricity_policy, state, mu; candidate_force=force_on_endpoint)
        forward, backward = copy(state.sc[1]), copy(state.sc[1])
        forward.vel .+= forces[1] / forward.mass
        backward.vel .-= forces[1] / backward.mass
        forward_eccentricity = norm(cross(forward.vel, cross(forward.pos, forward.vel)) / mu - forward.pos / norm(forward.pos))
        backward_eccentricity = norm(cross(backward.vel, cross(backward.pos, backward.vel)) / mu - backward.pos / norm(backward.pos))
        @test isapprox(model.linkgraph[key].state.score, (forward_eccentricity - backward_eccentricity) / 2; rtol=1e-5)
    end

    @testset "First-index SMA policy" begin
        spacecraft, state, model, key, runtime = endpoint_fixture()
        mu = runtime.args.environment_model.planet.μ
        policy = SchedulingPolicyModel(:gve_sma_first_index)
        @test_throws ArgumentError SchedulingPolicyModel(:gve_sma_higher_orbit)
        for (first_radius, second_radius) in (
                (7.01e6, 7e6), (7e6, 7.01e6), (7e6, 7e6)), velocity in (-7500.0, 7500.0)
            state.sc[1].pos .= [first_radius, 0.0, 0.0]
            state.sc[2].pos .= [0.0, second_radius, 0.0]
            state.sc[1].vel[2] = velocity
            model.linkgraph[key].state.available = true
            score_candidates!(model, policy, state, mu; candidate_force=force_on_endpoint)
            score = model.linkgraph[key].state.score
            target = state.sc[1]
            partner = state.sc[2]
            acceleration = force_on_endpoint(target, partner, model.linkgraph[key].parameters) / target.mass
            forward, backward = copy(target), copy(target)
            forward.vel .+= acceleration
            backward.vel .-= acceleration
            expected_rate = (semimajor_axis(forward, mu) - semimajor_axis(backward, mu)) / 2
            @test isapprox(score, expected_rate; atol=1e-7)
            @test select_interlinks!(model) == (score > 0 ? [key] : typeof(key)[])
            score_candidates!(model, SchedulingPolicyModel(:gve_sma_first_index; target_idx=2), state, mu;
                candidate_force=force_on_endpoint)
            @test model.linkgraph[key].state.score == score
        end
        model.linkgraph[key].state.available = false
        score_candidates!(model, policy, state, mu;
            candidate_force=(args...) -> error("Unavailable links must not be scored"))
        @test model.linkgraph[key].state.score == 0.0
        @test isempty(select_interlinks!(model))
        @test_throws ArgumentError score_candidates!(model, policy, state, 0.0; candidate_force=force_on_endpoint)
        model.linkgraph[key].state.available = true
        @test_throws ArgumentError score_candidates!(model, policy, state, mu;
            candidate_force=(args...) -> fill(NaN, 3))

        spacecraft = [SpacecraftModel(id=10 - index) for index in 1:3]
        model = InterLinkModel()
        key = register_candidate!(model, spacecraft, (3, 1), (2, 1))
        state = (sc=[copy(state.sc[1]), copy(state.sc[1]), copy(state.sc[2])],)
        model.linkgraph[key].state.available = true
        score_candidates!(model, policy, state, mu; candidate_force=force_on_endpoint)
        score = model.linkgraph[key].state.score
        score_candidates!(model, SchedulingPolicyModel(:gve_sma; target_idx=2), state, mu;
            candidate_force=force_on_endpoint)
        @test score == model.linkgraph[key].state.score
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