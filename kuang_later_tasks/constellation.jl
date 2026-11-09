const MOMENTUM_POWER_W = 10_000.0
const MOMENTUM_MAGNIFICATION = 100.0
const MOMENTUM_RANGE_M = 200e3
const MOMENTUM_FORCE_N = MOMENTUM_MAGNIFICATION * MOMENTUM_POWER_W / 299_792_458.0
const TARGETS = 31:40

function constellation_configuration(directory; linked=true)
    spacecraft = map(1:40) do index
        helper = index <= 30
        ring = helper ? (1:30) : TARGETS
        altitude = helper ? 1000e3 : 1100e3
        phase = 360.0 * (index - first(ring)) / length(ring) + (helper ? 0.0 : 1.0)
        bus = Link(root=true, m=227.0)
        initial = InitialCondition(a=PLANET.Rp_e + altitude, e=0.0,
            i=0.0, Ω=0.0, ν=phase)
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=initial, id=index, n_terminal=1)
    end
    model = linked ? InterlinkModel(record_available=true) : nothing
    if linked
        link = LaserMomentumLink(power_w=MOMENTUM_POWER_W,
            magnification=MOMENTUM_MAGNIFICATION, range_m=MOMENTUM_RANGE_M)
        for helper in 1:30, target in TARGETS
            register_link!(model, spacecraft, (helper, 1), (target, 1), link)
        end
    end
    duration = 0.25 * 2pi * sqrt((PLANET.Rp_e + 1100e3)^3 / PLANET.μ)
    return configuration(spacecraft, model, SemimajorAxisRatePolicy(targets=TARGETS);
        duration, dtmax=10.0, output=directory, data_rate=10.0)
end

semimajor_axis(state) = inv(2 / norm(state.pos) - dot(state.vel, state.vel) / PLANET.μ)

function momentum_candidates(state, time_s)
    return DataFrame(map([(i, j) for i in 1:30 for j in TARGETS]) do (i, j)
        target = state.sc[j]
        displacement = target.pos - state.sc[i].pos
        distance = norm(displacement)
        feasible = 0 < distance <= MOMENTUM_RANGE_M
        force = MOMENTUM_FORCE_N * displacement / distance
        score = 2 * semimajor_axis(target)^2 * dot(target.vel, force) /
            (PLANET.μ * target.mass)
        (; time_s, tx=i, tx_terminal=1, rx=j, rx_terminal=1, allowed=true,
            distance_m=distance, range_m=MOMENTUM_RANGE_M,
            geometrically_feasible=feasible, available=feasible,
            hypothetical_score_m_per_s=score, score_m_per_s=feasible ? score : 0.0,
            predicted_active=feasible && score > 0)
    end)
end

function save_orbit_states(solution, path)
    rows = NamedTuple[]
    for (time_s, state) in zip(solution.t, solution.u), i in 1:40
        satellite = state.sc[i]
        push!(rows, (; time_s, spacecraft=i, mass_kg=satellite.mass,
            x_m=satellite.pos[1], y_m=satellite.pos[2], z_m=satellite.pos[3],
            vx_m_per_s=satellite.vel[1], vy_m_per_s=satellite.vel[2],
            vz_m_per_s=satellite.vel[3], sma_m=semimajor_axis(satellite)))
    end
    CSV.write(path, DataFrame(rows))
end

function run_constellation()
    directory = joinpath(OUTPUT, "constellation")
    mkpath(directory)
    args = constellation_configuration(directory)
    initial = SE.build_initial_conditions(args)
    candidates = momentum_candidates(initial, 0.0)
    CSV.write(joinpath(directory, "candidates_before_run.csv"), candidates)
    CSV.write(joinpath(directory, "spacecraft_before_run.csv"), DataFrame(
        spacecraft=1:40, role=[i <= 30 ? "helper" : "target" for i in 1:40],
        mass_kg=fill(227.0, 40), terminal_capacity=ones(Int, 40),
        initial_sma_m=[sc.initial_condition.a for sc in args.dynamics_model.spacecraft],
        phase_deg=[i <= 30 ? 12.0 * (i - 1) : 36.0 * (i - 31) + 1.0 for i in 1:40]))
    write_toml(joinpath(directory, "settings.toml"), Dict(
        "duration_s" => args.mission_configuration.mission_time,
        "target_orbits" => 0.25, "spacecraft" => 40, "allowed_pairs" => 300,
        "target_indices" => collect(TARGETS), "terminal_capacity_each" => 1,
        "power_w" => MOMENTUM_POWER_W, "magnification" => MOMENTUM_MAGNIFICATION,
        "force_magnitude_n" => MOMENTUM_FORCE_N, "range_m" => MOMENTUM_RANGE_M,
        "policy" => "SemimajorAxisRatePolicy", "record_available" => true,
        "attitude" => false, "animation" => false, "occultation" => false,
        "feasibility" => "range-only ideal momentum link", "solver" => "Tsit5",
        "reltol" => 1e-10, "abstol" => 1e-8, "dtmax_s" => 10.0,
        "score_rtol" => SCORE_RTOL, "score_atol_m_per_s" => MOMENTUM_SCORE_ATOL,
        "control_sma_drift_atol_m" => 0.02))
    summary = nothing
    @testset "Quarter-orbit forty-spacecraft run" begin
        allowed = Set(row_key(row) for row in eachrow(candidates))
        @test length(allowed) == 300
        @test Set(keys(args.interlink_model.graph.edge_index_lookup)) == allowed
        println("Running 40 spacecraft for ", args.mission_configuration.mission_time,
            " s; animation disabled.")
        solution = run_checked(args)
        run_args = solution.prob.p.args
        history = run_args.interlink_model.history
        states = history_states(solution)
        audits = NamedTuple[]
        candidate_tables = DataFrame[]
        previous = Set{LinkKey}()
        maximum_score_error = 0.0
        for (sample, state) in zip(history, states)
            table = momentum_candidates(state, sample.time)
            feasible = Set(row_key(row) for row in eachrow(table) if row.available)
            expected = Set(row_key(row) for row in eachrow(table) if row.predicted_active)
            degree = zeros(Int, 40)
            for key in feasible, endpoint in key
                degree[endpoint[1]] += 1
            end
            @test maximum(degree) <= 1
            check_matching(sample, allowed, run_args.dynamics_model.spacecraft)
            @test Set(sample.available) == feasible
            @test Set(sample.active) == expected
            expected_scores = Dict(row_key(row) => row.score_m_per_s for row in eachrow(table))
            for (key, recorded_score) in zip(sample.active, sample.scores)
                @test isapprox(recorded_score, expected_scores[key];
                    rtol=SCORE_RTOL, atol=MOMENTUM_SCORE_ATOL)
                maximum_score_error = max(maximum_score_error, abs(recorded_score - expected_scores[key]))
            end
            optimum = sum(expected_scores[key] for key in expected; init=0.0)
            @test isapprox(sum(sample.scores), optimum;
                rtol=SCORE_RTOL, atol=MOMENTUM_SCORE_ATOL)
            selected = Set(sample.active)
            table.selected = [row_key(row) in selected for row in eachrow(table)]
            push!(candidate_tables, table)
            push!(audits, (; time_s=sample.time, allowed_count=length(allowed),
                feasible_count=length(feasible), selected_count=length(selected),
                matching=matching_label(selected), selected_score_m_per_s=sum(sample.scores),
                optimum_score_m_per_s=optimum, maximum_feasible_terminal_degree=maximum(degree),
                changed=sample.time > 0 && selected != previous))
            previous = selected
        end
        audit = DataFrame(audits)
        CSV.write(joinpath(directory, "candidate_scores.csv"), vcat(candidate_tables...))
        CSV.write(joinpath(directory, "accepted_steps.csv"), audit)
        schedule = save_schedule(solution, directory)
        @test Set(schedule.window) == Set(["contact", "active"])
        @test any(>(0), audit.selected_count)
        save_orbit_states(solution, joinpath(directory, "states.csv"))

        println("Running the matching 40-spacecraft no-link control.")
        control = run_checked(constellation_configuration(directory; linked=false))
        save_orbit_states(control, joinpath(directory, "control_states.csv"))
        comparison = DataFrame(map(1:40) do i
            initial_sma = semimajor_axis(initial.sc[i])
            linked_sma = semimajor_axis(solution.u[end].sc[i])
            control_sma = semimajor_axis(control.u[end].sc[i])
            (; spacecraft=i, role=i in TARGETS ? "target" : "helper", mass_kg=227.0,
                initial_sma_m=initial_sma, linked_sma_m=linked_sma, control_sma_m=control_sma,
                control_drift_m=control_sma - initial_sma, link_effect_m=linked_sma - control_sma)
        end)
        CSV.write(joinpath(directory, "spacecraft_comparison.csv"), comparison)
        @test maximum(abs, comparison.control_drift_m) <= 0.02
        @test minimum(comparison.link_effect_m[TARGETS]) >= -0.02
        @test maximum(comparison.link_effect_m[TARGETS]) > 0.1
        summary = (duration_s=solution.t[end], samples=length(history),
            switches=schedule_changes(history), maximum_active=maximum(audit.selected_count),
            maximum_score_error_m_per_s=maximum_score_error,
            minimum_target_effect_m=minimum(comparison.link_effect_m[TARGETS]),
            maximum_target_effect_m=maximum(comparison.link_effect_m[TARGETS]),
            maximum_control_drift_m=maximum(abs, comparison.control_drift_m))
        write_toml(joinpath(directory, "summary.toml"),
            Dict(string(key) => value for (key, value) in pairs(summary)))
    end
    return summary
end
