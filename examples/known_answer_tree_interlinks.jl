module KnownAnswerTreeInterlinks

using SpaceAGORA
using SpaceAGORA.SimulationModel
using CSV
using DataFrames
using LinearAlgebra
using Printf
using SciMLBase: successful_retcode
using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT = joinpath(REPO_ROOT, "output", "known_answer_tree_interlinks")
const REFERENCE = "dafb5f10de8c7497efe3e898e9e2630a458b7a7c"
const MASSES_KG = (180.0, 260.0, 200.0, 300.0)
const ALTITUDES_M = (1000e3, 1000e3, 1100e3, 1100e3)
const PHASES_DEG = (0.0, -0.4, 0.6, 0.2)
const ALLOWED_PAIRS = ((1, 3), (2, 3), (2, 4))
const PREDICTED_PAIRS = ((1, 3), (2, 4))
const TARGETS = 3:4
const DURATION_S = 60.0
const LINK = LaserMomentumLink(power_w=10_000.0, magnification=100.0, range_m=200e3)
const FORCE_N = LINK.magnification * LINK.power_w / 299_792_458.0

link_key(pair) = ((pair[1], 1), (pair[2], 1))
pair_label(pair) = "$(pair[1])-$(pair[2])"
matching_label(pairs) = isempty(pairs) ? "empty" : join(pair_label.(pairs), ";")

function configuration(; output=OUTPUT, duration::Real=DURATION_S)
    isfinite(duration) && duration > 0 || throw(ArgumentError("duration must be finite and positive."))
    planet = make_no_gram_planet(:earth)
    spacecraft = map(1:4) do index
        bus = Link(root=true, m=MASSES_KG[index])
        initial = InitialCondition(a=planet.Rp_e + ALTITUDES_M[index], e=0.0,
            i=0.0, Ω=0.0, ν=PHASES_DEG[index])
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=initial, id=index, n_terminal=1)
    end
    model = InterlinkModel(record_available=true)
    for pair in ALLOWED_PAIRS
        register_link!(model, spacecraft, link_key(pair)..., LINK)
    end
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=true, generate_plots=false,
            results_directory=output, save_csv=true),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=duration, orientation_sim=false, data_rate=1.0),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10,
            abstol_orbit=1e-12, dt_max_orbit=1.0),
        solver_config=SpaceAGORA.SimulationModel.SolverConfig(solver_mode=:tsit5),
        interlink_model=model,
        scheduling_policy_model=SemimajorAxisRatePolicy(targets=TARGETS))
end

# Independent oracle: do not call the scheduler's scorer or link-force hook.
function force_on_first(state, pair)
    displacement = state.sc[pair[1]].pos - state.sc[pair[2]].pos
    return FORCE_N * displacement / norm(displacement)
end

function candidate_table(args, state)
    mu = args.environment_model.planet.μ
    return DataFrame(map([(i, j) for i in 1:3 for j in (i + 1):4]) do pair
        i, j = pair
        distance = norm(state.sc[i].pos - state.sc[j].pos)
        force = force_on_first(state, pair)
        score = 0.0
        for (index, endpoint_force) in ((i, force), (j, -force))
            index in TARGETS || continue
            endpoint = state.sc[index]
            sma = inv(2 / norm(endpoint.pos) - dot(endpoint.vel, endpoint.vel) / mu)
            score += 2 * sma^2 * dot(endpoint.vel, endpoint_force) / (mu * endpoint.mass)
        end
        allowed = pair in ALLOWED_PAIRS
        geometric = 0 < distance <= LINK.range_m
        (; a_spacecraft=i, a_terminal=1, b_spacecraft=j, b_terminal=1,
            pair=pair_label(pair), allowed, distance_m=distance, range_m=LINK.range_m,
            geometrically_feasible=geometric, available=allowed && geometric,
            hypothetical_score_m_per_s=score,
            candidate_score_m_per_s=allowed ? (geometric ? score : 0.0) : missing,
            predicted_active=pair in PREDICTED_PAIRS)
    end)
end

function enumerate_matchings(candidates)
    allowed = candidates[candidates.allowed, :]
    return DataFrame(map(0:((1 << nrow(allowed)) - 1)) do mask
        indices = [index for index in 1:nrow(allowed) if !iszero(mask & (1 << (index - 1)))]
        edges = allowed[indices, :]
        endpoints = vcat(collect(zip(edges.a_spacecraft, edges.a_terminal)),
            collect(zip(edges.b_spacecraft, edges.b_terminal)))
        terminal_feasible = allunique(endpoints)
        geometric = all(edges.geometrically_feasible)
        positive = all(>(0), edges.candidate_score_m_per_s)
        (; mask, matching=isempty(indices) ? "empty" : join(edges.pair, ";"),
            terminal_feasible, geometrically_feasible=geometric, positive_scores=positive,
            feasible=terminal_feasible && geometric && positive,
            total_score_m_per_s=sum(edges.candidate_score_m_per_s; init=0.0))
    end)
end

function write_prediction(args, candidates, matchings, output)
    mkpath(output)
    spacecraft = DataFrame(index=1:4, role=["helper", "helper", "target", "target"],
        mass_kg=collect(MASSES_KG), altitude_m=collect(ALTITUDES_M),
        phase_deg=collect(PHASES_DEG),
        terminal_capacity=terminal_count.(args.dynamics_model.spacecraft))
    CSV.write(joinpath(output, "spacecraft_before_run.csv"), spacecraft)
    CSV.write(joinpath(output, "candidates_before_run.csv"), candidates)
    CSV.write(joinpath(output, "matchings_before_run.csv"), matchings)
    println("BEFORE RUN: path 1-3-2-4; predicted matching ", matching_label(PREDICTED_PAIRS))
    show(stdout, MIME"text/plain"(), candidates; allrows=true, allcols=true)
    println()
    show(stdout, MIME"text/plain"(), matchings[matchings.feasible, :]; allrows=true, allcols=true)
    println()
    flush(stdout)
end

function check_history(args, solution, history; fixed_prediction::Bool=true)
    rows = NamedTuple[]
    states = map(history) do sample
        index = searchsortedlast(solution.t, sample.time)
        @test index > 0 && solution.t[index] == sample.time
        solution.u[index]
    end
    for (sample, state) in zip(history, states)
        candidates = candidate_table(args, state)
        matchings = enumerate_matchings(candidates)
        feasible = matchings[matchings.feasible, :]
        best = maximum(feasible.total_score_m_per_s)
        expected_available = Set(link_key((row.a_spacecraft, row.b_spacecraft))
            for row in eachrow(candidates) if row.available)
        @test Set(sample.available) == expected_available
        if fixed_prediction
            @test Set(sample.active) == Set(link_key.(PREDICTED_PAIRS))
        end
        @test all(key -> key in expected_available, sample.active)
        endpoints = collect(Iterators.flatten(sample.active))
        @test allunique(endpoints)
        terminal_use = [count(endpoint -> endpoint[1] == i, endpoints) for i in 1:4]
        @test all(terminal_use .<= terminal_count.(args.dynamics_model.spacecraft))
        @test length(sample.scores) == length(sample.active)
        for (key, recorded_score) in zip(sample.active, sample.scores)
            @test recorded_score > 0
            row = only(eachrow(candidates[(candidates.a_spacecraft .== key[1][1]) .&
                (candidates.b_spacecraft .== key[2][1]), :]))
            @test isapprox(recorded_score, row.candidate_score_m_per_s; rtol=1e-12, atol=1e-14)
        end
        selected_total = sum(sample.scores)
        @test isapprox(selected_total, best; rtol=1e-12, atol=1e-14)
        excluded = only(eachrow(candidates[candidates.pair .== "1-4", :]))
        if fixed_prediction
            @test excluded.geometrically_feasible && excluded.hypothetical_score_m_per_s > 0
        end
        @test !(link_key((1, 4)) in sample.available)
        @test !(link_key((1, 4)) in sample.active)
        largest = maximum(candidates[candidates.allowed, :candidate_score_m_per_s])
        if fixed_prediction
            @test best > largest
        end
        pair_geometry = (; (Symbol("geometrically_feasible_$(row.a_spacecraft)_$(row.b_spacecraft)") =>
            row.geometrically_feasible for row in eachrow(candidates))...)
        push!(rows, (; time_s=sample.time,
            available=matching_label([(key[1][1], key[2][1]) for key in sample.available]),
            active=matching_label([(key[1][1], key[2][1]) for key in sample.active]),
            selected_total_m_per_s=selected_total, exhaustive_best_m_per_s=best,
            largest_edge_m_per_s=largest, advantage_over_largest_m_per_s=best - largest,
            maximum_terminal_use=maximum(terminal_use),
            excluded_1_4_distance_m=excluded.distance_m,
            excluded_1_4_geometrically_feasible=excluded.geometrically_feasible,
            excluded_1_4_active=link_key((1, 4)) in sample.active, pair_geometry...))
    end
    return DataFrame(rows), states
end

function check_endpoint_dv(history, states; fixed_prediction::Bool=true)
    impulse = [zeros(3) for _ in 1:4]
    active_duration = zeros(4)
    used_pairs = Set{Tuple{Int, Int}}()
    # Trapezoidal force integration over accepted steps, independent of the ODE's
    # interlink_dv state. The preceding sample's selection owns each interval.
    for index in 1:(length(history) - 1)
        dt = history[index + 1].time - history[index].time
        @test 0 < dt <= 1.0 + 1e-12
        active_endpoints = falses(4)
        for key in history[index].active
            pair = (key[1][1], key[2][1])
            push!(used_pairs, pair)
            active_endpoints[pair[1]] = active_endpoints[pair[2]] = true
            integral = dt / 2 * (force_on_first(states[index], pair) +
                force_on_first(states[index + 1], pair))
            impulse[pair[1]] .+= integral
            impulse[pair[2]] .-= integral
        end
        for endpoint in 1:4
            if active_endpoints[endpoint]
                active_duration[endpoint] += dt
            else
                @test states[index + 1].sc[endpoint].interlink_dv == states[index].sc[endpoint].interlink_dv
            end
        end
    end
    for (i, j) in used_pairs
        @test MASSES_KG[i] != MASSES_KG[j]
    end
    recorded_dv = [states[end].sc[index].interlink_dv - states[1].sc[index].interlink_dv for index in 1:4]
    rows = map(1:4) do index
        @test all(state -> state.sc[index].mass == MASSES_KG[index], states)
        expected = impulse[index] / MASSES_KG[index]
        recorded = recorded_dv[index]
        error = norm(recorded - expected)
        tolerance = 1e-12 + 2e-6 * norm(expected)
        @test error <= tolerance
        (; spacecraft=index, mass_kg=MASSES_KG[index], active_duration_s=active_duration[index],
            force_magnitude_n=FORCE_N,
            acceleration_magnitude_m_per_s2=FORCE_N / MASSES_KG[index],
            expected_dvx_m_per_s=expected[1], expected_dvy_m_per_s=expected[2],
            expected_dvz_m_per_s=expected[3],
            recorded_dvx_m_per_s=recorded[1], recorded_dvy_m_per_s=recorded[2],
            recorded_dvz_m_per_s=recorded[3], error_norm_m_per_s=error,
            tolerance_m_per_s=tolerance)
    end
    @test isapprox(sum(MASSES_KG[index] * recorded_dv[index] for index in 1:4),
        zeros(3); atol=1e-12, rtol=0)
    # After partner changes, accumulated momentum cancels across the constellation,
    # not necessarily within each initially selected pair.
    if fixed_prediction
        for (i, j) in PREDICTED_PAIRS
            @test isapprox(MASSES_KG[i] * recorded_dv[i] + MASSES_KG[j] * recorded_dv[j],
                zeros(3); atol=1e-12, rtol=0)
        end
    end
    return DataFrame(rows)
end

function write_checkin(args, candidates, matchings, schedule, comparison, dv, solution, output)
    revision = readchomp(`git -C $REPO_ROOT rev-parse HEAD`)
    open(joinpath(output, "checkin.md"), "w") do io
        println(io, "# Known-answer scheduler check\n\nReviewer: Nihal. Due: 2 November 2026.\n")
        println(io, "## Prediction, written before propagation\n")
        println(io, "Adapted from the forty-spacecraft example: two 1000 km helpers (1, 2), ",
            "two 1100 km targets (3, 4), circular equatorial orbits, inverse-square gravity, ",
            "no atmosphere. Phases are deliberately clustered, not the original 40-satellite spacing.")
        println(io, "Masses: $(MASSES_KG) kg; initial phases: $(PHASES_DEG) degrees. ",
            "Each spacecraft has **one virtual terminal**. Allowed graph: **1-3-2-4**; ",
            "only (1,3), (2,3), (2,4) are registered. No additional forbidden-pair filter is needed.\n")
        println(io, "The link is LaserMomentumLink: 10,000 W, magnification 100, range 200,000 m, ",
            "equal/opposite force magnitude $FORCE_N N. Its geometric test is `0 < distance <= range`; ",
            "it does not model pointing, occultation, beam loss or hardware power shedding.\n")
        println(io, "Policy: `SemimajorAxisRatePolicy(targets=3:4)`. Independently computed score: ",
            "`sum(2*a^2*dot(v,F)/(mu*m))` over target endpoints, in m/s of semimajor-axis rate ",
            "(not velocity change). Earth mu = $(args.environment_model.planet.μ) m^3/s^2.\n")
        println(io, "| Pair | Allowed | Distance (km) | Geometrically feasible | Candidate score (m/s) |")
        println(io, "|---|---|---:|---|---:|")
        for row in eachrow(candidates)
            score = ismissing(row.candidate_score_m_per_s) ? "not a candidate" :
                @sprintf("%.12g", row.candidate_score_m_per_s)
            @printf(io, "| %s | %s | %.6f | %s | %s |\n",
                row.pair, string(row.allowed), row.distance_m / 1000,
                string(row.geometrically_feasible), score)
        end
        println(io, "\nAll five feasible matchings (the other three subsets violate a terminal capacity):\n")
        println(io, "| Matching | Initial total score (m/s) |\n|---|---:|")
        for row in eachrow(matchings[matchings.feasible, :])
            @printf(io, "| %s | %.12g |\n", row.matching, row.total_score_m_per_s)
        end
        println(io, "\n**Predicted matching: (1,3) + (2,4)**. Edge (2,3) is the individually ",
            "largest, but blocks both outer edges; their sum is larger. The independent enumeration ",
            "is saved before `run_simulation` is called.\n")
        println(io, "## Recorded windows and comparison\n")
        println(io, "Ran $(DURATION_S) s with Tsit5, dtmax 1 s, reltol 1e-10, abstol 1e-12; ",
            "$(nrow(comparison)) scheduler samples including initialization. ",
            "Solver return code: $(solution.retcode). Schedule history is automatic; ",
            "`record_available=true`, `results=true` and `save_csv=true` also save contact windows ",
            "and the schedule CSV. History was read from the returned run, not the input configuration.\n")
        println(io, "| Window | Pair | Start (s) | End (s) | Duration (s) | Samples |")
        println(io, "|---|---|---:|---:|---:|---:|")
        for row in eachrow(schedule)
            @printf(io, "| %s | %d-%d | %.6f | %.6f | %.6f | %d |\n",
                row.window, row.a_spacecraft, row.b_spacecraft, row.start_s,
                row.end_s, row.duration_s, row.steps)
        end
        println(io, "\n**No scheduler disagreement.** At every recorded scheduling instant, ",
            "the selected set equals the prediction, its recorded total agrees with exhaustive ",
            "enumeration, and no terminal is used twice. Edge (2,3) is allowed and geometrically ",
            "feasible but idle because the two outer edges give the better total.")
        @printf(io, " The minimum advantage over the largest single edge is %.12g m/s.\n",
            minimum(comparison.advantage_over_largest_m_per_s))
        @printf(io, "\nExcluded pair (1,4) stays geometrically feasible at all %d samples ",
            nrow(comparison))
        @printf(io, "(distance %.6f to %.6f km) and would have a positive score, ",
            minimum(comparison.excluded_1_4_distance_m) / 1000,
            maximum(comparison.excluded_1_4_distance_m) / 1000)
        println(io, "but is never available to the scheduler or active because it was not registered. ",
            "All six pairs happen to be geometrically feasible in this clustered case; ",
            "only the three allowed pairs are scheduling candidates. Windows have accepted-step ",
            "resolution, not continuous event-location precision.\n")
        println(io, "## Unequal-mass acceleration check\n")
        println(io, "For each endpoint, compare the recorded `interlink_dv` with the independently ",
            "trapezoid-integrated force impulse divided by **that endpoint's own mass**. ",
            "This is the link-only inertial velocity contribution, not total velocity change ",
            "(which includes gravity), nor a subtraction of a separately propagated orbit. ",
            "The 1 s sampling uses tolerance `1e-12 + 2e-6*norm(expected_dv)` m/s; ",
            "small residuals are quadrature error, not matching disagreements. ",
            "Each selected pair also passes equal/opposite mass-weighted impulse conservation.\n")
        println(io, "| Spacecraft | Mass (kg) | F/m (m/s^2) | DV error norm (m/s) | Allowed error (m/s) |")
        println(io, "|---|---:|---:|---:|---:|")
        for row in eachrow(dv)
            @printf(io, "| %d | %.1f | %.12g | %.12g | %.12g |\n",
                row.spacecraft, row.mass_kg, row.acceleration_magnitude_m_per_s2,
                row.error_norm_m_per_s, row.tolerance_m_per_s)
        end
        println(io, "\n## Reproduce and inspect\n")
        println(io, "From the repository root:\n\n```sh\n",
            "julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'\n",
            "julia --startup-file=no --project=. examples/known_answer_tree_interlinks.jl\n```\n")
        println(io, "Verified with Julia $VERSION, SpaceAGORA $(pkgversion(SpaceAGORA)), ",
            "checkout `$revision`. No scheduler or force implementation was changed.\n")
        println(io, "- [Spacecraft and capacities](spacecraft_before_run.csv)\n",
            "- [Separate pre-run candidate list, including excluded pairs](candidates_before_run.csv)\n",
            "- [All eight subsets and feasibility](matchings_before_run.csv)\n",
            "- [Recorded link schedule](simulation_results_link_schedule.csv)\n",
            "- [Every accepted-step comparison](accepted_step_comparison.csv)\n",
            "- [Endpoint velocity vectors and mass checks](endpoint_dv_check.csv)\n")
        println(io, "Generated outputs are local and git-ignored; rerun the script or attach these files ",
            "when sharing the check-in. The schedule lists windows and selected scores, not a ",
            "complete rejection explanation; the separate candidate and matching tables provide that.\n")
        println(io, "References: [interlink documentation](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/",
            REFERENCE, "/docs/src/user/interlinks.md), ",
            "[forty-spacecraft example](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/",
            REFERENCE, "/examples/targeted_40sat_sma_interlinks.jl). ",
            "The reference example is unchanged in this checkout; current documentation has expanded.\n")
        println(io, "## Separate follow-up: cycles\n\n",
            "This first check uses only Julia's exact tree matcher. PythonCall/NetworkX was not ",
            "installed or loaded for it. A cyclic candidate graph belongs in a separate environment ",
            "with the documented PythonCall extension and NetworkX >=3.4,<4; it is not part of this check.")
    end
end

function run_check(; output=OUTPUT)
    args = configuration(; output)
    initial = SpaceAGORA.SimulationEngine.build_initial_conditions(args)
    candidates = candidate_table(args, initial)
    matchings = enumerate_matchings(candidates)
    write_prediction(args, candidates, matchings, output)
    evidence = nothing
    @testset "Known-answer tree interlinks" begin
        @test Set(keys(args.interlink_model.graph.edge_index_lookup)) == Set(link_key.(ALLOWED_PAIRS))
        @test terminal_count.(args.dynamics_model.spacecraft) == ones(Int, 4)
        @test nrow(matchings) == 8
        feasible = matchings[matchings.feasible, :]
        @test nrow(feasible) == 5
        @test feasible.matching[argmax(feasible.total_score_m_per_s)] == matching_label(PREDICTED_PAIRS)
        allowed = candidates[candidates.allowed, :]
        @test allowed.pair[argmax(allowed.candidate_score_m_per_s)] == "2-3"
        @test maximum(feasible.total_score_m_per_s) > maximum(allowed.candidate_score_m_per_s)
        @test Base.get_extension(SpaceAGORA, :SpaceAGORAPythonCallExt) === nothing
        solution = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "true",
            "SPACEAGORA_RHS_CALIBRATE" => "off", "SPACEAGORA_RHS_IDENTIFY" => "0",
            "SPACEAGORA_RHS_EXECUTION_MODE" => "serial", "SPACEAGORA_VISUALIZATION" => "0") do
            run_simulation(args; return_solution=true)
        end
        @test successful_retcode(solution)
        @test solution.t[end] == DURATION_S
        @test all(isfinite, solution.u[end])
        run_args = solution.prob.p.args
        history = run_args.interlink_model.history
        @test isempty(args.interlink_model.history)
        @test length(history) == solution.stats.naccept + 1
        @test first(history).time == 0.0 && last(history).time == DURATION_S
        comparison, states = check_history(run_args, solution, history)
        dv = check_endpoint_dv(history, states)
        schedule = link_schedule(run_args.interlink_model; policy=run_args.scheduling_policy_model)
        @test nrow(schedule) == 5
        @test count(==("active"), schedule.window) == 2
        @test count(==("contact"), schedule.window) == 3
        @test all(==(0.0), schedule.start_s)
        @test all(==(DURATION_S), schedule.end_s)
        @test all(==(DURATION_S), schedule.duration_s)
        @test Base.get_extension(SpaceAGORA, :SpaceAGORAPythonCallExt) === nothing
        recorded = CSV.read(joinpath(output, "simulation_results_link_schedule.csv"), DataFrame)
        @test isequal(recorded, schedule)
        CSV.write(joinpath(output, "accepted_step_comparison.csv"), comparison)
        CSV.write(joinpath(output, "endpoint_dv_check.csv"), dv)
        evidence = (; schedule, comparison, dv, solution)
    end
    write_checkin(args, candidates, matchings, evidence.schedule, evidence.comparison,
        evidence.dv, evidence.solution, output)
    println("Check-in and verified evidence saved to ", output)
    return output
end

end # module KnownAnswerTreeInterlinks

if abspath(PROGRAM_FILE) == @__FILE__
    KnownAnswerTreeInterlinks.run_check()
end
