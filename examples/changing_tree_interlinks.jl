module ChangingTreeInterlinks

using SpaceAGORA
using CSV
using DataFrames
using Printf
using SciMLBase: successful_retcode
using Test
include(joinpath(@__DIR__, "known_answer_tree_interlinks.jl"))
const K = KnownAnswerTreeInterlinks
const OUTPUT = joinpath(K.REPO_ROOT, "output", "changing_tree_interlinks")
const DURATION_S = 2400.0

function status_changes(comparison)
    isempty(comparison) && throw(ArgumentError("Cannot find changes in an empty recording."))
    indices = [1; [index for index in 2:nrow(comparison)
        if comparison.active[index] != comparison.active[index - 1] ||
            comparison.available[index] != comparison.available[index - 1]]]
    return comparison[indices, :]
end

function write_summary(changes, schedule, dv, sample_count, output)
    open(joinpath(output, "status_changes.md"), "w") do io
        println(io, "# Changing-link demonstration\n")
        println(io, "Ran **2400 s (40 minutes)** with the original four spacecraft, initial phases, ",
            "unequal masses, one terminal per spacecraft and allowed tree `1-3-2-4`. ",
            "Only the duration changed; this example writes to a separate directory from the 60-second run.\n")
        println(io, "The 1000 km helpers orbit faster than the 1100 km targets. Their relative ",
            "phase drifts by approximately 0.00114107 degrees/s in the unforced circular model. ",
            "This changes the semimajor-axis-rate scores and later takes links beyond 200 km. ",
            "The table below is from the actual forced simulation, not that approximation.\n")
        println(io, "| First recorded time (s) | Allowed, available links | Active matching | Explanation |")
        println(io, "|---:|---|---|---|")
        for (index, row) in enumerate(eachrow(changes))
            reason = if index == 1
                "Initial maximum-score matching."
            elseif row.available != changes.available[index - 1]
                "Separations crossed the 200 km range limit."
            elseif row.active == "empty"
                "No available candidate has a positive target SMA-rate score."
            else
                "The single middle edge now beats the sum of the two outer edges."
            end
            @printf(io, "| %.6f | %s | %s | %s |\n",
                row.time_s, row.available, row.active, reason)
        end
        println(io, "\nEvery one of the $sample_count scheduling samples was checked against exhaustive ",
            "matching enumeration using independently calculated scores. All selected links are ",
            "allowed, geometrically available, positive-scoring and within the one-terminal limit. ",
            "All excluded pairs stay unregistered. Time resolution is the accepted solver step ",
            "(at most 1 s); these are recorded sample times, not precisely located crossing times.\n")
        println(io, "## Recorded windows\n")
        println(io, "| Window | Pair | Start (s) | End (s) | Duration (s) |")
        println(io, "|---|---|---:|---:|---:|")
        for row in eachrow(schedule)
            @printf(io, "| %s | %d-%d | %.6f | %.6f | %.6f |\n",
                row.window, row.a_spacecraft, row.b_spacecraft,
                row.start_s, row.end_s, row.duration_s)
        end
        println(io, "\n## Equatorial geometry and endpoint velocity check\n")
        println(io, "All four initial orbits have inclination 0 degrees and eccentricity 0: ",
            "helpers start at 1000 km and targets at 1100 km. All $sample_count accepted-state ",
            "samples remain in the equatorial plane, with zero out-of-plane position and velocity.\n")
        println(io, "For each endpoint, independent trapezoidal force integration over the recorded ",
            "active intervals, divided by its own mass, agrees with its link-only velocity change. ",
            "The matching at the start of each interval owns that interval, including partner changes. ",
            "Every idle interval was checked for zero additional link-only velocity change. ",
            "Accumulated interlink momentum balances across all four spacecraft, rather than ",
            "within the original pairs after their partners switch.\n")
        println(io, "| Spacecraft | Mass (kg) | Active time (s) | Velocity error (m/s) | Tolerance (m/s) |")
        println(io, "|---|---:|---:|---:|---:|")
        for row in eachrow(dv)
            @printf(io, "| %d | %.0f | %.6f | %.6e | %.6e |\n", row.spacecraft, row.mass_kg,
                row.active_duration_s, row.error_norm_m_per_s, row.tolerance_m_per_s)
        end
        println(io, "\nThe [endpoint velocity CSV](endpoint_dv_check.csv) contains all expected and ",
            "recorded vector components. Its force and acceleration magnitudes apply while linked, ",
            "not during idle periods. These checks cover 0-2400 s; they do not compare total orbital ",
            "velocity change, which also includes gravity.\n")
        println(io, "\n## Animation and evidence\n")
        println(io, "- [Changing graph animation](interlinks.gif)\n",
            "- [Allowed-link status-change times](status_changes.csv)\n",
            "- [Candidate scores and distances at each allowed-link change](candidates_at_changes.csv)\n",
            "- [Every accepted-step comparison](accepted_step_comparison.csv)\n",
            "- [Recorded active/contact windows](simulation_results_link_schedule.csv)\n",
            "- [40-minute endpoint velocity and own-mass checks](endpoint_dv_check.csv)\n")
        println(io, "The animation preserves the simplified square and right-side legend. ",
            "Orange means active, green means geometrically available but idle, solid gray means ",
            "not in the allowed graph but geometrically feasible, and dashed light gray means ",
            "any pair is geometrically unavailable, whether allowed or excluded. ",
            "The accepted-step comparison records geometry for all six pairs. ",
            "Excluded-pair range crossings change the animation but not the allowed-link status table. ",
            "A green link can have a nonpositive score: geometry alone does not make it worth selecting. ",
            "The 241 uniformly spaced frames play at 12 fps, compressing 40 minutes to about 20 seconds; ",
            "the CSVs retain all scheduling samples. Nodes are schematic, not moving orbital positions. ",
            "PythonCall/NetworkX is not needed.\n")
        println(io, "From the repository root:\n\n```sh\n",
            "julia --startup-file=no --project=. examples/changing_tree_interlinks.jl\n",
            "julia --startup-file=no --project=examples -e 'include(\"examples/animate_known_answer_tree_interlinks.jl\"); ",
            "KnownAnswerTreeAnimation.save_animation(input=\"output/changing_tree_interlinks\", frame_count=241)'\n```\n",
            "Prepare the plotting environment with `julia --startup-file=no examples/setup.jl` if needed. ",
            "Generated outputs are local and git-ignored.")
    end
end

function run_example(; output=OUTPUT)
    args = K.configuration(; output, duration=DURATION_S)
    initial = SpaceAGORA.SimulationEngine.build_initial_conditions(args)
    candidates = K.candidate_table(args, initial)
    K.write_prediction(args, candidates, K.enumerate_matchings(candidates), output)
    println("Long-run prediction: outer edges -> middle edge -> no active links -> range losses.")
    evidence = nothing
    @testset "Changing tree interlinks (40 minutes)" begin
        for (index, spacecraft) in enumerate(args.dynamics_model.spacecraft)
            @test spacecraft.initial_condition.i == 0.0
            @test spacecraft.initial_condition.e == 0.0
            @test spacecraft.initial_condition.a - args.environment_model.planet.Rp_e == K.ALTITUDES_M[index]
        end
        solution = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "true",
            "SPACEAGORA_RHS_CALIBRATE" => "off", "SPACEAGORA_RHS_IDENTIFY" => "0",
            "SPACEAGORA_RHS_EXECUTION_MODE" => "serial", "SPACEAGORA_VISUALIZATION" => "0") do
            run_simulation(args; return_solution=true)
        end
        @test successful_retcode(solution)
        @test solution.t[end] == DURATION_S
        @test all(isfinite, solution.u[end])
        @test isempty(args.interlink_model.history)
        run_args = solution.prob.p.args
        history = run_args.interlink_model.history
        times = getproperty.(history, :time)
        @test length(history) == solution.stats.naccept + 1
        @test first(times) == 0.0 && last(times) == DURATION_S
        @test all(dt -> 0 < dt <= 1.0 + 1e-12, diff(times))
        @test Base.get_extension(SpaceAGORA, :SpaceAGORAPythonCallExt) === nothing
        comparison, states = K.check_history(run_args, solution, history; fixed_prediction=false)
        maximum_z_m = maximum(abs(state.sc[index].pos[3]) for state in states for index in 1:4)
        maximum_vz_m_per_s = maximum(abs(state.sc[index].vel[3]) for state in states for index in 1:4)
        @test maximum_z_m == 0.0
        @test maximum_vz_m_per_s == 0.0
        dv = K.check_endpoint_dv(history, states; fixed_prediction=false)
        @test Set(comparison.active) == Set(["1-3;2-4", "2-3", "empty"])
        @test first(comparison.active) == "1-3;2-4"
        @test last(comparison.active) == last(comparison.available) == "empty"
        @test length(unique(comparison.available)) >= 3
        changes = status_changes(comparison)
        @test nrow(changes) >= 5
        schedule = link_schedule(run_args.interlink_model; policy=run_args.scheduling_policy_model)
        @test count(==("active"), schedule.window) == 3
        @test count(==("contact"), schedule.window) == 3
        @test isequal(CSV.read(joinpath(output, "simulation_results_link_schedule.csv"), DataFrame), schedule)
        CSV.write(joinpath(output, "accepted_step_comparison.csv"), comparison)
        CSV.write(joinpath(output, "endpoint_dv_check.csv"), dv)
        CSV.write(joinpath(output, "status_changes.csv"), changes)
        event_candidates = map(changes.time_s) do time
            index = searchsortedfirst(times, time)
            table = K.candidate_table(run_args, states[index])
            rename!(table, :predicted_active => :initially_predicted_active)
            insertcols!(table, 1, :time_s => fill(time, nrow(table)))
            table.selected_active = [K.link_key((row.a_spacecraft, row.b_spacecraft)) in history[index].active
                for row in eachrow(table)]
            table
        end
        CSV.write(joinpath(output, "candidates_at_changes.csv"), vcat(event_candidates...))
        evidence = (; changes, schedule, dv, sample_count=nrow(comparison), maximum_z_m, maximum_vz_m_per_s)
    end
    write_summary(evidence.changes, evidence.schedule, evidence.dv, evidence.sample_count, output)
    show(stdout, MIME"text/plain"(), evidence.changes[:, [:time_s, :available, :active]]; allrows=true)
    println("\nEquatorial geometry: maximum |z| = $(evidence.maximum_z_m) m; ",
        "maximum |vz| = $(evidence.maximum_vz_m_per_s) m/s.")
    show(stdout, MIME"text/plain"(), evidence.dv[:, [:spacecraft, :mass_kg, :active_duration_s,
        :error_norm_m_per_s, :tolerance_m_per_s]]; allrows=true)
    println("\nVerified changing-link evidence saved to ", output)
    return output
end

end # module ChangingTreeInterlinks

if abspath(PROGRAM_FILE) == @__FILE__
    ChangingTreeInterlinks.run_example()
end
