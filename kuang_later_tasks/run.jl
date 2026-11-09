module KuangLaterTasks

using SpaceAGORA
using SpaceAGORA.SimulationModel
using CSV
using DataFrames
using Dates
using LinearAlgebra
using Printf
using SHA
using TOML
using Test
using SciMLBase: successful_retcode

const ROOT = realpath(joinpath(@__DIR__, ".."))
const OUTPUT = joinpath(@__DIR__, "results")
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PLANET = make_no_gram_planet(:earth)
const LinkKey = Tuple{Tuple{Int, Int}, Tuple{Int, Int}}
const BIT_ATOL = 1.0
const SCORE_RTOL = 1e-12
const DATA_SCORE_ATOL = 1e-5
const MOMENTUM_SCORE_ATOL = 1e-14

key_label(key) = "$(key[1][1]).$(key[1][2])->$(key[2][1]).$(key[2][2])"
matching_label(keys) = isempty(keys) ? "empty" : join(key_label.(sort(collect(keys))), ";")
row_key(row) = ((row.tx, row.tx_terminal), (row.rx, row.rx_terminal))

function write_toml(path, values)
    open(path, "w") do io
        TOML.print(io, values; sorted=true)
    end
end

function configuration(spacecraft, model, policy; duration, dtmax, output, data_rate=1.0)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, generate_plots=false,
            results_directory=output, save_csv=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=duration, orientation_sim=false, data_rate=data_rate),
        environment_model=make_no_gram_environment(planet=deepcopy(PLANET)),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10,
            abstol_orbit=1e-8, dt_max_orbit=dtmax),
        solver_config=SM.SolverConfig(solver_mode=:tsit5),
        interlink_model=model,
        scheduling_policy_model=policy)
end

function run_checked(args)
    solution = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "true",
        "SPACEAGORA_RHS_CALIBRATE" => "off", "SPACEAGORA_RHS_IDENTIFY" => "0",
        "SPACEAGORA_RHS_EXECUTION_MODE" => "serial", "SPACEAGORA_VISUALIZATION" => "0") do
        run_simulation(args; return_solution=true)
    end
    @test successful_retcode(solution)
    @test solution.t[end] == args.mission_configuration.mission_time
    @test all(state -> all(isfinite, state), solution.u)
    @test maximum(diff(solution.t)) <= args.integration_tolerances.dt_max_orbit + 1e-9
    return solution
end

function history_states(solution)
    history = solution.prob.p.args.interlink_model.history
    @test !isempty(history)
    @test first(history).time == 0.0
    @test last(history).time == solution.t[end]
    @test all(>(0.0), diff([sample.time for sample in history]))
    return map(history) do sample
        index = searchsortedlast(solution.t, sample.time)
        @test index > 0 && solution.t[index] == sample.time
        solution.u[index]
    end
end

function check_matching(sample, allowed, spacecraft)
    @test issubset(Set(sample.available), allowed)
    @test issubset(Set(sample.active), Set(sample.available))
    @test allunique(sample.active)
    endpoints = collect(Iterators.flatten(sample.active))
    @test allunique(endpoints)
    @test all(i -> count(endpoint -> endpoint[1] == i, endpoints) <=
        terminal_count(spacecraft[i]), eachindex(spacecraft))
    @test length(sample.active) == length(sample.scores)
    @test all(>(0.0), sample.scores)
end

function save_schedule(solution, directory)
    args = solution.prob.p.args
    schedule = link_schedule(args.interlink_model; policy=args.scheduling_policy_model)
    path = joinpath(directory, "link_schedule.csv")
    CSV.write(path, schedule)
    @test isequal(CSV.read(path, DataFrame), schedule)
    return schedule
end

function schedule_changes(history)
    return count(i -> Set(history[i].active) != Set(history[i - 1].active), 2:length(history))
end

include("constellation.jl")
include("backlog.jl")

function write_report(constellation, backlog, provenance, elapsed)
    open(joinpath(@__DIR__, "REPORT.md"), "w") do io
        println(io, "# Kuang: later-task results\n")
        println(io, "All assertions passed. Reviewer: Nihal. Generated: ", provenance["started_utc"],
            ". Runtime: ", round(elapsed; digits=2), " s.\n")
        println(io, "## Shortened forty-spacecraft momentum run\n")
        println(io, "Adapted from the existing forty-spacecraft example: 30 helpers at 1000 km, ",
            "10 targets at 1100 km, circular equatorial orbits, original ring spacing and ",
            "1 degree target phase offset. All masses remain 227 kg and each spacecraft has ",
            "one virtual terminal. Allowed pairs are all 300 helper-to-target pairs; other ",
            "pairs are not registered. This run does **not** replace the earlier unequal-mass check.\n")
        @printf(io, "- Duration: %.9f s (0.25 target orbits, rather than 20).\n", constellation.duration_s)
        @printf(io, "- Recorded scheduling samples: %d; selected-set changes: %d; maximum active links: %d.\n",
            constellation.samples, constellation.switches, constellation.maximum_active)
        @printf(io, "- Target SMA change relative to a separately propagated no-link control: %.9g to %.9g m.\n",
            constellation.minimum_target_effect_m, constellation.maximum_target_effect_m)
        @printf(io, "- Maximum no-link SMA drift across all 40 spacecraft: %.9g m (limit 0.02 m).\n\n",
            constellation.maximum_control_drift_m)
        println(io, "`LaserMomentumLink`: 10,000 W, magnification 100, range 200 km; ",
            "ideal equal/opposite force `100 * 10000 / 299792458` N. ",
            "`SemimajorAxisRatePolicy(targets=31:40)` is checked against the independent ",
            "energy-rate formula `2*a^2*dot(v,F)/(mu*m)` at every recorded state. ",
            "All feasible terminal components were verified to be disjoint edges, so selecting ",
            "all positive-score edges is the exact optimum even though the full allowed graph is cyclic. ",
            "No Python backend was loaded.\n")
        println(io, "Attitude, animation, plotting and atmosphere are disabled. The ideal momentum ",
            "model uses range-only feasibility: no mounted-beam capture, pointing or occultation ",
            "is claimed here. Gaussian capture and pointing loss exist in the separate beam model; ",
            "their ORACLE adequacy belongs with Nihal's physics work.\n")
        println(io, "Evidence: [settings](results/constellation/settings.toml), ",
            "[pre-run candidates](results/constellation/candidates_before_run.csv), ",
            "[every candidate at every accepted step](results/constellation/candidate_scores.csv), ",
            "[matching checks](results/constellation/accepted_steps.csv), ",
            "[recorded contact/active windows](results/constellation/link_schedule.csv), ",
            "[SMA comparison](results/constellation/spacecraft_comparison.csv), ",
            "[linked states](results/constellation/states.csv), ",
            "[no-link states](results/constellation/control_states.csv).\n")
        println(io, "## Separate data-link experiments\n")
        println(io, "Circular co-orbital spacecraft at 700 km, separated by 10 km of arc. ",
            "All endpoints have finite data stores; the final spacecraft is a finite-store ",
            "**destination**, not a ground station. Directed `DataLink`s use `FixedRateBudget`, ",
            "virtual terminals and zero transmit power to isolate scheduling from RF/power physics. ",
            "Data-link planet occultation is enabled by the built-in model and independently checked. ",
            "All runs last 20 s with backlog horizon 5 s.\n")
        println(io, "| Case | dtmax (s) | Samples | Set changes | Delivered to destination (bit) | Production (bit) | Permitted production overflow (bit) | Max node balance residual (bit) |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|")
        for row in eachrow(backlog.summary)
            @printf(io, "| [%s](results/backlog/%s/accepted_steps.csv) | %.3f | %d | %d | %.9f | %.9f | %.9f | %.6g |\n",
                row.case, row.case, row.dtmax_s, row.samples, row.switches, row.delivered_bits,
                row.produced_bits, row.lost_bits, row.maximum_node_balance_error_bits)
        end
        println(io, "\n### What the fixtures establish\n")
        println(io, "1. **Empty producer control:** a 200,000 bit/s source remains selected and ",
            "delivers exactly 4,000,000 bit over 20 s, to tolerance.\n",
            "2. **Competition:** that same empty producer initially has a positive score, but ",
            "loses the destination's single terminal to another sender with 8,000,000 bit. ",
            "Later selections and residual backlog are recorded: ongoing production is not ",
            "a selection guarantee, and backlog score is not a fairness guarantee.\n",
            "3. **Busy relay:** source 1 starts with 80,000,000 bit; relay 2 starts full at ",
            "20,000,000 bit and produces 1,000,000 bit/s; destination 3 starts empty. ",
            "Rates are 1,000,000 bit/s for source-to-relay and relay-to-destination, and ",
            "500,000 bit/s direct. Terminal counts are (1,2,2). ",
            "The terminal graph is a forest: relay input and output use distinct terminals, ",
            "as do the destination's direct and relay inputs.\n")
        println(io, "Initially the full relay has no **previous** outflow credit, so the scheduler ",
            "selects direct + relay-out. At the next accepted step, its previously active ",
            "feasible relay-out supplies 1,000,000 bit/s credit; relay-in + relay-out scores ",
            "2,000,000 bit/s, beating direct + relay-out at 1,500,000 bit/s. ",
            "Actual delivery then falls to 1,000,000 bit/s, while the relay's own production ",
            "overflows at 1,000,000 bit/s. Incoming traffic displaces local production; ",
            "a full receiver not overfilling is not a no-loss proof.\n")
        @printf(io, "The same-backlog-policy direct-route control explicitly disables relay-in, ",
            )
        @printf(io, "delivering %.9f bit more and losing %.9f bit less than the busy-relay run. ",
            backlog.delivery_advantage_bits, backlog.loss_reduction_bits)
        println(io, "This is a feasible routing comparator, **not** a claim of a globally optimal ",
            "routing policy. Store counters do not track packet origin, so destination totals ",
            "must not be described as source-1-only delivery.\n")
        println(io, "Scores are checked independently as ",
            "`min(rate, held_tx/5 + production_tx, room_rx/5 + previous_feasible_outflow_rx)`. ",
            "All matching subsets are enumerated at every scheduling sample. A read-only ",
            "scorer replay also checks rejected candidates, using the previous selection and ",
            "store flags reconstructed at that state, never the newly selected links. ",
            "The separate [outflow probes](results/backlog/outflow_probes.csv) demonstrate ",
            "zero credit with no previous link, positive credit with a previous feasible link, ",
            "and zero credit when that previous link is out of range.\n")
        println(io, "Every case saves its spacecraft/stores, pre-run candidates and matching subsets, ",
            "full candidate-score audit, selected-set history, contact/active windows, ",
            "per-node accounting and network accounting. ",
            "[Analytic comparisons](results/backlog/hand_comparisons.csv) integrate the simple ",
            "constant-rate modes over the **recorded accepted-step** intervals, not a fictional ",
            "continuously located scheduling event. ",
            "[Step-refinement results](results/backlog/step_refinement.csv) compare dtmax ",
            "0.25 s with 0.125 s for competition and relay selection. Numerical solver error ",
            "and scheduler time-resolution sensitivity are different quantities.\n")
        println(io, "## Accounting, frames and tolerances\n")
        println(io, "- Positions/velocities are planet-centered inertial Cartesian SI values; ",
            "initial angles are degrees. No attitude or body torques are propagated.\n",
            "- Solver: Tsit5; reltol 1e-10, abstol 1e-8; dtmax 10 s for the constellation, ",
            "0.25 s for data cases, 0.125 s for refined cases. Every accepted state is retained.\n",
            "- Momentum scores: rtol 1e-12, atol 1e-14 m/s. Data scores: rtol 1e-12, ",
            "atol 1e-5 bit/s. Matching ties may select any equal-optimal set.\n",
            "- Data counter/hand-comparison absolute tolerance: 1 bit. Store bounds and ",
            "monotonic counters are checked at every saved accepted state. Store scoring ",
            "uses the engine's documented full/empty threshold of 1e-9 times capacity.\n",
            "- Per node: `stored = initial + produced + received - sent - lost`. ",
            "Internal sent and received totals must agree. Network accounting includes ",
            "the destination store, so delivery is **not subtracted again**: ",
            "`sum(initial) + sum(produced) = sum(stored) + sum(lost)`.\n",
            "- Loss is permitted only at relay 2 in the busy-relay cases, with its quantity ",
            "independently predicted. Other fixtures require zero loss to 1 bit. ",
            "These bounded checks make no global no-loss claim.\n",
            "- `record_available=true`; the saved schedule is the engine's history converted ",
            "by `link_schedule`, not a schedule inferred from raw force/rate queries. ",
            "Windows and switches have accepted-step resolution.\n")
        println(io, "## Coordination and separate work\n")
        println(io, "Give Nihal the per-node/network accounting and hand-comparison tables for ",
            "review of production, transfers and permitted overflow; Kuang owns the matching ",
            "and delivery comparison. Ground-station elevation/planet geometry, viewer/animation ",
            "validation, a genuinely cyclic feasible-graph backend test, mounted-beam physics ",
            "and resumed laser history are **not tested** here. Resumed laser history remains ",
            "Nihal's later task. Fixed rates, no acquisition time and no hardware power budget ",
            "are deliberate scheduler fixtures, not a validated ORACLE communication model.\n")
        println(io, "## Reproducibility\n")
        println(io, "Julia ", VERSION, "; SpaceAGORA ", pkgversion(SpaceAGORA),
            "; checkout `", provenance["git_revision"], "`. ",
            "[Provenance and hashes](results/provenance.toml), [run status](results/status.toml), ",
            "[data summary](results/backlog/summary.csv), [run instructions](README.md). ",
            "No existing repository files were edited and no dependencies were installed. ",
            "All generated deliverables are inside this folder.")
    end
end

function main()
    realpath(dirname(dirname(pathof(SpaceAGORA)))) == ROOT ||
        error("Use the parent repository project; the loaded SpaceAGORA is from another checkout.")
    mkpath(OUTPUT)
    started = time()
    provenance = Dict(
        "started_utc" => string(now(UTC)),
        "julia_version" => string(VERSION),
        "spaceagora_version" => string(pkgversion(SpaceAGORA)),
        "spaceagora_source" => pathof(SpaceAGORA),
        "git_revision" => readchomp(`git -C $ROOT rev-parse HEAD`),
        "git_status_before_run" => readchomp(`git -C $ROOT status --short`),
        "project_sha256" => bytes2hex(sha256(read(joinpath(ROOT, "Project.toml")))),
        "manifest_sha256" => bytes2hex(sha256(read(joinpath(ROOT, "Manifest.toml")))),
        "code_sha256" => Dict(name => bytes2hex(sha256(read(joinpath(@__DIR__, name))))
            for name in ("run.jl", "constellation.jl", "backlog.jl")))
    write_toml(joinpath(OUTPUT, "provenance.toml"), provenance)
    status = Dict{String, Any}("status" => "running", "started_utc" => provenance["started_utc"])
    write_toml(joinpath(OUTPUT, "status.toml"), status)
    passed = false
    try
        constellation = nothing
        backlog = nothing
        @testset "Kuang later-task validation" begin
            @test Base.get_extension(SpaceAGORA, :SpaceAGORAPythonCallExt) === nothing
            constellation = run_constellation()
            backlog = run_backlog_suite()
            @test Base.get_extension(SpaceAGORA, :SpaceAGORAPythonCallExt) === nothing
        end
        write_report(constellation, backlog, provenance, time() - started)
        passed = true
    finally
        status["status"] = passed ? "passed" : "failed"
        status["completed_utc"] = string(now(UTC))
        status["elapsed_s"] = time() - started
        write_toml(joinpath(OUTPUT, "status.toml"), status)
    end
    println("Verified results and REPORT.md saved in ", @__DIR__)
    return nothing
end

end # module KuangLaterTasks

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) || error("This runner takes no arguments; see README.md for focused reruns.")
    KuangLaterTasks.main()
end
