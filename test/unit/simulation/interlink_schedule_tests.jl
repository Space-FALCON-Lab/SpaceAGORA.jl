module InterlinkScheduleTests

using Test
using Arrow
using CSV
using DataFrames
using DiffEqBase
using SHA
using TOML
using SpaceAGORA
include(joinpath(@__DIR__, "..", "..", "..", "examples", "gve_sma_interlinks.jl"))

const SE = SpaceAGORA.SimulationEngine
const SP = SM.SchedulingPolicies

function with_settings(args; overrides...)
    settings = args.simulation_settings
    names = fieldnames(typeof(settings))
    fields = NamedTuple{names}(map(name -> getfield(settings, name), names))
    return SM.SimConfig._with_configuration(args;
        simulation_settings=SimulationSettings(; merge(fields, (; overrides...))...))
end
with_policy(args, policy; model=args.interlink_model) =
    SM.SimConfig._with_configuration(args; scheduling_policy_model=SchedulingPolicyModel(policy), interlink_model=model)

# The gve example's two candidate laser links, on a model of one's own.
function laser_model(spacecraft; kwargs...)
    model = InterlinkModel(; kwargs...)
    for partner in 2:3
        register_link!(model, spacecraft, (1, 1), (partner, 1),
            LaserMomentumLink(power_w=10_000.0, magnification=100.0, range_m=200e3))
    end
    return model
end

# Weights and scales that are not powers of two, so that any change in the order of
# the mix's arithmetic changes its bits.
const MIX = WeightedPolicy(PolicyTerm(SemimajorAxisRatePolicy(targets=1); scale=0.3, weight=0.7),
    PolicyTerm(SemimajorAxisRatePolicy(targets=:both); scale=1.7, weight=-0.35))

bits(x) = reinterpret(UInt64, collect(Float64, x))

@testset "Interlink link schedule" begin
    @testset "history samples" begin
        key = ((1, 1), (2, 1))
        @test InterLinkHistory(1.0, [key]).scores == Float64[]
        @test size(InterLinkHistory(1.0, [key]).terms) == (0, 1)
        @test size(InterLinkHistory(1.0, [key], InterLinkKey[], [2.0]).terms) == (0, 1)
        @test InterLinkHistory(1.0, [key], InterLinkKey[], [2.0], [3.0 ; 4.0 ;;]).terms == [3.0; 4.0;;]
        @test_throws ArgumentError InterLinkHistory(1.0, [key], InterLinkKey[], [1.0, 2.0])
        @test_throws ArgumentError InterLinkHistory(1.0, [key], InterLinkKey[], [1.0], zeros(2, 2))
    end

    @testset "each sample records the score that selected each link" begin
        args = gve_sma_configuration(duration=600.0)
        checked = Ref(0)
        mismatches = Ref(0)
        probe = DiffEqBase.DiscreteCallback((u, t, integrator) -> true, function (integrator)
            run_args = integrator.p.args
            model = run_args.interlink_model
            sample = model.history[end]
            work = run_args.scheduling_policy_model.matching
            for (key, score) in zip(sample.active, sample.scores)
                checked[] += 1
                score === work.score[model.graph.edge_index_lookup[key]] || (mismatches[] += 1)
            end
            length(sample.scores) == length(sample.active) || (mismatches[] += 1)
        end; save_positions=(false, false))
        run_simulation(args; extra_callbacks=(probe,))
        @test checked[] > 0
        @test mismatches[] == 0
    end

    @testset "a weighted policy's terms" begin
        args = gve_sma_configuration(duration=20.0)
        spacecraft = args.dynamics_model.spacecraft
        state = SE.build_initial_conditions(args)
        mu = args.environment_model.planet.μ
        scored = map((false, true)) do record_terms
            model = laser_model(spacecraft; record_terms)
            scheduler = SchedulingPolicyModel(MIX)
            update_availability!(model, spacecraft, state)
            score_candidates!(model, scheduler, state, mu; spacecraft)
            (model=model, work=scheduler.matching)
        end
        plain, recorded = scored
        @test plain.work.terms === nothing
        @test size(recorded.work.terms) == (2, 2)
        # Recording the terms leaves the mix's scores unchanged, bit for bit.
        @test bits(plain.work.score) == bits(recorded.work.score)
        @test any(!iszero, recorded.work.score)
        # Each row is its term's own score.
        for (k, term) in enumerate(MIX.terms)
            alone = SchedulingPolicyModel(term.policy)
            score_candidates!(recorded.model, alone, state, mu; spacecraft)
            @test bits(recorded.work.terms[k, :]) == bits(alone.matching.score)
        end
        # A nested mix is one term, its own mixed score.
        nested = WeightedPolicy(PolicyTerm(MIX; scale=1.0), PolicyTerm(DataThroughputPolicy(); scale=1e6))
        model = laser_model(spacecraft; record_terms=true)
        scheduler = SchedulingPolicyModel(nested)
        update_availability!(model, spacecraft, state)
        score_candidates!(model, scheduler, state, mu; spacecraft)
        @test size(scheduler.matching.terms) == (2, 2)
        @test bits(scheduler.matching.terms[1, :]) == bits(recorded.work.score)
        @test all(iszero, scheduler.matching.terms[2, :])

        # In a run: the history and the table carry the terms, and the run is the
        # same with and without recording them.
        runs = map((false, true)) do record_terms
            run_args = with_policy(args, MIX; model=laser_model(spacecraft; record_terms))
            run_simulation(run_args; return_results=true, return_solution=true)
        end
        @test bits(runs[1].solution.u[end]) == bits(runs[2].solution.u[end])
        history = runs[2].configuration.interlink_model.history
        @test all(sample -> size(sample.terms) == (2, length(sample.active)), history)
        @test all(sample -> size(sample.terms, 1) == 0, runs[1].configuration.interlink_model.history)
        table = link_schedule(runs[2])
        @test "term1_start" in names(table) && "term2_start" in names(table) && !("term3_start" in names(table))
        @test !("term1_start" in names(link_schedule(runs[1])))
        @test all(==("active"), table.window)   # no contact rows without recorded availability
        sample = history[findfirst(sample -> sample.time == table.start_s[1], history)]
        index = findfirst(==(((table.a_spacecraft[1], table.a_terminal[1]), (table.b_spacecraft[1], table.b_terminal[1]))),
            sample.active)
        @test table.term1_start[1] === sample.terms[1, index]
        @test table.term2_start[1] === sample.terms[2, index]
        @test DataFrames.metadata(table)["policy"] == "WeightedPolicy" && DataFrames.metadata(table)["score_unit"] == "1"
        # Each term column says which term it is.
        @test DataFrames.colmetadata(table, "term1_start") ==
            Dict("policy" => "SemimajorAxisRatePolicy", "score_unit" => "m/s", "weight" => 0.7, "scale" => 0.3)
        @test DataFrames.colmetadata(table, "term2_start")["weight"] == -0.35

        # Only a weighted policy has terms to record.
        refused = with_policy(args, SemimajorAxisRatePolicy(); model=laser_model(spacecraft; record_terms=true))
        @test_throws ArgumentError run_simulation(refused)
    end

    @testset "the table from a hand-built history" begin
        spacecraft = [(n_terminal=2,) for _ in 1:3]
        model = InterlinkModel()
        laser = register_link!(model, spacecraft, (1, 1), (2, 1),
            LaserMomentumLink(power_w=12.0, magnification=3.0, range_m=100.0))
        data = register_link!(model, spacecraft, (2, 2), (3, 1),
            DataLink(transmit_power_w=1.0, budget=FixedRateBudget(1e6), range_m=1e6))
        @test laser < data
        push!(model.history,
            InterLinkHistory(0.0, [laser], [laser, data], [1.0]),
            InterLinkHistory(10.0, [laser, data], [laser, data], [2.0, 5.0]),
            InterLinkHistory(20.0, [data], [data], [6.0]),
            InterLinkHistory(30.0, InterLinkKey[], [data], Float64[]),
            InterLinkHistory(40.0, [laser], InterLinkKey[], [3.0]))
        table = link_schedule(model; policy=SchedulingPolicyModel(PowerDeliveryPolicy()))
        @test names(table) == ["window", "link", "kind", "directed", "a_spacecraft", "a_terminal", "b_spacecraft",
            "b_terminal", "start_s", "end_s", "duration_s", "steps", "score_start", "score_end", "score_min", "score_max"]
        @test eltype.(eachcol(table)) == [String, Int, String, Bool, Int, Int, Int, Int, Float64, Float64, Float64,
            Int, Float64, Float64, Float64, Float64]
        @test table.window == ["active", "contact", "contact", "active", "active"]
        @test table.link == [1, 1, 2, 2, 1]
        @test table.kind == ["LaserMomentumLink", "LaserMomentumLink", "DataLink{FixedRateBudget}", "DataLink{FixedRateBudget}",
            "LaserMomentumLink"]
        @test table.directed == [false, false, true, true, false]
        @test [(row.a_spacecraft, row.a_terminal, row.b_spacecraft, row.b_terminal) for row in eachrow(table)] ==
            [(1, 1, 2, 1), (1, 1, 2, 1), (2, 2, 3, 1), (2, 2, 3, 1), (1, 1, 2, 1)]
        @test table.start_s == [0.0, 0.0, 0.0, 10.0, 40.0]
        @test table.end_s == [20.0, 20.0, 40.0, 30.0, 40.0]
        @test table.duration_s == [20.0, 20.0, 40.0, 20.0, 0.0]
        @test table.steps == [2, 2, 4, 2, 1]
        active = table.window .== "active"
        @test table.score_start[active] == [1.0, 5.0, 3.0]
        @test table.score_end[active] == [2.0, 6.0, 3.0]
        @test table.score_min[active] == [1.0, 5.0, 3.0]
        @test table.score_max[active] == [2.0, 6.0, 3.0]
        @test all(isnan, Matrix(table[.!active, [:score_start, :score_end, :score_min, :score_max]]))
        @test DataFrames.metadata(table) == Dict("policy" => "PowerDeliveryPolicy", "score_unit" => "W",
            "history_start_s" => 0.0, "history_end_s" => 40.0, "record_available" => true)
        # The same spans as contact_windows.
        for (rows, windows) in ((table[active, :], contact_windows(model.history; active=true)),
                                (table[.!active, :], contact_windows(model.history)))
            @test [(row.start_s, row.end_s, row.duration_s) for row in eachrow(rows)] ==
                [(window.start_s, window.end_s, window.duration_s) for window in windows]
        end
        # Samples without scores give NaN; terms give their columns.
        bare = InterlinkModel()
        register_link!(bare, spacecraft, (1, 1), (2, 1), LaserMomentumLink(power_w=1.0, magnification=1.0, range_m=1.0))
        push!(bare.history, InterLinkHistory(0.0, [laser]), InterLinkHistory(5.0, InterLinkKey[]))
        bare_table = link_schedule(bare)
        @test nrow(bare_table) == 1 && isnan(bare_table.score_start[1]) && bare_table.end_s[1] == 5.0
        @test DataFrames.metadata(bare_table) == Dict("policy" => "", "score_unit" => "",
            "history_start_s" => 0.0, "history_end_s" => 5.0, "record_available" => false)
        push!(bare.history, InterLinkHistory(6.0, [laser], InterLinkKey[], [4.0], [7.0; 8.0;;]))
        termed = link_schedule(bare; policy=MIX)
        @test isnan(termed.term1_start[1]) && termed.term1_start[2] == 7.0
        @test isnan(termed.term2_start[1]) && termed.term2_start[2] == 8.0
        @test DataFrames.metadata(termed)["score_unit"] == "1"
        # An empty history gives the columns and no rows, and covers no span.
        empty_table = link_schedule(InterlinkModel())
        @test size(empty_table) == (0, 16)
        @test isnan(DataFrames.metadata(empty_table)["history_start_s"]) && isnan(DataFrames.metadata(empty_table)["history_end_s"])
    end

    @testset "files, manifest and returned results" begin
        mktempdir() do directory
            args = with_settings(gve_sma_configuration(duration=600.0); results=true, results_directory=directory,
                save_csv=true)
            result = run_simulation(args; return_results=true)
            table = link_schedule(result)
            @test nrow(table) > 0
            csv_path = joinpath(directory, "simulation_results_link_schedule.csv")
            feather_path = joinpath(directory, "simulation_results_link_schedule.feather")
            @test csv_path in result.files && feather_path in result.files
            from_csv = CSV.read(csv_path, DataFrame)
            from_feather = DataFrame(Arrow.Table(feather_path))
            for frame in (from_csv, from_feather)
                @test names(frame) == names(table)
                @test frame.window == table.window && frame.link == table.link && frame.steps == table.steps
                @test frame.start_s == table.start_s && frame.end_s == table.end_s
                @test isequal(frame.score_start, table.score_start) && isequal(frame.score_max, table.score_max)
            end
            manifest = TOML.parsefile(joinpath(directory, "simulation_results.manifest.toml"))
            @test manifest["files"]["link_schedule_csv"]["sha256"] == bytes2hex(sha256(read(csv_path)))
            @test manifest["files"]["link_schedule_feather"]["path"] == feather_path
            @test manifest["files"]["link_schedule_feather"]["sha256"] == bytes2hex(sha256(read(feather_path)))
            history = result.configuration.interlink_model.history
            @test manifest["interlinks"] == Dict("policy" => "SemimajorAxisRatePolicy", "score_unit" => "m/s",
                "rows" => nrow(table), "record_available" => false, "record_terms" => false,
                "history_start_s" => first(history).time, "history_end_s" => last(history).time)
            @test (first(history).time, last(history).time) == (0.0, 600.0)
            # Without a CSV, the feather file alone.
            rm(csv_path)
            run_simulation(with_settings(args; save_csv=false))
            @test !isfile(csv_path) && isfile(feather_path)
            manifest = TOML.parsefile(joinpath(directory, "simulation_results.manifest.toml"))
            @test !haskey(manifest["files"], "link_schedule_csv") && haskey(manifest["files"], "link_schedule_feather")
        end
        # In memory only: no files, and the table still comes back.
        in_memory = run_simulation(gve_sma_configuration(duration=60.0); return_results=true)
        @test isempty(in_memory.files) && nrow(link_schedule(in_memory)) > 0
        # A run without interlinks writes no schedule and returns the empty table.
        mktempdir() do directory
            plain_args = with_settings(gve_sma_configuration(duration=60.0, with_interlinks=false); results=true,
                results_directory=directory, save_csv=true)
            plain = run_simulation(plain_args; return_results=true)
            @test !isfile(joinpath(directory, "simulation_results_link_schedule.csv"))
            @test !isfile(joinpath(directory, "simulation_results_link_schedule.feather"))
            @test !haskey(TOML.parsefile(joinpath(directory, "simulation_results.manifest.toml")), "interlinks")
            @test size(link_schedule(plain)) == (0, 16)
        end
        # The schedule is not thinned with the results rows.
        thinned = withenv("SPACEAGORA_RESULTS_FINAL_ONLY" => "1") do
            run_simulation(gve_sma_configuration(duration=600.0); return_results=true)
        end
        @test nrow(thinned.table) == 1
        @test link_schedule(thinned).steps == link_schedule(run_simulation(gve_sma_configuration(duration=600.0);
            return_results=true)).steps
    end

    @testset "long and streamed results retain their link schedule" begin
        for checkpoints in (false, true), csv in (false, true)
            mktempdir() do directory
                args = with_settings(gve_sma_configuration(duration=30.0); results=true,
                    results_directory=directory, results_layout=:long, save_csv=csv,
                    checkpoint_enabled=checkpoints, checkpoint_interval_s=10.0,
                    checkpoint_directory=directory)
                result = run_simulation(args; return_results=true)
                expected = link_schedule(result)
                prefix = joinpath(directory, "simulation_results")
                schedule_path = prefix * "_link_schedule.feather"
                manifest_path = prefix * ".manifest.toml"
                manifest = TOML.parsefile(manifest_path)
                @test manifest["layout"] == "long"
                @test nrow(expected) > 0
                @test isequal(DataFrame(Arrow.Table(schedule_path)), expected)
                @test schedule_path in result.files
                @test haskey(manifest["files"], "spacecraft_feather")
                @test manifest["files"]["link_schedule_feather"]["sha256"] == bytes2hex(sha256(read(schedule_path)))
                @test manifest["interlinks"]["rows"] == nrow(expected)
                @test haskey(manifest["files"], "link_schedule_csv") == csv
                @test (prefix * "_link_schedule.csv" in result.files) == csv
                if checkpoints
                    # A finished resume must not replace or claim old schedule files.
                    before = read.( (schedule_path, manifest_path) )
                    resumed = run_simulation(with_settings(args; resume_from_checkpoint=true); return_results=true)
                    @test read.( (schedule_path, manifest_path) ) == before
                    @test !(schedule_path in resumed.files)
                    @test isequal(link_schedule(resumed), expected)
                    # Final-write recovery restores validated history and regenerates
                    # a schedule owned by the recovered manifest.
                    unfinished = TOML.parsefile(manifest_path)
                    unfinished["streaming"] = true
                    open(io -> TOML.print(io, unfinished), manifest_path, "w")
                    rm(prefix * "_spacecraft.feather")
                    recovered = run_simulation(with_settings(args; resume_from_checkpoint=true); return_results=true)
                    recovered_manifest = TOML.parsefile(manifest_path)
                    @test read(schedule_path) == before[1]
                    @test haskey(recovered_manifest, "interlinks")
                    @test haskey(recovered_manifest["files"], "link_schedule_feather")
                    @test haskey(recovered_manifest["files"], "link_schedule_csv") == csv
                    @test schedule_path in recovered.files
                    @test isequal(recovered.table, result.table)
                end
            end
        end
        # Disabled file output still returns the long table and schedule in memory.
        memory = run_simulation(with_settings(gve_sma_configuration(duration=20.0); results_layout=:long);
            return_results=true)
        @test isempty(memory.files) && nrow(link_schedule(memory)) > 0
        for checkpoints in (false, true)
            mktempdir() do directory
                args = with_settings(gve_sma_configuration(duration=30.0); results=true,
                    results_directory=directory, results_layout=:long, checkpoint_enabled=checkpoints,
                    checkpoint_interval_s=10.0, checkpoint_directory=directory)
                crash = DiffEqBase.DiscreteCallback((u, t, integrator) -> t > 15.0,
                    integrator -> error("link schedule integration interruption"); save_positions=(false, false))
                @test_throws Exception run_simulation(args; extra_callbacks=(crash,))
                prefix = joinpath(directory, "simulation_results")
                manifest = TOML.parsefile(prefix * ".manifest.toml")
                partial = DataFrame(Arrow.Table(prefix * "_link_schedule.feather"))
                @test nrow(partial) > 0 && maximum(partial.end_s) < 30.0
                @test manifest["files"]["link_schedule_feather"]["sha256"] ==
                    bytes2hex(sha256(read(prefix * "_link_schedule.feather")))
            end
        end
        # CSV-only output has no bundle; no stale feather is claimed.
        mktempdir() do directory
            args = with_settings(gve_sma_configuration(duration=20.0); results=true, results_directory=directory,
                results_layout=:long, save_csv=true)
            result = withenv("SPACEAGORA_SAVE_BUNDLE" => "0") do
                run_simulation(args; return_results=true)
            end
            @test joinpath(directory, "simulation_results_link_schedule.csv") in result.files
            @test !isfile(joinpath(directory, "simulation_results_link_schedule.feather"))
            @test !isfile(joinpath(directory, "simulation_results.manifest.toml"))
        end
    end

    @testset "a resumed run's schedule holds the whole run" begin
        checkpointed(directory, duration; kwargs...) = with_settings(gve_sma_configuration(duration=duration);
            results=true, results_directory=directory, save_csv=true, checkpoint_enabled=true,
            checkpoint_interval_s=10.0, checkpoint_directory=directory, kwargs...)
        extend(args, duration) = SM.SimConfig._with_configuration(with_settings(args; resume_from_checkpoint=true);
            mission_configuration=gve_sma_configuration(duration=duration).mission_configuration)
        feather(directory) = read(joinpath(directory, "simulation_results_link_schedule.feather"))
        manifest(directory) = TOML.parsefile(joinpath(directory, "simulation_results.manifest.toml"))
        # A run checkpointed every 10 s and never interrupted, to 30 s.
        straight_dir = mktempdir()
        straight = link_schedule(run_simulation(checkpointed(straight_dir, 30.0); return_results=true))
        @test isfile(joinpath(straight_dir, "simulation_checkpoint.interlinks.bin"))
        # The same run stopped at 20 s and resumed to 30 s: the same schedule,
        # files and span, and the manifest names the resume.
        mktempdir() do directory
            run_simulation(checkpointed(directory, 20.0))
            resumed = run_simulation(extend(checkpointed(directory, 20.0), 30.0); return_results=true)
            table = link_schedule(resumed)
            @test isequal(table, straight)
            @test DataFrames.metadata(table) == DataFrames.metadata(straight)
            @test DataFrames.metadata(table)["history_start_s"] == 0.0
            @test feather(directory) == feather(straight_dir)
            @test manifest(directory)["resumed_from_s"] == 20.0
            @test manifest(directory)["interlinks"] == manifest(straight_dir)["interlinks"]
            @test joinpath(directory, "simulation_checkpoint.interlinks.bin") in resumed.files
        end
        # A run that fails after its 20 s checkpoint resumes to the same schedule:
        # the samples after the checkpoint are dropped and computed again.
        mktempdir() do directory
            stop = DiffEqBase.DiscreteCallback((u, t, integrator) -> t > 25.0, integrator -> error("stop here");
                save_positions=(false, false))
            @test_throws Exception run_simulation(checkpointed(directory, 30.0); extra_callbacks=(stop,))
            resumed = run_simulation(extend(checkpointed(directory, 30.0), 30.0); return_results=true)
            @test isequal(link_schedule(resumed), straight)
        end
        # Without its record of the history (a checkpoint written before
        # SpaceAGORA kept it, a missing file, another run's file), a resume
        # warns, and its schedule covers the resumed interval and says so.
        for spoil! in (directory -> rm(joinpath(directory, "simulation_checkpoint.interlinks.bin")),
                       directory -> cp(joinpath(straight_dir, "simulation_checkpoint.interlinks.bin"),
                           joinpath(directory, "simulation_checkpoint.interlinks.bin"); force=true))
            mktempdir() do directory
                run_simulation(checkpointed(directory, 20.0))
                spoil!(directory)
                resumed = @test_logs (:warn, r"link schedule cannot include the windows before its checkpoint") match_mode=:any begin
                    run_simulation(extend(checkpointed(directory, 20.0), 30.0); return_results=true)
                end
                table = link_schedule(resumed)
                @test minimum(table.start_s) == 20.0 && maximum(table.end_s) == 30.0
                @test DataFrames.metadata(table)["history_start_s"] == 20.0
                @test manifest(directory)["interlinks"]["history_start_s"] == 20.0 == manifest(directory)["resumed_from_s"]
                @test resumed.table.time[1] == 0.0                 # the results still hold the whole run
            end
        end
        # A sample the scheduler replaces after it was kept (at a segment's start,
        # at the same time) is kept again, not left as it was.
        mktempdir() do directory
            path = joinpath(directory, "simulation_checkpoint.interlinks.bin")
            key = ((1, 1), (2, 1))
            model = InterlinkModel()
            push!(model.history, InterLinkHistory(0.0, [key], InterLinkKey[], [1.0]),
                InterLinkHistory(10.0, [key], InterLinkKey[], [2.0]))
            log = SE._InterlinkHistoryLog(path, "token", true, 0, nothing, 0, true, 0.0, false)
            first_entry = SE._append_interlink_history!(log, model)
            @test first_entry.sha256 == bytes2hex(sha256(read(path)))
            model.history[end] = InterLinkHistory(10.0, InterLinkKey[], InterLinkKey[], Float64[])
            push!(model.history, InterLinkHistory(20.0, [key], InterLinkKey[], [3.0]))
            entry = SE._append_interlink_history!(log, model)
            @test entry.sha256 == bytes2hex(sha256(read(path)))
            stream = SM.IOSerialization._read_checkpoint_rows(path, entry.bytes)
            @test entry.samples == 3 && [sample.time for sample in stream.data] == [0.0, 10.0, 20.0]
            @test isempty(stream.data[2].active) && stream.data[3].scores == [3.0]
            # The stream is taken back only when it reaches the checkpoint's time
            # with the number of samples the checkpoint records.
            kept, reason = SE._checkpoint_interlink_history(entry, path, 20.0, model)
            @test kept !== nothing && reason === nothing
            for (doctored, t) in ((merge(entry, (samples=entry.samples + 1,)), 20.0), (entry, 25.0))
                refused, why = SE._checkpoint_interlink_history(doctored, path, t, model)
                @test refused === nothing && occursin("does not reach the checkpoint's time", why)
            end
            @test occursin("another run", SE._checkpoint_interlink_history(merge(entry, (token="other",)), path, 20.0, model)[2])
            @test SE._checkpoint_interlink_history(nothing, path, 20.0, model)[1] === nothing
            @test SE._checkpoint_interlink_history(merge(entry, (version=1,)), path, 20.0, model)[1] === nothing
            for malformed in (nothing, missing, "", "a"^63, "z"^64, 123)
                @test occursin("malformed", SE._checkpoint_interlink_history(
                    merge(entry, (sha256=malformed,)), path, 20.0, model)[2])
            end
            @test occursin("byte length", SE._checkpoint_interlink_history(
                merge(entry, (bytes=-1,)), path, 20.0, model)[2])
            original = read(path)
            # Flip a bit in a score, leaving a readable stream of the same size.
            score_bytes = collect(reinterpret(UInt8, [3.0]))
            locations = [i for i in 1:(length(original)-7) if original[i:i+7] == score_bytes]
            @test length(locations) == 1
            corrupt = copy(original)
            corrupt[only(locations)] ⊻= 0x01
            write(path, corrupt)
            readable = SM.IOSerialization._read_checkpoint_rows(path, entry.bytes)
            @test readable !== nothing && readable.data[3].scores != [3.0]
            @test occursin("checksum", SE._checkpoint_interlink_history(entry, path, 20.0, model)[2])
            settings = (checkpoint_directory=directory, results_directory=directory)
            fresh_args() = (simulation_settings=settings, interlink_model=InterlinkModel())
            no_solve = fresh_args()
            @test_logs (:warn, r"checksum.*This resume solves nothing") SE._start_interlink_history!(
                no_solve, true, 20.0, entry; solves=false)
            @test isempty(no_solve.interlink_model.history) && read(path) == corrupt
            resumed_args = fresh_args()
            replacement = @test_logs (:warn, r"checksum") SE._start_interlink_history!(resumed_args, true, 20.0, entry)
            @test replacement.start_s == 20.0 && read(path) == corrupt
            push!(resumed_args.interlink_model.history, InterLinkHistory(20.0, [key], InterLinkKey[], [3.0]))
            replaced_entry = SE._append_interlink_history!(replacement, resumed_args.interlink_model)
            archives = filter(p -> occursin("interlinks.before_resume", p), readdir(directory; join=true))
            @test length(archives) == 1 && read(only(archives)) == corrupt
            @test replaced_entry.sha256 == bytes2hex(sha256(read(path)))

            # Legacy metadata retains readable history with an explicit warning;
            # a subsequent checkpoint protects the retained prefix and new bytes.
            write(path, original)
            legacy = Base.structdiff(entry, NamedTuple{(:sha256,)})
            legacy_args = fresh_args()
            legacy_log = @test_logs (:warn, r"has no checksum") SE._start_interlink_history!(legacy_args, true, 20.0, legacy)
            @test length(legacy_args.interlink_model.history) == entry.samples
            push!(legacy_args.interlink_model.history, InterLinkHistory(30.0, [key], InterLinkKey[], [4.0]))
            upgraded = SE._append_interlink_history!(legacy_log, legacy_args.interlink_model)
            @test upgraded.sha256 == bytes2hex(sha256(read(path)))
            @test SE._checkpoint_interlink_history(upgraded, path, 30.0, legacy_args.interlink_model)[1] !== nothing

            # A no-solve resume keeps a crash tail; a continuing one removes only
            # that tail and initializes the digest for the next append correctly.
            write(path, original)
            open(io -> write(io, "after checkpoint"), path, "a")
            with_tail = read(path)
            valid_no_solve = fresh_args()
            SE._start_interlink_history!(valid_no_solve, true, 20.0, entry; solves=false)
            @test read(path) == with_tail && length(valid_no_solve.interlink_model.history) == entry.samples
            valid_resume = fresh_args()
            resumed_log = SE._start_interlink_history!(valid_resume, true, 20.0, entry)
            @test read(path) == original
            push!(valid_resume.interlink_model.history, InterLinkHistory(30.0, [key], InterLinkKey[], [4.0]))
            continued = SE._append_interlink_history!(resumed_log, valid_resume.interlink_model)
            @test continued.sha256 == bytes2hex(sha256(read(path)))
            @test SE._checkpoint_interlink_history(continued, path, 30.0, valid_resume.interlink_model)[1] !== nothing
            write(path, original)
            open(io -> write(io, "after checkpoint"), path, "a")
            @test all(all(isequal(getfield(a, field), getfield(b, field)) for field in fieldnames(typeof(a))) for (a, b) in zip(
                SE._checkpoint_interlink_history(entry, path, 20.0, model)[1].data, stream.data))
            open(io -> write(io, original[1:end-1]), path, "w")
            @test SE._checkpoint_interlink_history(entry, path, 20.0, model)[1] === nothing

        end
        # Weighted terms and available contacts survive resume together.
        mktempdir() do directory
            original = with_policy(checkpointed(directory, 20.0), MIX)
            original.interlink_model.record_available = true
            original.interlink_model.record_terms = true
            first = run_simulation(original; return_results=true)
            kept = deepcopy(first.configuration.interlink_model.history)
            resumed = run_simulation(extend(original, 30.0); return_results=true)
            history = resumed.configuration.interlink_model.history
            @test history[1].time == 0.0 && history[end].time == 30.0
            @test all(all(isequal(getfield(a, field), getfield(b, field)) for field in fieldnames(typeof(a)))
                for (a, b) in zip(kept[1:end-1], history[1:length(kept)-1]))
            @test any(!isempty(sample.available) for sample in history)
            @test all(size(sample.terms, 1) == length(MIX.terms) for sample in history)
        end
        # Invalid old bytes remain until replacement, including a failed first segment.
        mktempdir() do directory
            args = checkpointed(directory, 20.0)
            run_simulation(args)
            path = SE._checkpoint_paths(args).interlinks
            write(path, "unusable earlier history")
            original = read(path)
            stop = DiffEqBase.DiscreteCallback((u, t, integrator) -> t > 20.0,
                integrator -> error("first resumed segment failed"); save_positions=(false, false))
            @test_throws Exception run_simulation(extend(args, 30.0); extra_callbacks=(stop,))
            @test read(path) == original
            run_simulation(extend(args, 30.0))
            kept = filter(p -> occursin("interlinks.before_resume", p), readdir(directory; join=true))
            @test length(kept) == 1 && read(only(kept)) == original
        end
        # Recording identity cannot be changed across resume; refusal must not
        # truncate any checkpoint stream, even trailing crash bytes.
        for flag in (:record_available, :record_terms), old_value in (false, true)
            mktempdir() do directory
                args = with_policy(checkpointed(directory, 20.0), MIX)
                setproperty!(args.interlink_model, flag, old_value)
                run_simulation(args)
                paths = SE._checkpoint_paths(args)
                open(io -> write(io, "crash tail"), paths.interlinks, "a")
                paths_to_check = filter(isfile, collect(values(paths)))
                before = Dict(path => read(path) for path in paths_to_check)
                changed = extend(with_policy(checkpointed(directory, 20.0), MIX), 30.0)
                setproperty!(changed.interlink_model, flag, !old_value)
                @test_throws r"different record_available or record_terms" run_simulation(changed)
                @test all(read(path) == bytes for (path, bytes) in before)
                rm(paths.interlinks)
                delete!(before, paths.interlinks)
                @test_throws r"different record_available or record_terms" run_simulation(changed)
                @test all(read(path) == bytes for (path, bytes) in before)
                @test !isfile(paths.interlinks)
            end
        end
        # A complete no-op returns the retained schedule without changing files;
        # missing history gives an empty schedule and still changes no files.
        mktempdir() do directory
            args = checkpointed(directory, 20.0)
            result = run_simulation(args; return_results=true)
            paths = filter(isfile, readdir(directory; join=true))
            before = Dict(path => read(path) for path in paths)
            resumed = run_simulation(with_settings(args; resume_from_checkpoint=true); return_results=true)
            @test isequal(link_schedule(resumed), link_schedule(result))
            @test all(read(path) == bytes for (path, bytes) in before)
            @test isempty(resumed.files)
            rm(SE._checkpoint_paths(args).interlinks)
            empty_resume = @test_logs (:warn, r"This resume solves nothing") match_mode=:any begin
                run_simulation(with_settings(args; resume_from_checkpoint=true); return_results=true)
            end
            @test isempty(link_schedule(empty_resume)) && isempty(empty_resume.files)
        end
        # Recovery without valid history preserves but does not claim prior schedules.
        mktempdir() do directory
            args = with_settings(checkpointed(directory, 20.0); results_layout=:long)
            result = run_simulation(args; return_results=true)
            prefix = joinpath(directory, "simulation_results")
            schedule_path = prefix * "_link_schedule.feather"
            before = read(schedule_path)
            unfinished = manifest(directory)
            unfinished["streaming"] = true
            open(io -> TOML.print(io, unfinished), prefix * ".manifest.toml", "w")
            rm(prefix * "_spacecraft.feather")
            rm(SE._checkpoint_paths(args).interlinks)
            recovered = @test_logs (:warn, r"This resume solves nothing") match_mode=:any begin
                run_simulation(with_settings(args; resume_from_checkpoint=true); return_results=true)
            end
            @test read(schedule_path) == before
            @test !haskey(manifest(directory), "interlinks")
            @test !haskey(manifest(directory)["files"], "link_schedule_feather")
            @test !(schedule_path in recovered.files)
            @test isempty(link_schedule(recovered))
            @test isequal(recovered.table, result.table)
        end
        # A run without interlinks keeps no history beside its checkpoint.
        mktempdir() do directory
            run_simulation(with_settings(gve_sma_configuration(duration=20.0, with_interlinks=false); results=false,
                checkpoint_enabled=true, checkpoint_interval_s=10.0, checkpoint_directory=directory))
            @test !isfile(joinpath(directory, "simulation_checkpoint.interlinks.bin"))
            @test !hasproperty(SE._load_checkpoint(with_settings(gve_sma_configuration(duration=20.0,
                with_interlinks=false); checkpoint_enabled=true, checkpoint_directory=directory)).runtime_state,
                :interlink_history)
        end
    end

    @testset "edge cases of the written schedule" begin
        # A link that never closes writes a header-only CSV and an empty feather file.
        mktempdir() do directory
            args = gve_sma_configuration(duration=60.0)
            never = InterlinkModel()
            register_link!(never, args.dynamics_model.spacecraft, (1, 1), (2, 1),
                LaserMomentumLink(power_w=1.0, magnification=1.0, range_m=1.0))
            run_simulation(with_settings(with_policy(args, SemimajorAxisRatePolicy(); model=never); results=true,
                results_directory=directory, save_csv=true))
            csv = CSV.read(joinpath(directory, "simulation_results_link_schedule.csv"), DataFrame)
            feather = DataFrame(Arrow.Table(joinpath(directory, "simulation_results_link_schedule.feather")))
            @test nrow(csv) == 0 && length(names(csv)) == 16
            @test nrow(feather) == 0 && names(feather) == names(link_schedule(InterlinkModel()))
            @test TOML.parsefile(joinpath(directory, "simulation_results.manifest.toml"))["interlinks"]["rows"] == 0
        end
        # A solve that fails part-way still writes the schedule up to its last sample.
        mktempdir() do directory
            stop = DiffEqBase.DiscreteCallback((u, t, integrator) -> t > 10.0, integrator -> error("stop here");
                save_positions=(false, false))
            args = with_settings(gve_sma_configuration(duration=60.0); results=true, results_directory=directory,
                save_csv=true)
            @test_throws Exception run_simulation(args; extra_callbacks=(stop,))
            partial = DataFrame(Arrow.Table(joinpath(directory, "simulation_results_link_schedule.feather")))
            @test nrow(partial) > 0 && minimum(partial.start_s) == 0.0 && 0 < maximum(partial.end_s) < 11.0
        end
        # Checkpoint segments within one run do not split its windows.
        unbroken = link_schedule(run_simulation(gve_sma_configuration(duration=60.0); return_results=true))
        mktempdir() do directory
            segmented = link_schedule(run_simulation(with_settings(gve_sma_configuration(duration=60.0); results=false,
                checkpoint_enabled=true, checkpoint_interval_s=10.0, checkpoint_directory=directory); return_results=true))
            @test segmented[:, [:window, :link, :start_s, :end_s]] == unbroken[:, [:window, :link, :start_s, :end_s]]
        end
        # Of a directed pair, the schedule records the selected direction and its score.
        args = gve_sma_configuration(duration=60.0)
        spacecraft = args.dynamics_model.spacecraft
        pair = InterlinkModel()
        register_link!(pair, spacecraft, (1, 1), (2, 1), DataLink(transmit_power_w=1.0, budget=FixedRateBudget(1e6), range_m=1e6))
        register_link!(pair, spacecraft, (2, 1), (1, 1), DataLink(transmit_power_w=1.0, budget=FixedRateBudget(2e6), range_m=1e6))
        directed = link_schedule(run_simulation(with_policy(args, DataThroughputPolicy(); model=pair); return_results=true))
        active = directed[directed.window .== "active", :]
        @test nrow(active) == 1 && active.a_spacecraft[1] == 2 && active.b_spacecraft[1] == 1
        @test active.score_start[1] == 2e6 && active.score_min[1] == 2e6 && active.directed[1]
    end
end

end # module InterlinkScheduleTests
