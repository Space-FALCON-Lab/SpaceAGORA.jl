module BenchmarkWrapperColumnsTests
using Test
using CSV
using DataFrames
using Dates
using Statistics

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RUNTIME = joinpath(REPO, "benchmarks", "studies", "performance_runtime_analysis")
const PIPELINE = joinpath(REPO, "benchmarks", "scripts", "performance_paper_pipeline.jl")
const STATIC = joinpath(REPO, "benchmarks", "studies", "performance_static_vs_parallel.jl")
const REPORT_FIXTURE = joinpath(@__DIR__, "benchmark_report_tests.jl")

# Load selected real definitions and reuse the existing synthetic row schema.
# The wrappers' entrypoint includes would otherwise load the native benchmark.
function load_definition(file::String, wanted::Symbol)
    nodes = Meta.parseall(read(file, String)).args
    for expr in nodes
        if expr isa Expr && expr.head == :module
            nodes = expr.args[end].args
            break
        end
    end
    for expr in nodes
        expr isa Expr || continue
        node = expr.head == :macrocall ? expr.args[end] : expr
        node isa Expr || continue
        name = if node.head == :struct
            node.args[2]
        elseif node.head in (:function, :(=))
            signature = node.args[1]
            signature isa Expr || continue
            signature = signature.head == :(::) ? signature.args[1] : signature
            signature.head == :call ? signature.args[1] : nothing
        else
            nothing
        end
        if name == wanted
            Core.eval(@__MODULE__, expr)
            return
        end
    end
    error("Required wrapper definition $wanted not found in $file")
end
const HELPERS = joinpath(RUNTIME, "case_catalog", "profile_types_and_helpers.jl")
load_definition(HELPERS, :ProfileSpec)
load_definition(HELPERS, :_safe_unique_join)
struct BenchmarkCase end
const PERF_BASELINE_SCENARIO = "baseline_fixture"
load_definition(REPORT_FIXTURE, :raw_row)
include(joinpath(RUNTIME, "reporting", "summaries.jl"))
for name in (:PipelineConfig, :ModeRunArtifacts, :_fmt_md, :_write_markdown_table,
             :build_mode_overview, :write_pipeline_report)
    load_definition(PIPELINE, name)
end
for name in (:StaticVsParallelConfig, :ArmSpec, :PolicyMatrixSpec, :ArmPassResult,
             :_latest_artifact_path, :_latest_artifact_path_optional, :_same_run_artifact_path,
             :_stage_elapsed_s, :_arm_result_artifacts, :_tag_arm_column, :_write_aggregate_arm_report,
             :_aggregate_arm_artifacts, :_write_static_vs_parallel_report)
    load_definition(STATIC, name)
end
const SPEC = ProfileSpec(name="quick", repeats=1, warmup=0, max_attempts=2,
    mission_short_s=120.0, mission_long_s=600.0, montecarlo_samples=0,
    montecarlo_mission_s=120.0)
const FALLBACK = (machine_label="fallback-machine", hardware_class="fallback-class",
    host_name="fallback-host", cpu_name="fallback-cpu", cpu_threads=3,
    julia_threads=1, os="fallback-os", arch="fallback-arch")
_machine_label() = FALLBACK.machine_label
_hardware_class_name() = FALLBACK.hardware_class
_runtime_hardware_snapshot() = FALLBACK
const ARM = ArmSpec(mode=:auto, label="auto_fixture", backend="none", adaptive=false)
const MATRIX = PolicyMatrixSpec(key=:outer_pinned, label="fixture", density_mode="off",
    control_mode="off", multibody_mode="off", effector_mode="off")
function static_config(dir)
    StaticVsParallelConfig(profile=SPEC, outdir=dir, clean=false, include_process=false,
        process_workers=nothing, passes=2, randomize_arm_order=false, random_seed=17,
        policy_matrices=[:outer_pinned], include_control_stress_per_orbit=false,
        control_stress_repeats_full=1, control_stress_warmup_full=0)
end
function recorded_rows(; failed_only=false)
    raw = DataFrame([raw_row(solve_success=false, total_time_s=100.0, solve_retcode="Unstable"),
        raw_row(solve_success=!failed_only, total_time_s=2.0, attempt=2,
                solve_retcode=failed_only ? "FixtureFailure" : "Success")])
    raw[!, :machine_label] = fill("recorded-machine", 2)
    raw[!, :hardware_class] = fill("recorded-class", 2)
    raw[!, :cpu_threads] = fill(64, 2)
    raw[!, :julia_threads] = fill(8, 2)
    return raw
end
function entry_rows(pass_idx=1)
    reference = 100.0 * pass_idx
    error = pass_idx == 1 ? 3.0 : 7.0
    second_error = pass_idx == 1 ? 1.0 : 9.0
    rows = DataFrame([raw_row(category="entry", scenario="entry_fixture",
        description="Synthetic saved entry-duration fixture", repeat=1,
        entry_atmospheric_interface_count=count, entry_run_role=role,
        solve_success=success, attempt=attempt, is_terminal_attempt=terminal,
        solve_retcode=success ? "Success" : "Unstable", total_time_s=cost,
        entry_passage_duration_s=duration, entry_wall_time_per_passage_s=wall,
        entry_event_time_abs_error_s=event_error,
        entry_reference_terminal_time_s=reference_time, terminal_time_s=event_time)
        for (count, role, success, attempt, terminal, cost, duration, wall, event_error, reference_time, event_time) in (
            (1,"reference",true,1,true,4.0,30.0,4.0,0.0,reference,reference),
            (1,"measured",false,1,false,1000.0 * pass_idx,999.0,999.0,999.0,reference,999.0),
            (1,"measured",true,2,true,2.0 * pass_idx,20.0 * pass_idx,2.0 * pass_idx,error,reference,reference + error),
            (3,"reference",true,1,true,6.0,60.0,2.0,0.0,reference + 10.0,reference + 10.0),
            (3,"measured",true,1,true,9.0,90.0,3.0,second_error,reference + 10.0,reference + 10.0 + second_error),
            (5,"measured",false,2,true,50.0,missing,missing,missing,missing,missing))])
    return rows
end
function artifact(dir; raw=recorded_rows(), gate=DataFrame(pass_all=[true, false]),
                  entry_raw=entry_rows(), bench=2.0, split=1.0, orbit=4.0, entry=3.0,
                  total=bench+split+orbit+entry)
    mkpath(dir)
    orbit_raw = copy(raw)
    orbit_raw[!, :orbit_count] = fill(1, nrow(raw))
    orbit_raw[!, :mission_time_multiplier] = fill(1, nrow(raw))
    orbit_raw[!, :orbital_period_s] = fill(120.0, nrow(raw))
    orbit_path = joinpath(dir, "orbit_raw.csv")
    CSV.write(orbit_path, orbit_raw)
    summary = nrow(raw) == 0 ? summarize_results(recorded_rows())[1:0, :] : summarize_results(raw)
    entry_raw_path = joinpath(dir, "entry_raw.csv")
    entry_summary_path = joinpath(dir, "entry_summary.csv")
    CSV.write(entry_raw_path, entry_raw)
    CSV.write(entry_summary_path, summarize_entry_duration_results(copy(entry_raw)))
    ModeRunArtifacts(mode=:auto, backend="none", elapsed_s=total,
        bench_elapsed_s=bench, split_gate_elapsed_s=split, orbit_elapsed_s=orbit,
        entry_duration_elapsed_s=entry, raw_path=joinpath(dir,"raw.csv"),
        summary_path=joinpath(dir,"summary.csv"), orbit_raw_path=orbit_path,
        orbit_summary_path=joinpath(dir,"orbit_summary.csv"),
        entry_duration_raw_path=entry_raw_path,
        entry_duration_summary_path=entry_summary_path,
        report_path=joinpath(dir,"source.md"), split_gate_df=gate,
        raw_df=raw, summary_df=summary, orbit_summary_df=DataFrame())
end
function reports(dir, artifacts; overview=build_mode_overview(artifacts))
    comparison = DataFrame(scenario=String[])
    before = deepcopy((overview, comparison, [a.raw_df for a in artifacts], [a.split_gate_df for a in artifacts]))
    pipeline_path = joinpath(dir, "pipeline.md")
    static_path = joinpath(dir, "static.md")
    write_pipeline_report(pipeline_path, PipelineConfig(profile=SPEC,modes=[:auto],outdir=dir), overview, comparison, artifacts)
    _write_static_vs_parallel_report(static_path, static_config(dir), MATRIX, [ARM],
        overview, comparison, artifacts, DataFrame(pass=[1],arm=[ARM.label]))
    @test isequal((overview, comparison, [a.raw_df for a in artifacts], [a.split_gate_df for a in artifacts]), before)
    return read(pipeline_path, String), read(static_path, String)
end

@testset "stage timings require recorded complete finite costs" begin
    mktempdir() do dir
        path = joinpath(dir, "stages.csv")
        @test_throws ArgumentError _stage_elapsed_s(path)
        valid = DataFrame(stage=["run_per_orbit", "run_benchmarks", "run_entry_duration_sweep", "total"],
                          elapsed_s=[4.0,2.0,3.0,12.0])
        for (table, expected) in (
            (valid,(2.0,0.0,4.0,3.0,12.0)),
            (vcat(valid, DataFrame(stage=["run_split_rollout_gate", "run_multirate_rollout_gate"],elapsed_s=[1.0,2.0])),(2.0,1.0,4.0,3.0,12.0)),
            (DataFrame(stage=copy(valid.stage),elapsed_s=zeros(4)),(0.0,0.0,0.0,0.0,0.0)))
            CSV.write(path, table)
            original = read(path)
            @test _stage_elapsed_s(path) == expected
            @test read(path) == original
        end
        malformed = DataFrame[DataFrame(), DataFrame(stage=["total"]), DataFrame(elapsed_s=[1.0]),
            vcat(valid, DataFrame(stage=[missing], elapsed_s=[1.0])),
            DataFrame(stage=String[],elapsed_s=Float64[])]
        for name in valid.stage
            push!(malformed, filter(row -> row.stage != name, valid))
            push!(malformed, vcat(valid, filter(row -> row.stage == name, valid)))
        end
        push!(malformed, vcat(valid, DataFrame(stage=fill("run_split_rollout_gate",2),elapsed_s=[1.0,2.0])))
        for name in vcat(valid.stage, ["run_split_rollout_gate"]), bad in (-1.0, Inf, NaN, missing, "invalid")
            table = vcat(valid, DataFrame(stage=["run_split_rollout_gate"], elapsed_s=[1.0]))
            table[!, :elapsed_s] = Any[table.elapsed_s...]
            table[findfirst(==(name),table.stage), :elapsed_s] = bad
            push!(malformed, table)
        end
        for table in malformed
            CSV.write(path, table)
            original = read(path)
            @test_throws ArgumentError _stage_elapsed_s(path)
            @test read(path) == original
        end
    end
end

# Files have the same names and stamp relationship as the real runtime CLI.
# Gate output deliberately has its own timestamp, as its producer does.
function save_runtime_run(dir; stamp="20261007_120000", entry_raw=entry_rows(), hardware=true)
    mkpath(dir)
    path(prefix, suffix=".csv") = joinpath(dir,"$(prefix)_$(SPEC.name)_$(stamp)$(suffix)")
    raw = recorded_rows()
    orbit = copy(raw)
    orbit[!, :orbit_count] = fill(1,nrow(raw))
    orbit[!, :mission_time_multiplier] = fill(1,nrow(raw))
    orbit[!, :orbital_period_s] = fill(120.0,nrow(raw))
    for (prefix, table) in (("runtime_raw",raw), ("runtime_summary",summarize_results(raw)),
        ("runtime_per_orbit_raw",orbit), ("runtime_per_orbit_summary",summarize_per_orbit_results(orbit)),
        ("runtime_entry_duration_raw",entry_raw),
        ("runtime_entry_duration_summary",summarize_entry_duration_results(copy(entry_raw))),
        ("runtime_stage_timing",DataFrame(stage=["run_benchmarks","run_per_orbit","run_entry_duration_sweep","run_multirate_rollout_gate","total"],elapsed_s=[2.0,4.0,3.0,4.0,13.0])))
        CSV.write(path(prefix),table)
    end
    hardware && CSV.write(path("runtime_hardware_info"),DataFrame(machine_label=["same-run-machine"]))
    write(path("runtime_report",".md"),"Recorded runtime report")
    CSV.write(joinpath(dir,"split_rollout_gate_$(SPEC.name)_20261007_115959.csv"),DataFrame(pass_all=[true,false]))
    write(joinpath(dir,"split_rollout_gate_$(SPEC.name)_20261007_115959.md"),"Recorded gate report")
    return path
end

@testset "arm constructor uses complete same-run entry artifacts" begin
    mktempdir() do dir
        path = save_runtime_run(dir)
        # Newer unrelated companions must not replace any selected run's files.
        for prefix in ("runtime_summary","runtime_per_orbit_raw","runtime_per_orbit_summary",
                       "runtime_entry_duration_raw","runtime_entry_duration_summary",
                       "runtime_stage_timing","runtime_hardware_info")
            CSV.write(joinpath(dir,"$(prefix)_$(SPEC.name)_other_stamp.csv"),DataFrame(contamination=[true]))
        end
        write(joinpath(dir,"runtime_report_$(SPEC.name)_other_stamp.md"),"Wrong run")
        a = _arm_result_artifacts(ARM,static_config(dir),dir)
        @test a isa ModeRunArtifacts
        @test a.raw_path == path("runtime_raw")
        @test a.summary_path == path("runtime_summary")
        @test a.orbit_raw_path == path("runtime_per_orbit_raw")
        @test a.orbit_summary_path == path("runtime_per_orbit_summary")
        @test a.entry_duration_raw_path == path("runtime_entry_duration_raw")
        @test a.entry_duration_summary_path == path("runtime_entry_duration_summary")
        @test a.stage_timing_path == path("runtime_stage_timing")
        @test a.hardware_info_path == path("runtime_hardware_info")
        @test a.report_path == path("runtime_report",".md")
        @test (a.bench_elapsed_s,a.split_gate_elapsed_s,a.orbit_elapsed_s,a.entry_duration_elapsed_s,a.elapsed_s) == (2.0,0.0,4.0,3.0,13.0)
        @test isequal(CSV.read(a.entry_duration_raw_path,DataFrame),entry_rows())
        @test a.split_gate_df.pass_all == [true,false]
        @test endswith(a.split_gate_csv_path,"20261007_115959.csv")
        @test endswith(a.split_gate_report_path,"20261007_115959.md")
        rm(path("runtime_hardware_info"))
        @test isempty(_arm_result_artifacts(ARM,static_config(dir),dir).hardware_info_path)
        for prefix in ("runtime_entry_duration_raw","runtime_entry_duration_summary",
                       "runtime_stage_timing","runtime_summary","runtime_per_orbit_raw","runtime_per_orbit_summary")
            file = path(prefix)
            bytes = read(file)
            rm(file)
            @test_throws ArgumentError _arm_result_artifacts(ARM,static_config(dir),dir)
            write(file,bytes)
        end
        file = path("runtime_report",".md")
        bytes = read(file)
        rm(file)
        @test_throws ArgumentError _arm_result_artifacts(ARM,static_config(dir),dir)
        write(file,bytes)
        rm(path("runtime_raw"))
        @test_throws ErrorException _arm_result_artifacts(ARM,static_config(dir),dir)
    end
end

@testset "saved wrapper artifacts preserve provenance failures and gate counts" begin
    mktempdir() do dir
        a = artifact(joinpath(dir,"mixed"))
        overview = build_mode_overview([a])
        @test overview.machine_label == ["recorded-machine"]
        @test overview.hardware_class == ["recorded-class"]
        @test overview.failed_rows == [1] && overview.unstable_rows == [1]
        @test overview.rows_raw == [2]
        @test overview.baseline_mean_s == [2.0]
        @test overview.split_gate_rows == [2]
        @test overview.split_gate_pass_rows == [1]
        @test overview.split_gate_pass_rate_pct == [50.0]
        @test overview.orbit_share_pct == [40.0] && overview.entry_duration_share_pct == [30.0]
        pipeline, static = reports(dir, [a]; overview=overview)
        @test occursin("Hardware classes observed: `recorded-class`", pipeline)
        @test occursin("Machine labels observed: `recorded-machine`", pipeline)
        @test occursin("Split rollout gate rows: `2`; pass rows: `1` (`50.0%`)", pipeline)
        @test occursin("Solver-success samples across modes: `1/2`", pipeline)
        @test occursin("| auto | recorded-machine | recorded-class | 64 | 8 |", static)
        @test occursin("Split rollout gate rows (aggregated): `1/2` pass (`50.0%`)", static)
        @test !occursin("fallback-machine", pipeline)
        failed = artifact(joinpath(dir,"failed"); raw=recorded_rows(failed_only=true),gate=DataFrame(pass_all=[false,false]))
        failed_overview = build_mode_overview([failed])
        @test failed_overview.failed_rows == [2]
        @test ismissing(only(failed_overview.baseline_mean_s))
        @test failed_overview.split_gate_pass_rows == [0]
        failed_pipeline, failed_static = reports(dir, [failed])
        @test occursin("Solver-success samples across modes: `0/2`", failed_pipeline)
        @test occursin("Split rollout gate rows (aggregated): `0/2` pass (`0.0%`)", failed_static)
        both_pipeline, both_static = reports(dir, [a, failed])
        @test occursin("Split rollout gate rows: `4`; pass rows: `1` (`25.0%`)", both_pipeline)
        @test occursin("Solver-success samples across modes: `1/4`", both_pipeline)
        @test occursin("Split rollout gate rows (aggregated): `1/4` pass (`25.0%`)", both_static)
        for gate in (nothing, DataFrame(pass_all=Bool[]), DataFrame(other=[true]))
            legacy = artifact(joinpath(dir,"legacy"); raw=select(recorded_rows(),Not([:machine_label,:hardware_class,:cpu_threads,:julia_threads])), gate=gate)
            legacy_overview = build_mode_overview([legacy])
            @test legacy_overview.machine_label == ["fallback-machine"]
            @test legacy_overview.hardware_class == ["fallback-class"]
            @test legacy_overview.split_gate_rows == [0] && legacy_overview.split_gate_pass_rows == [0]
            @test ismissing(only(legacy_overview.split_gate_pass_rate_pct))
            p, s = reports(dir,[legacy])
            @test occursin("Split rollout gate rows: none",p)
            @test occursin("| auto | n/a | n/a | n/a | n/a |",s)
        end
        empty_artifact = artifact(joinpath(dir,"empty");raw=recorded_rows()[1:0,:],gate=DataFrame(pass_all=Bool[]))
        empty_overview = build_mode_overview([empty_artifact])
        @test empty_overview.machine_label == ["fallback-machine"]
        @test empty_overview.failed_rows == [0] && empty_overview.rows_raw == [0]
        p,s = reports(dir,[empty_artifact])
        @test occursin("Solver-success samples across modes: `0/0`",p)
        @test occursin("| auto | n/a | n/a | n/a | n/a |",s)
        absent = select(overview,Not([:machine_label,:hardware_class,:failed_rows]))
        p,_ = reports(dir,[a];overview=absent)
        @test !occursin("Hardware classes observed:",p)
        @test !occursin("Machine labels observed:",p)
        @test !occursin("Solver-success samples across modes:",p)
        @test isempty(build_mode_overview(ModeRunArtifacts[]))
    end
end

@testset "aggregate returns complete artifacts and preserves supplied entry rows" begin
    mktempdir() do dir
        first_pass = artifact(joinpath(dir,"pass1");entry_raw=entry_rows(1),bench=2.0,split=1.0,orbit=4.0,entry=3.0,total=13.0)
        second_pass = artifact(joinpath(dir,"pass2");entry_raw=entry_rows(2),gate=DataFrame(pass_all=[false]),bench=4.0,split=3.0,orbit=6.0,entry=7.0,total=27.0)
        runs = [ArmPassResult(pass=1,arm=ARM,artifact=first_pass),ArmPassResult(pass=2,arm=ARM,artifact=second_pass)]
        originals = deepcopy([r.artifact.raw_df for r in runs])
        original_entry = [read(r.artifact.entry_duration_raw_path) for r in runs]
        a = _aggregate_arm_artifacts(dir,MATRIX,static_config(dir),ARM,runs)
        @test a isa ModeRunArtifacts
        @test isequal([r.artifact.raw_df for r in runs],originals)
        @test [read(r.artifact.entry_duration_raw_path) for r in runs] == original_entry
        @test a.elapsed_s == 20.0
        @test a.entry_duration_elapsed_s == 5.0
        @test all(isfile,(a.raw_path,a.summary_path,a.orbit_raw_path,a.orbit_summary_path,
            a.entry_duration_raw_path,a.entry_duration_summary_path,a.stage_timing_path,a.report_path))
        report = read(a.report_path,String)
        @test occursin("Split rollout gate pass rows: `1/3`",report)
        @test occursin("Mean total elapsed: `20.0 s`",report)
        @test occursin("Mean entry-duration elapsed: `5.0 s`",report)
        for path in (a.entry_duration_raw_path,a.entry_duration_summary_path,
                     first_pass.entry_duration_raw_path,second_pass.entry_duration_raw_path)
            @test occursin(path,report)
        end
        gate = CSV.read(a.split_gate_csv_path,DataFrame)
        @test gate.pass_all == [true,false,false]
        @test gate.pass == [1,1,2]
        raw = CSV.read(a.raw_path,DataFrame)
        @test raw.solve_success == [false,true,false,true]
        @test raw.pass == [1,1,2,2]
        @test raw.arm == fill(ARM.label,4)
        summary = CSV.read(a.summary_path,DataFrame)
        @test only(summary.samples_failed) == 2 && only(summary.samples_success) == 2
        @test only(summary.total_time_mean_s) == 2.0
        @test only(summary.total_time_all_attempts_s) == 204.0
        @test _stage_elapsed_s(a.stage_timing_path) == (3.0,2.0,5.0,5.0,20.0)
        entry = CSV.read(a.entry_duration_raw_path,DataFrame)
        expected = vcat(entry_rows(1),entry_rows(2))
        @test isequal(select(entry,names(expected)),expected)
        @test entry.pass == vcat(fill(1,6),fill(2,6))
        @test entry.arm == fill(ARM.label,12)
        @test entry.policy_matrix == fill(String(MATRIX.key),12)
        @test entry.is_terminal_attempt == repeat([true,false,true,true,true,true],2)
        @test entry.attempt == repeat([1,1,2,1,1,2],2)
        entry_summary = CSV.read(a.entry_duration_summary_path,DataFrame)
        @test entry_summary.entry_atmospheric_interface_count == [1,1,3,3,5]
        @test entry_summary.entry_run_role == ["measured","reference","measured","reference","measured"]
        @test entry_summary.samples_total == [4,2,2,2,2]
        @test entry_summary.samples_success == [2,2,2,2,0]
        @test entry_summary.samples_failed == [2,0,0,0,2]
        measured = only(eachrow(filter(row -> row.entry_atmospheric_interface_count == 1 && row.entry_run_role == "measured",entry_summary)))
        @test measured.total_time_mean_s == 3.0
        @test measured.passage_duration_mean_s == 30.0
        @test measured.event_time_abs_error_mean_s == 5.0
        @test measured.event_time_abs_error_max_s == 7.0
        @test measured.reference_event_time_mean_s == 150.0
        @test measured.event_time_mean_s == 155.0
        failed = only(eachrow(filter(row -> row.entry_atmospheric_interface_count == 5,entry_summary)))
        @test failed.success_rate == 0.0
        for column in (:samples,:passage_duration_mean_s,:wall_time_per_passage_mean_s,
                       :event_time_abs_error_mean_s,:reference_event_time_mean_s,:total_time_mean_s)
            @test ismissing(failed[column])
        end
        overview = build_mode_overview([a])
        @test overview.orbit_share_pct == [25.0]
        @test overview.entry_duration_share_pct == [25.0]
        pipeline, static = reports(dir,[a];overview=overview)
        for path in (a.entry_duration_raw_path,a.entry_duration_summary_path)
            @test occursin(path,pipeline)
            @test occursin(path,static)
        end
    end
end

@testset "empty entry selections remain empty without invented measurements" begin
    for empty_entry in (entry_rows()[1:0,:],DataFrame())
        mktempdir() do dir
            saved = joinpath(dir,"saved")
            save_runtime_run(saved;entry_raw=empty_entry,hardware=false)
            loaded = _arm_result_artifacts(ARM,static_config(dir),saved)
            @test isempty(CSV.read(loaded.entry_duration_raw_path,DataFrame))
            @test isempty(CSV.read(loaded.entry_duration_summary_path,DataFrame))
            runs = [ArmPassResult(pass=1,arm=ARM,artifact=artifact(joinpath(dir,"pass1");entry_raw=empty_entry)),
                    ArmPassResult(pass=2,arm=ARM,artifact=artifact(joinpath(dir,"pass2");entry_raw=empty_entry))]
            a = _aggregate_arm_artifacts(dir,MATRIX,static_config(dir),ARM,runs)
            entry = CSV.read(a.entry_duration_raw_path,DataFrame)
            summary = CSV.read(a.entry_duration_summary_path,DataFrame)
            @test isempty(entry) && isempty(summary)
            @test all(column -> column in propertynames(entry),(:arm,:pass,:policy_matrix))
            @test a.entry_duration_elapsed_s == 3.0
            @test a.elapsed_s == 10.0
        end
    end
end

@testset "wrapper fixture has no native worker or simulation imports" begin
    for name in (:SpaceAGORA,:SimulationModel,:GRAMSuite,:SPICE,:Plots,:Distributed)
        @test !isdefined(@__MODULE__,name)
    end
end
end # module
