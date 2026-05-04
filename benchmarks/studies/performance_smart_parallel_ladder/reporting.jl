function _write_aggregate_rung_report(
    path::String,
    config::SmartLadderConfig,
    rung::LadderRungSpec,
    runs::Vector{LadderPassResult},
    bench_elapsed_s::Float64,
    orbit_elapsed_s::Float64,
    entry_duration_elapsed_s::Float64,
    total_elapsed_s::Float64;
    solver_label::Union{Nothing, String}=nothing,
    solver_mode::Union{Nothing, String}=nothing,
    split_gate_csv_path::Union{Nothing, String}=nothing,
    split_gate_pass_rows::Int=0,
    split_gate_total_rows::Int=0,
    multirate_gate_csv_path::Union{Nothing, String}=nothing,
    multirate_gate_pass_rows::Int=0,
    multirate_gate_total_rows::Int=0
)
    open(path, "w") do io
        println(io, "# Smart Parallel Ladder Aggregated Rung Report")
        println(io)
        println(io, "- Generated (UTC): $(now(UTC))")
        println(io, "- Profile: `$(config.profile.name)`")
        println(io, "- Rung: `$(rung.label)` (`$(rung.mode)`) ")
        println(io, "- Description: $(rung.description)")
        println(io, "- Matrix: `$(rung.matrix)`")
        println(io, "- Backend: `$(rung.backend)`")
        println(io, "- Inner adaptive: `$(rung.inner_adaptive)`")
        println(io, "- Outer-route adaptive: `$(rung.outer_route_adaptive)`")
        if !(solver_label === nothing)
            println(io, "- Solver label: `$(solver_label)`")
        end
        if !(solver_mode === nothing)
            println(io, "- Solver mode override: `$(solver_mode)`")
        end
        println(io, "- Aggregated passes: `$(length(runs))`")
        println(io, "- Mean run_benchmarks elapsed: `$(round(bench_elapsed_s; digits=3)) s`")
        println(io, "- Mean per-orbit elapsed: `$(round(orbit_elapsed_s; digits=3)) s`")
        println(io, "- Mean entry-duration elapsed: `$(round(entry_duration_elapsed_s; digits=3)) s`")
        println(io, "- Mean total elapsed: `$(round(total_elapsed_s; digits=3)) s`")
        if split_gate_total_rows > 0
            println(io, "- Split rollout gate pass rows: `$(split_gate_pass_rows)/$(split_gate_total_rows)`")
        else
            println(io, "- Split rollout gate rows: `0`")
        end
        if !(split_gate_csv_path === nothing)
            println(io, "- Split rollout gate aggregated CSV: `$(split_gate_csv_path)`")
        end
        if multirate_gate_total_rows > 0
            println(io, "- Multirate rollout gate pass rows: `$(multirate_gate_pass_rows)/$(multirate_gate_total_rows)`")
        else
            println(io, "- Multirate rollout gate rows: `0`")
        end
        if !(multirate_gate_csv_path === nothing)
            println(io, "- Multirate rollout gate aggregated CSV: `$(multirate_gate_csv_path)`")
        end
        println(io)
        println(io, "## Source Runs")
        println(io)
        for run in runs
            println(io, "- pass=$(run.pass), solver=$(run.solver_label): raw=`$(run.artifact.raw_path)`, per-orbit raw=`$(run.artifact.orbit_raw_path)`, report=`$(run.artifact.report_path)`")
        end
    end
end

function _aggregate_rung_artifacts(
    config::SmartLadderConfig,
    rung::LadderRungSpec,
    runs::Vector{LadderPassResult};
    solver_split::Bool=false,
    solver_label::Union{Nothing, String}=nothing,
    solver_mode::Union{Nothing, String}=nothing
)::ModeRunArtifacts
    isempty(runs) && error("No pass data collected for rung '$(rung.label)'.")

    raw_df = DataFrame()
    orbit_raw_df = DataFrame()
    entry_duration_raw_df = DataFrame()
    split_gate_df = DataFrame()
    multirate_gate_df = DataFrame()
    for run in runs
        run_solver_mode = run.solver_mode === nothing ? nothing : String(run.solver_mode)
        raw_df = vcat(
            raw_df,
            _tag_rung_column(run.artifact.raw_df, rung; pass_idx=run.pass, solver_label=run.solver_label, solver_mode=run_solver_mode);
            cols=:union
        )
        local_orbit_raw = CSV.read(run.artifact.orbit_raw_path, DataFrame)
        orbit_raw_df = vcat(
            orbit_raw_df,
            _tag_rung_column(local_orbit_raw, rung; pass_idx=run.pass, solver_label=run.solver_label, solver_mode=run_solver_mode);
            cols=:union
        )
        local_entry_duration_raw = CSV.read(run.artifact.entry_duration_raw_path, DataFrame)
        entry_duration_raw_df = vcat(
            entry_duration_raw_df,
            _tag_rung_column(local_entry_duration_raw, rung; pass_idx=run.pass, solver_label=run.solver_label, solver_mode=run_solver_mode);
            cols=:union
        )
        if !(run.artifact.split_gate_df === nothing)
            split_gate_df = vcat(
                split_gate_df,
                _tag_rung_column(run.artifact.split_gate_df, rung; pass_idx=run.pass, solver_label=run.solver_label, solver_mode=run_solver_mode);
                cols=:union
            )
        end
        run_outdir = dirname(run.artifact.raw_path)
        run_multirate_csv_path = _latest_artifact_path_optional(run_outdir, "multirate_rollout_gate", config.profile.name, ".csv")
        if !(run_multirate_csv_path === nothing)
            run_multirate_df = CSV.read(run_multirate_csv_path, DataFrame)
            multirate_gate_df = vcat(
                multirate_gate_df,
                _tag_rung_column(run_multirate_df, rung; pass_idx=run.pass, solver_label=run.solver_label, solver_mode=run_solver_mode);
                cols=:union
            )
        end
    end

    summary_df = summarize_results(raw_df)
    orbit_summary_df = summarize_per_orbit_results(orbit_raw_df)
    entry_duration_summary_df = summarize_entry_duration_results(entry_duration_raw_df)

    bench_elapsed_s = mean([run.artifact.bench_elapsed_s for run in runs])
    split_gate_elapsed_s = mean([run.artifact.split_gate_elapsed_s for run in runs])
    orbit_elapsed_s = mean([run.artifact.orbit_elapsed_s for run in runs])
    entry_duration_elapsed_s = mean([run.artifact.entry_duration_elapsed_s for run in runs])
    total_elapsed_s = bench_elapsed_s + split_gate_elapsed_s + orbit_elapsed_s + entry_duration_elapsed_s

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    agg_outdir = solver_split && !(solver_label === nothing) ?
        joinpath(config.outdir, "aggregate", "solver_$(solver_label)", rung.label) :
        joinpath(config.outdir, "aggregate", rung.label)
    mkpath(agg_outdir)

    raw_path = joinpath(agg_outdir, "runtime_raw_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    summary_path = joinpath(agg_outdir, "runtime_summary_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    orbit_raw_path = joinpath(agg_outdir, "runtime_per_orbit_raw_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    orbit_summary_path = joinpath(agg_outdir, "runtime_per_orbit_summary_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    entry_duration_raw_path = joinpath(agg_outdir, "runtime_entry_duration_raw_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    entry_duration_summary_path = joinpath(agg_outdir, "runtime_entry_duration_summary_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    stage_timing_path = joinpath(agg_outdir, "runtime_stage_timing_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    hardware_info_path = joinpath(agg_outdir, "runtime_hardware_info_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    split_gate_csv_path = nrow(split_gate_df) > 0 ? joinpath(agg_outdir, "split_rollout_gate_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv") : nothing
    multirate_gate_csv_path = nrow(multirate_gate_df) > 0 ? joinpath(agg_outdir, "multirate_rollout_gate_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv") : nothing
    report_path = joinpath(agg_outdir, "runtime_report_agg_$(config.profile.name)_$(rung.label)_$(stamp).md")

    split_gate_pass_rows = (nrow(split_gate_df) > 0 && _has_column(split_gate_df, :pass_all)) ? count(Bool.(split_gate_df.pass_all)) : 0
    split_gate_total_rows = nrow(split_gate_df)
    multirate_gate_pass_rows = (nrow(multirate_gate_df) > 0 && _has_column(multirate_gate_df, :pass_all)) ? count(Bool.(multirate_gate_df.pass_all)) : 0
    multirate_gate_total_rows = nrow(multirate_gate_df)

    hw = _runtime_hardware_snapshot()
    hardware_info_df = DataFrame([
        (
            profile=config.profile.name,
            rung=rung.label,
            rung_mode=String(rung.mode),
            matrix=String(rung.matrix),
            backend=rung.backend,
            inner_adaptive=rung.inner_adaptive,
            outer_route_adaptive=rung.outer_route_adaptive,
            solver_label=(solver_label === nothing ? missing : solver_label),
            solver_mode=(solver_mode === nothing ? missing : solver_mode),
            machine_label=hw.machine_label,
            hardware_class=hw.hardware_class,
            host_name=hw.host_name,
            cpu_name=hw.cpu_name,
            cpu_threads=hw.cpu_threads,
            julia_threads=hw.julia_threads,
            os=hw.os,
            arch=hw.arch
        )
    ])
    stage_timing_df = DataFrame(
        stage=["run_benchmarks", "run_split_rollout_gate", "run_per_orbit", "run_entry_duration_sweep", "total"],
        elapsed_s=[bench_elapsed_s, split_gate_elapsed_s, orbit_elapsed_s, entry_duration_elapsed_s, total_elapsed_s]
    )

    CSV.write(raw_path, raw_df)
    CSV.write(summary_path, summary_df)
    CSV.write(orbit_raw_path, orbit_raw_df)
    CSV.write(orbit_summary_path, orbit_summary_df)
    CSV.write(entry_duration_raw_path, entry_duration_raw_df)
    CSV.write(entry_duration_summary_path, entry_duration_summary_df)
    CSV.write(stage_timing_path, stage_timing_df)
    CSV.write(hardware_info_path, hardware_info_df)
    if !(split_gate_csv_path === nothing)
        CSV.write(split_gate_csv_path, split_gate_df)
    end
    if !(multirate_gate_csv_path === nothing)
        CSV.write(multirate_gate_csv_path, multirate_gate_df)
    end

    _write_aggregate_rung_report(
        report_path,
        config,
        rung,
        runs,
        bench_elapsed_s,
        orbit_elapsed_s,
        entry_duration_elapsed_s,
        total_elapsed_s;
        solver_label=solver_label,
        solver_mode=solver_mode,
        split_gate_csv_path=split_gate_csv_path,
        split_gate_pass_rows=split_gate_pass_rows,
        split_gate_total_rows=split_gate_total_rows,
        multirate_gate_csv_path=multirate_gate_csv_path,
        multirate_gate_pass_rows=multirate_gate_pass_rows,
        multirate_gate_total_rows=multirate_gate_total_rows
    )

    backend_label =
        "$(rung.backend)|matrix=$(rung.matrix)|inner=$(rung.inner_adaptive ? "adaptive" : "static")|outer=$(rung.outer_route_adaptive ? "adaptive" : "static")"
    return ModeRunArtifacts(
        mode=rung.mode,
        backend=backend_label,
        elapsed_s=total_elapsed_s,
        bench_elapsed_s=bench_elapsed_s,
        orbit_elapsed_s=orbit_elapsed_s,
        entry_duration_elapsed_s=entry_duration_elapsed_s,
        raw_path=raw_path,
        summary_path=summary_path,
        orbit_raw_path=orbit_raw_path,
        orbit_summary_path=orbit_summary_path,
        entry_duration_raw_path=entry_duration_raw_path,
        entry_duration_summary_path=entry_duration_summary_path,
        report_path=report_path,
        stage_timing_path=stage_timing_path,
        hardware_info_path=hardware_info_path,
        split_gate_elapsed_s=split_gate_elapsed_s,
        split_gate_csv_path=split_gate_csv_path,
        split_gate_df=(split_gate_total_rows > 0 ? split_gate_df : nothing),
        raw_df=raw_df,
        summary_df=summary_df,
        orbit_summary_df=orbit_summary_df
    )
end

function _build_solver_factorial_overview_df(
    config::SmartLadderConfig,
    rungs::Vector{LadderRungSpec},
    artifacts_by_solver::Dict{String, Vector{ModeRunArtifacts}},
    solver_mode_by_label::Dict{String, Union{Nothing, String}},
    primary_solver_label::String
)::DataFrame
    rows = DataFrame()
    rung_map = Dict(r.mode => r for r in rungs)
    solver_labels = sort(collect(keys(artifacts_by_solver)))
    for solver_label in solver_labels
        artifacts = artifacts_by_solver[solver_label]
        local_df = build_mode_overview(artifacts)
        local_df.mode = [
            haskey(rung_map, Symbol(mode_name)) ? rung_map[Symbol(mode_name)].label : String(mode_name)
            for mode_name in String.(local_df.mode)
        ]
        local_df[!, :solver_label] = fill(solver_label, nrow(local_df))
        local_df[!, :solver_mode] = fill(get(solver_mode_by_label, solver_label, nothing), nrow(local_df))
        local_df[!, :solver_primary] = fill(solver_label == primary_solver_label, nrow(local_df))
        rows = vcat(rows, local_df; cols=:union)
    end
    if nrow(rows) > 0
        sort!(rows, [:solver_label, :mode])
    end
    return rows
end

function _build_solver_factorial_comparison_df(
    artifacts_by_solver::Dict{String, Vector{ModeRunArtifacts}},
    solver_mode_by_label::Dict{String, Union{Nothing, String}},
    primary_solver_label::String
)::DataFrame
    rows = DataFrame()
    solver_labels = sort(collect(keys(artifacts_by_solver)))
    for solver_label in solver_labels
        artifacts = artifacts_by_solver[solver_label]
        local_df = build_comparison_table(artifacts)
        local_df[!, :solver_label] = fill(solver_label, nrow(local_df))
        local_df[!, :solver_mode] = fill(get(solver_mode_by_label, solver_label, nothing), nrow(local_df))
        local_df[!, :solver_primary] = fill(solver_label == primary_solver_label, nrow(local_df))
        rows = vcat(rows, local_df; cols=:union)
    end
    if nrow(rows) > 0
        sort!(rows, [:solver_label, :scenario])
    end
    return rows
end

function _write_smart_ladder_report(
    path::String,
    config::SmartLadderConfig,
    rungs::Vector{LadderRungSpec},
    overview_df::DataFrame,
    comparison_df::DataFrame,
    speedup_vs_r0_df::DataFrame,
    mission_family_df::DataFrame,
    layer_attribution_speedup_df::DataFrame,
    layer_attribution_family_df::DataFrame,
    layer_attribution_synergy_df::DataFrame,
    thermal_contribution_df::DataFrame,
    fidelity_df::DataFrame,
    accuracy_df::DataFrame,
    deep_accuracy_df::DataFrame,
    montecarlo_parity_df::DataFrame,
    route_mix_df::DataFrame,
    artifacts::Vector{ModeRunArtifacts},
    run_order_df::DataFrame,
    protocol_summary_df::DataFrame;
    solver_variants::Vector{NamedTuple{(:label, :solver_mode), Tuple{String, Union{Nothing, String}}}}=NamedTuple{(:label, :solver_mode), Tuple{String, Union{Nothing, String}}}[],
    primary_solver_label::String="",
    primary_solver_mode::Union{Nothing, String}=nothing,
    comparison_plot_paths::Vector{String}=String[],
    comparison_csv_path::String="",
    overview_csv_path::String="",
    speedup_vs_r0_csv_path::String="",
    mission_family_csv_path::String="",
    layer_attribution_speedup_csv_path::String="",
    layer_attribution_family_csv_path::String="",
    layer_attribution_synergy_csv_path::String="",
    thermal_contribution_csv_path::String="",
    fidelity_csv_path::String="",
    accuracy_csv_path::String="",
    deep_accuracy_csv_path::String="",
    montecarlo_parity_csv_path::String="",
    route_mix_csv_path::String="",
    merged_raw_csv_path::String="",
    merged_summary_csv_path::String="",
    run_order_csv_path::String="",
    protocol_summary_csv_path::String="",
    solver_factorial_overview_csv_path::String="",
    solver_factorial_comparison_csv_path::String=""
)
    generated = string(now(UTC))
    nthreads = Threads.nthreads()
    cpu_threads = Sys.CPU_THREADS
    rung_labels = join([r.label for r in rungs], ", ")

    open(path, "w") do io
        println(io, "# SpaceAGORA Smart Parallel Ladder")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Julia: `$(VERSION)`")
        println(io, "- Threads in wrapper process: `$nthreads`")
        println(io, "- Detected CPU threads: `$cpu_threads`")
        println(io, "- Profile: `$(config.profile.name)`")
        println(io, "- Rungs executed: `$rung_labels`")
        println(io, "- Passes: `$(config.passes)`")
        println(io, "- Randomized rung order: `$(config.randomize_rung_order)` (seed=`$(config.random_seed)`)")
        println(io, "- Outer-only backend: `$(config.outer_only_backend)`")
        println(io, "- Include layer-attribution matrix: `$(config.include_layer_attribution)`")
        println(io, "- Include control_stress in per-orbit stage: `$(config.include_control_stress_per_orbit)`")
        println(io, "- Control_stress full schedule override: warmup=`$(config.control_stress_warmup_full)`, repeats=`$(config.control_stress_repeats_full)`")
        println(io, "- Solver axis mode: `$(config.solver_axis)`")
        println(io, "- Primary solver label: `$(primary_solver_label)`")
        println(io, "- Primary solver mode override: `$(primary_solver_mode === nothing ? "inherit" : primary_solver_mode)`")
        if !isempty(solver_variants)
            variant_tokens = String[]
            for variant in solver_variants
                mode_token = variant.solver_mode === nothing ? "inherit" : String(variant.solver_mode)
                push!(variant_tokens, "$(variant.label):$(mode_token)")
            end
            println(io, "- Solver variants: `$(join(variant_tokens, ", "))`")
        end
        println(io)

        println(io, "## Ladder Definition")
        println(io)
        definition_df = DataFrame(
            rung=[r.label for r in rungs],
            mode=[String(r.mode) for r in rungs],
            matrix=[String(r.matrix) for r in rungs],
            backend=[r.backend for r in rungs],
            inner_adaptive=[r.inner_adaptive for r in rungs],
            outer_route_adaptive=[r.outer_route_adaptive for r in rungs],
            description=[r.description for r in rungs]
        )
        _write_markdown_table(io, definition_df)
        println(io)

        attribution_rungs = [r for r in rungs if _is_layer_attribution_mode(r.mode)]
        if !isempty(attribution_rungs)
            println(io, "## Layer Attribution Matrix")
            println(io)
            attribution_df = DataFrame(
                rung=[r.label for r in attribution_rungs],
                mode=[String(r.mode) for r in attribution_rungs],
                layer_kind=[_layer_kind_from_mode(r.mode) for r in attribution_rungs],
                layer_set=[_layer_set_from_mode(r.mode) for r in attribution_rungs],
                matrix=[String(r.matrix) for r in attribution_rungs],
                backend=[r.backend for r in attribution_rungs],
                description=[r.description for r in attribution_rungs]
            )
            _write_markdown_table(io, attribution_df)
            println(io)
        end

        println(io, "## Run Order By Pass")
        println(io)
        _write_markdown_table(io, run_order_df)
        println(io)

        println(io, "## Experimental Protocol Summary")
        println(io)
        println(io, "_Per-rung statistics summarize wrapper elapsed runtime across passes (mean/std/CV/95% CI). `start_mode` is derived from cache/state reset knobs (`clean`, `SPACEAGORA_PERF_OUTER_ROUTE_STATE_RESET`, `SPACEAGORA_PARALLEL_POLICY_STATE_RESET`)._")
        println(io)
        if nrow(protocol_summary_df) > 0
            protocol_columns = Symbol[
                :solver_label,
                :solver_mode,
                :solver_primary,
                :rung,
                :pass_count,
                :start_mode,
                :cache_mode,
                :total_elapsed_mean_s,
                :total_elapsed_std_s,
                :total_elapsed_cv_pct,
                :total_elapsed_ci95_low_s,
                :total_elapsed_ci95_high_s
            ]
            protocol_table = select(protocol_summary_df, protocol_columns)
            _write_markdown_table(io, protocol_table)
        else
            println(io, "- No protocol summary rows were produced.")
        end
        println(io)

        println(io, "## Mode Overview (Aggregated Across Passes)")
        println(io)
        _write_markdown_table(io, overview_df)
        println(io)

        println(io, "## Merged Comparison (Scenario Means + CI)")
        println(io)
        println(io, "_Note: `serial_total_time_mean_s` corresponds to rung `r0_true_serial`._")
        println(io)
        _write_markdown_table(io, comparison_df)
        println(io)

        println(io, "## Speedup vs R0 (Wall Time and Stages)")
        println(io)
        _write_markdown_table(io, speedup_vs_r0_df)
        println(io)

        println(io, "## Mission-Family Speedup Distribution vs R0")
        println(io)
        println(io, "_Family rules include dedicated buckets for `thermal_entry`, `srp_heavy`, `articulated_multibody`, `multi_sat_control`, `long_constellation`, and `effector_stress`; otherwise `mc_short`/`mc_long` by mission time <=/> 7200 s, `multi_sat` for satellites >= 2, `high_fidelity_nbody` for N-body/harmonics, `long` for mission_time > 7200 s, then `short_light`._")
        println(io)
        if nrow(mission_family_df) > 0
            _write_markdown_table(io, mission_family_df)
        else
            println(io, "- No mission-family speedup samples available.")
        end
        println(io)

        if nrow(layer_attribution_speedup_df) > 0
            println(io, "## Layer Attribution Speedup vs Outer-Only")
            println(io)
            _write_markdown_table(io, layer_attribution_speedup_df)
            println(io)
        end

        if nrow(layer_attribution_family_df) > 0
            println(io, "## Layer Attribution Mission-Family Speedup vs Outer-Only")
            println(io)
            _write_markdown_table(io, layer_attribution_family_df)
            println(io)
        end

        if nrow(layer_attribution_synergy_df) > 0
            println(io, "## Layer Pairwise Synergy")
            println(io)
            println(io, "_`synergy_ratio > 1` indicates pairwise performance better than the multiplicative independent prediction from the corresponding single-layer runs._")
            println(io)
            _write_markdown_table(io, layer_attribution_synergy_df)
            println(io)
        end

        println(io, "## Thermal Contribution vs R0")
        println(io)
        println(io, "_Thermal rows include scenarios tagged `thermal_stress` or scenario names containing `thermal`._")
        println(io)
        if nrow(thermal_contribution_df) > 0
            _write_markdown_table(io, thermal_contribution_df)
        else
            println(io, "- No thermal contribution samples available.")
        end
        println(io)

        println(io, "## Fidelity Parity")
        println(io)
        _write_markdown_table(io, fidelity_df)
        println(io)

        println(io, "## Accuracy Parity vs R0")
        println(io)
        println(io, "_Metrics are relative errors (%) on terminal primary-spacecraft norms (position, velocity, mass) against `r0_true_serial`._")
        println(io)
        if nrow(accuracy_df) > 0
            _write_markdown_table(io, accuracy_df)
        else
            println(io, "- No accuracy parity samples available.")
        end
        println(io)

        println(io, "## Deep Accuracy Parity vs R0")
        println(io)
        println(io, "_Trajectory metrics report RMS/max relative position/velocity error; event-time parity reports periapsis and atmosphere-interface timing error; control parity reports propellant and equivalent impulse parity from primary-spacecraft mass consumption; callback-state parity reports event-count agreement._")
        println(io)
        if nrow(deep_accuracy_df) > 0
            _write_markdown_table(io, deep_accuracy_df)
        else
            println(io, "- Deep parity probes disabled or no deep parity samples available.")
        end
        println(io)

        println(io, "## Monte Carlo Distribution Parity vs R0")
        println(io)
        println(io, "_Distribution parity uses scenario-level mean/std/p90 relative deltas and KS distance for final-state and event-time metrics._")
        println(io)
        if nrow(montecarlo_parity_df) > 0
            _write_markdown_table(io, montecarlo_parity_df)
        else
            println(io, "- No Monte Carlo distribution parity samples available.")
        end
        println(io)

        println(io, "## Route Mix")
        println(io)
        _write_markdown_table(io, route_mix_df)
        println(io)

        println(io, "## Plots")
        println(io)
        println(io, "- Plot artifacts generated: `$(length(comparison_plot_paths))`")
        for plot_path in comparison_plot_paths
            println(io, "- `$(plot_path)`")
        end
        println(io)

        println(io, "## Output Files")
        println(io)
        println(io, "- Canonical comparison/overview outputs use primary solver label: `$(primary_solver_label)`")
        println(io, "- Comparison CSV: `$(comparison_csv_path)`")
        println(io, "- Overview CSV: `$(overview_csv_path)`")
        println(io, "- Speedup vs R0 CSV: `$(speedup_vs_r0_csv_path)`")
        println(io, "- Mission-family speedup CSV: `$(mission_family_csv_path)`")
        println(io, "- Layer-attribution speedup CSV: `$(layer_attribution_speedup_csv_path)`")
        println(io, "- Layer-attribution mission-family CSV: `$(layer_attribution_family_csv_path)`")
        println(io, "- Layer-attribution synergy CSV: `$(layer_attribution_synergy_csv_path)`")
        println(io, "- Thermal contribution CSV: `$(thermal_contribution_csv_path)`")
        println(io, "- Fidelity parity CSV: `$(fidelity_csv_path)`")
        println(io, "- Accuracy parity CSV: `$(accuracy_csv_path)`")
        println(io, "- Deep accuracy parity CSV: `$(deep_accuracy_csv_path)`")
        println(io, "- Monte Carlo distribution parity CSV: `$(montecarlo_parity_csv_path)`")
        println(io, "- Route mix CSV: `$(route_mix_csv_path)`")
        println(io, "- Merged raw CSV (all passes): `$(merged_raw_csv_path)`")
        println(io, "- Merged summary CSV (all passes): `$(merged_summary_csv_path)`")
        println(io, "- Run order CSV: `$(run_order_csv_path)`")
        println(io, "- Protocol summary CSV: `$(protocol_summary_csv_path)`")
        if !isempty(solver_factorial_overview_csv_path)
            println(io, "- Solver factorial mode overview CSV: `$(solver_factorial_overview_csv_path)`")
        end
        if !isempty(solver_factorial_comparison_csv_path)
            println(io, "- Solver factorial comparison CSV: `$(solver_factorial_comparison_csv_path)`")
        end
        println(io, "- Per-rung metadata CSV: `pass_XX[/solver_<label>]/<rung>/smart_parallel_ladder_rung_metadata_<profile>_<rung>.csv`")
        println(io, "- Per-rung env CSV: `pass_XX[/solver_<label>]/<rung>/smart_parallel_ladder_rung_env_<profile>_<rung>.csv`")
        for artifact in artifacts
            println(io, "- Rung `$(artifact.mode)` aggregated raw: `$(artifact.raw_path)`")
            println(io, "- Rung `$(artifact.mode)` aggregated summary: `$(artifact.summary_path)`")
            println(io, "- Rung `$(artifact.mode)` aggregated per-orbit summary: `$(artifact.orbit_summary_path)`")
            if !isempty(artifact.stage_timing_path)
                println(io, "- Rung `$(artifact.mode)` aggregated stage timing: `$(artifact.stage_timing_path)`")
            end
            if !isempty(artifact.hardware_info_path)
                println(io, "- Rung `$(artifact.mode)` aggregated hardware info: `$(artifact.hardware_info_path)`")
            end
            if !(artifact.split_gate_csv_path === nothing)
                println(io, "- Rung `$(artifact.mode)` aggregated split gate CSV: `$(artifact.split_gate_csv_path)`")
            end
            println(io, "- Rung `$(artifact.mode)` aggregated report: `$(artifact.report_path)`")
        end
        println(io)

        println(io, "## Reproducibility")
        println(io, "```bash")
        process_suffix = isnothing(config.process_workers) ? "" : " --process-workers=$(config.process_workers)"
        rung_suffix = isempty(config.rung_filter) ? "" : " --rungs=$(join(config.rung_filter, ","))"
        project_flag = Base.shell_escape(SMART_LADDER_PROJECT)
        println(
            io,
            "JULIA_NUM_THREADS=$(Threads.nthreads()) julia --project=$(project_flag) benchmarks/studies/performance_smart_parallel_ladder.jl " *
            "--profile=$(config.profile.name) --outdir=$(config.outdir) --clean=1 --passes=$(config.passes) " *
            "--randomize-rung-order=$(config.randomize_rung_order ? 1 : 0) --seed=$(config.random_seed) " *
            "--outer-only-backend=$(config.outer_only_backend)$(process_suffix) " *
            "--layer-attribution=$(config.include_layer_attribution ? 1 : 0)$(rung_suffix) " *
            "--include-control-stress-per-orbit=$(config.include_control_stress_per_orbit ? 1 : 0) " *
            "--control-stress-repeats-full=$(config.control_stress_repeats_full) " *
            "--control-stress-warmup-full=$(config.control_stress_warmup_full) " *
            "--solver-axis=$(config.solver_axis) " *
            "--solver-mode=$(config.solver_mode) " *
            "--solver-factors=$(join(config.solver_factor_modes, ","))"
        )
        println(io, "```")
        println(io)
        println(io, "Cross-machine replay/merge helper:")
        println(io, "```bash")
        cross_outdir = joinpath(dirname(config.outdir), "smart_parallel_ladder_cross_machine")
        println(
            io,
            "julia --project=.AGORA benchmarks/studies/performance_smart_parallel_ladder_cross_machine.jl " *
            "--profile=$(config.profile.name) --outdir=$(cross_outdir) " *
            "--input=machine_a:/path/to/machine_a/ladder_outdir --input=machine_b:/path/to/machine_b/ladder_outdir"
        )
        println(io, "```")
    end
end

function main_smart_parallel_ladder()
    config = parse_smart_ladder_cli()
    if config.clean
        rm(config.outdir; recursive=true, force=true)
    end
    mkpath(config.outdir)

    rungs = _ladder_rungs(config)
    rng = MersenneTwister(config.random_seed)
    solver_variants = _solver_variants(config)
    primary_solver = _primary_solver_variant(config)
    primary_solver_label = primary_solver.label
    primary_solver_mode = primary_solver.solver_mode
    solver_split = length(solver_variants) > 1

    println("Smart parallel ladder profile: $(config.profile.name)")
    println("Outdir: $(config.outdir)")
    println("Project path: $(SMART_LADDER_PROJECT)")
    println("Clean outdir: $(config.clean)")
    println("Rungs: $(join([r.label for r in rungs], ", "))")
    println("Passes: $(config.passes)")
    println("Randomized rung order: $(config.randomize_rung_order) (seed=$(config.random_seed))")
    println("Outer-only backend: $(config.outer_only_backend)")
    println("Include layer-attribution matrix: $(config.include_layer_attribution)")
    println("Include control_stress in per-orbit stage: $(config.include_control_stress_per_orbit)")
    println("Control_stress full schedule override: warmup=$(config.control_stress_warmup_full), repeats=$(config.control_stress_repeats_full)")
    println("Solver axis: $(config.solver_axis)")
    println("Primary solver: label=$(primary_solver_label), mode=$(primary_solver_mode === nothing ? "inherit" : primary_solver_mode)")
    println("Solver variants: $(join(["$(v.label):$(v.solver_mode === nothing ? "inherit" : String(v.solver_mode))" for v in solver_variants], ", "))")
    println("Wrapper Threads.nthreads()=$(Threads.nthreads()), Sys.CPU_THREADS=$(Sys.CPU_THREADS)")

    pass_results = LadderPassResult[]
    order_rows = NamedTuple[]
    for pass_idx in 1:config.passes
        pass_outdir = joinpath(config.outdir, "pass_$(lpad(pass_idx, 2, '0'))")
        mkpath(pass_outdir)
        ordered = _ordered_rungs(rungs, rng, config.randomize_rung_order)
        order_line = join([r.label for r in ordered], " -> ")
        println("[smart-ladder] pass=$(pass_idx) rung-order=$(order_line)")
        for solver_variant in solver_variants
            solver_label = solver_variant.label
            solver_mode = solver_variant.solver_mode
            solver_outdir = solver_split ? joinpath(pass_outdir, "solver_$(solver_label)") : pass_outdir
            mkpath(solver_outdir)
            println(
                "[smart-ladder] pass=$(pass_idx) solver=$(solver_label) mode=$(solver_mode === nothing ? "inherit" : String(solver_mode))"
            )
            for (order_idx, rung) in enumerate(ordered)
                push!(order_rows, (
                    pass=pass_idx,
                    solver_label=solver_label,
                    solver_mode=(solver_mode === nothing ? missing : String(solver_mode)),
                    order_index=order_idx,
                    rung=rung.label,
                    mode=String(rung.mode),
                    matrix=String(rung.matrix),
                    backend=rung.backend,
                    inner_adaptive=rung.inner_adaptive,
                    outer_route_adaptive=rung.outer_route_adaptive
                ))
                artifact = run_rung(
                    rung,
                    config,
                    pass_idx,
                    solver_outdir;
                    solver_label=solver_label,
                    solver_mode=solver_mode
                )
                push!(
                    pass_results,
                    LadderPassResult(
                        pass=pass_idx,
                        solver_label=solver_label,
                        solver_mode=solver_mode,
                        rung=rung,
                        artifact=artifact
                    )
                )
            end
        end
    end

    solver_mode_by_label = Dict{String, Union{Nothing, String}}(variant.label => variant.solver_mode for variant in solver_variants)
    artifacts_by_solver = Dict{String, Vector{ModeRunArtifacts}}()
    for solver_variant in solver_variants
        solver_label = solver_variant.label
        solver_mode = solver_variant.solver_mode
        solver_runs = [result for result in pass_results if result.solver_label == solver_label]
        solver_artifacts = ModeRunArtifacts[]
        for rung in rungs
            rung_runs = [result for result in solver_runs if result.rung.label == rung.label]
            push!(
                solver_artifacts,
                _aggregate_rung_artifacts(
                    config,
                    rung,
                    rung_runs;
                    solver_split=solver_split,
                    solver_label=solver_label,
                    solver_mode=solver_mode
                )
            )
        end
        artifacts_by_solver[solver_label] = solver_artifacts
    end
    aggregated_artifacts = artifacts_by_solver[primary_solver_label]

    comparison_df = build_comparison_table(aggregated_artifacts)
    overview_df = build_mode_overview(aggregated_artifacts)
    rung_map = Dict(r.mode => r for r in rungs)
    rung_label_by_mode = Dict(r.mode => r.label for r in rungs)
    overview_df.mode = [
        haskey(rung_map, Symbol(mode_name)) ? rung_map[Symbol(mode_name)].label : String(mode_name)
        for mode_name in String.(overview_df.mode)
    ]
    speedup_vs_r0_df = _build_vs_r0_speedup_table(aggregated_artifacts, rung_label_by_mode)
    mission_family_df = _build_mission_family_speedup_table(aggregated_artifacts, rung_label_by_mode)
    layer_attribution_speedup_df = _build_layer_attribution_speedup_table(aggregated_artifacts, rung_label_by_mode)
    layer_attribution_family_df = _build_layer_attribution_mission_family_speedup_table(aggregated_artifacts, rung_label_by_mode)
    layer_attribution_synergy_df = _build_layer_pairwise_synergy_table(aggregated_artifacts)
    thermal_contribution_df = _build_thermal_contribution_table(aggregated_artifacts, rung_label_by_mode)
    fidelity_df = _build_fidelity_parity_table(aggregated_artifacts, rung_label_by_mode, config)
    accuracy_df = _build_accuracy_parity_table(aggregated_artifacts, rung_label_by_mode)
    deep_accuracy_df = _build_deep_accuracy_parity_table(aggregated_artifacts, rungs, rung_label_by_mode, config)
    montecarlo_parity_df = _build_montecarlo_distribution_parity_table(aggregated_artifacts, rung_label_by_mode)
    route_mix_df = _build_route_mix_table(aggregated_artifacts, rung_label_by_mode)

    merged_raw_df = DataFrame()
    merged_summary_df = DataFrame()
    for result in pass_results
        result_solver_mode = result.solver_mode === nothing ? nothing : String(result.solver_mode)
        merged_raw_df = vcat(
            merged_raw_df,
            _tag_rung_column(
                result.artifact.raw_df,
                result.rung;
                pass_idx=result.pass,
                solver_label=result.solver_label,
                solver_mode=result_solver_mode
            );
            cols=:union
        )
        merged_summary_df = vcat(
            merged_summary_df,
            _tag_rung_column(
                result.artifact.summary_df,
                result.rung;
                pass_idx=result.pass,
                solver_label=result.solver_label,
                solver_mode=result_solver_mode
            );
            cols=:union
        )
    end

    run_order_df = DataFrame(order_rows)
    protocol_summary_df = _build_protocol_summary_df(config, rungs, pass_results; primary_solver_label=primary_solver_label)
    solver_factorial_overview_df = DataFrame()
    solver_factorial_comparison_df = DataFrame()
    if solver_split
        solver_factorial_overview_df = _build_solver_factorial_overview_df(
            config,
            rungs,
            artifacts_by_solver,
            solver_mode_by_label,
            primary_solver_label
        )
        solver_factorial_comparison_df = _build_solver_factorial_comparison_df(
            artifacts_by_solver,
            solver_mode_by_label,
            primary_solver_label
        )
    end

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    comparison_path = joinpath(config.outdir, "smart_parallel_ladder_compare_summary_$(config.profile.name)_$(stamp).csv")
    overview_path = joinpath(config.outdir, "smart_parallel_ladder_mode_overview_$(config.profile.name)_$(stamp).csv")
    speedup_vs_r0_path = joinpath(config.outdir, "smart_parallel_ladder_speedup_vs_r0_$(config.profile.name)_$(stamp).csv")
    mission_family_path = joinpath(config.outdir, "smart_parallel_ladder_mission_family_speedup_$(config.profile.name)_$(stamp).csv")
    layer_attribution_speedup_path = joinpath(config.outdir, "smart_parallel_ladder_layer_attribution_speedup_$(config.profile.name)_$(stamp).csv")
    layer_attribution_family_path = joinpath(config.outdir, "smart_parallel_ladder_layer_attribution_mission_family_$(config.profile.name)_$(stamp).csv")
    layer_attribution_synergy_path = joinpath(config.outdir, "smart_parallel_ladder_layer_attribution_synergy_$(config.profile.name)_$(stamp).csv")
    thermal_contribution_path = joinpath(config.outdir, "smart_parallel_ladder_thermal_contribution_$(config.profile.name)_$(stamp).csv")
    fidelity_path = joinpath(config.outdir, "smart_parallel_ladder_fidelity_parity_$(config.profile.name)_$(stamp).csv")
    accuracy_path = joinpath(config.outdir, "smart_parallel_ladder_accuracy_parity_$(config.profile.name)_$(stamp).csv")
    deep_accuracy_path = joinpath(config.outdir, "smart_parallel_ladder_deep_accuracy_parity_$(config.profile.name)_$(stamp).csv")
    montecarlo_parity_path = joinpath(config.outdir, "smart_parallel_ladder_montecarlo_distribution_parity_$(config.profile.name)_$(stamp).csv")
    route_mix_path = joinpath(config.outdir, "smart_parallel_ladder_route_mix_$(config.profile.name)_$(stamp).csv")
    merged_raw_path = joinpath(config.outdir, "smart_parallel_ladder_raw_merged_$(config.profile.name)_$(stamp).csv")
    merged_summary_path = joinpath(config.outdir, "smart_parallel_ladder_summary_merged_$(config.profile.name)_$(stamp).csv")
    run_order_path = joinpath(config.outdir, "smart_parallel_ladder_run_order_$(config.profile.name)_$(stamp).csv")
    protocol_summary_path = joinpath(config.outdir, "smart_parallel_ladder_protocol_summary_$(config.profile.name)_$(stamp).csv")
    solver_factorial_overview_path = joinpath(
        config.outdir,
        "smart_parallel_ladder_solver_factorial_mode_overview_$(config.profile.name)_$(stamp).csv"
    )
    solver_factorial_comparison_path = joinpath(
        config.outdir,
        "smart_parallel_ladder_solver_factorial_compare_summary_$(config.profile.name)_$(stamp).csv"
    )
    report_path = joinpath(config.outdir, "smart_parallel_ladder_report_$(config.profile.name)_$(stamp).md")
    comparison_plot_paths = generate_pipeline_comparison_plots(
        config.outdir,
        config.profile,
        "smart_parallel_ladder_$(stamp)",
        overview_df,
        comparison_df
    )

    CSV.write(comparison_path, comparison_df)
    CSV.write(overview_path, overview_df)
    CSV.write(speedup_vs_r0_path, speedup_vs_r0_df)
    CSV.write(mission_family_path, mission_family_df)
    CSV.write(layer_attribution_speedup_path, layer_attribution_speedup_df)
    CSV.write(layer_attribution_family_path, layer_attribution_family_df)
    CSV.write(layer_attribution_synergy_path, layer_attribution_synergy_df)
    CSV.write(thermal_contribution_path, thermal_contribution_df)
    CSV.write(fidelity_path, fidelity_df)
    CSV.write(accuracy_path, accuracy_df)
    CSV.write(deep_accuracy_path, deep_accuracy_df)
    CSV.write(montecarlo_parity_path, montecarlo_parity_df)
    CSV.write(route_mix_path, route_mix_df)
    CSV.write(merged_raw_path, merged_raw_df)
    CSV.write(merged_summary_path, merged_summary_df)
    CSV.write(run_order_path, run_order_df)
    CSV.write(protocol_summary_path, protocol_summary_df)
    if solver_split
        CSV.write(solver_factorial_overview_path, solver_factorial_overview_df)
        CSV.write(solver_factorial_comparison_path, solver_factorial_comparison_df)
    end
    _write_smart_ladder_report(
        report_path,
        config,
        rungs,
        overview_df,
        comparison_df,
        speedup_vs_r0_df,
        mission_family_df,
        layer_attribution_speedup_df,
        layer_attribution_family_df,
        layer_attribution_synergy_df,
        thermal_contribution_df,
        fidelity_df,
        accuracy_df,
        deep_accuracy_df,
        montecarlo_parity_df,
        route_mix_df,
        aggregated_artifacts,
        run_order_df,
        protocol_summary_df;
        solver_variants=solver_variants,
        primary_solver_label=primary_solver_label,
        primary_solver_mode=primary_solver_mode,
        comparison_plot_paths=comparison_plot_paths,
        comparison_csv_path=comparison_path,
        overview_csv_path=overview_path,
        speedup_vs_r0_csv_path=speedup_vs_r0_path,
        mission_family_csv_path=mission_family_path,
        layer_attribution_speedup_csv_path=layer_attribution_speedup_path,
        layer_attribution_family_csv_path=layer_attribution_family_path,
        layer_attribution_synergy_csv_path=layer_attribution_synergy_path,
        thermal_contribution_csv_path=thermal_contribution_path,
        fidelity_csv_path=fidelity_path,
        accuracy_csv_path=accuracy_path,
        deep_accuracy_csv_path=deep_accuracy_path,
        montecarlo_parity_csv_path=montecarlo_parity_path,
        route_mix_csv_path=route_mix_path,
        merged_raw_csv_path=merged_raw_path,
        merged_summary_csv_path=merged_summary_path,
        run_order_csv_path=run_order_path,
        protocol_summary_csv_path=protocol_summary_path,
        solver_factorial_overview_csv_path=(solver_split ? solver_factorial_overview_path : ""),
        solver_factorial_comparison_csv_path=(solver_split ? solver_factorial_comparison_path : "")
    )

    println()
    println("Smart parallel ladder complete.")
    println("comparison summary: $(comparison_path)")
    println("mode overview: $(overview_path)")
    println("speedup vs r0: $(speedup_vs_r0_path)")
    println("mission-family speedup: $(mission_family_path)")
    println("layer-attribution speedup: $(layer_attribution_speedup_path)")
    println("layer-attribution mission-family: $(layer_attribution_family_path)")
    println("layer-attribution synergy: $(layer_attribution_synergy_path)")
    println("thermal contribution: $(thermal_contribution_path)")
    println("fidelity parity: $(fidelity_path)")
    println("accuracy parity: $(accuracy_path)")
    println("deep accuracy parity: $(deep_accuracy_path)")
    println("monte carlo distribution parity: $(montecarlo_parity_path)")
    println("route mix: $(route_mix_path)")
    println("merged raw: $(merged_raw_path)")
    println("merged summary: $(merged_summary_path)")
    println("run order: $(run_order_path)")
    println("protocol summary: $(protocol_summary_path)")
    if solver_split
        println("solver factorial mode overview: $(solver_factorial_overview_path)")
        println("solver factorial comparison: $(solver_factorial_comparison_path)")
    end
    println("comparison plots: $(length(comparison_plot_paths))")
    println("report: $(report_path)")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_smart_parallel_ladder()
end
