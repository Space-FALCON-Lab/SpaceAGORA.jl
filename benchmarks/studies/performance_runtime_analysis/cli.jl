function parse_cli()
    profile_name = lowercase(get(ENV, "SPACEAGORA_PERF_PROFILE", "quick"))
    outdir = get(ENV, "SPACEAGORA_PERF_OUTDIR", DEFAULT_OUTPUT_DIR)

    for arg in ARGS
        if arg == "smoke"
            profile_name = "quick"
            ENV["SPACEAGORA_PERF_SMOKE"] = "1"
        elseif arg in ("quick", "full")
            profile_name = arg
        elseif startswith(arg, "--profile=")
            profile_name = lowercase(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]
        else
            throw(ArgumentError("Unknown argument '$arg'. Supported: [quick|full|smoke], --profile=..., --outdir=..."))
        end
    end

    return _profile_from_name(profile_name), abspath(outdir)
end

function main()
    spec, outdir = parse_cli()
    mkpath(outdir)
    outer_route_state = _load_outer_route_history!(spec, outdir)

    solver_mode_default = _perf_default_solver_mode(spec.name)
    solver_mode_effective = _perf_solver_mode_env(spec.name)
    sat_threshold = _priority_inner_sat_threshold()
    link_threshold = _priority_inner_link_threshold()
    outer_light_sat = _priority_outer_light_sat_threshold()
    outer_light_link = _priority_outer_light_link_threshold()
    outer_light_mission_s = _priority_outer_light_mission_threshold_s()

    println("Performance runtime analysis profile: $(spec.name)")
    println("Output directory: $outdir")
    println("Smoke mode: $(_perf_smoke_mode() ? "on" : "off")")
    println("Solver mode default: $(solver_mode_default)")
    println("Solver mode effective: $(solver_mode_effective)")
    println(
        "Outer-route adaptive=$(_outer_route_adaptive_enabled() ? "on" : "off"), " *
        "min_samples=$(_outer_route_min_samples()), " *
        "mc_process_min_samples=$(_outer_route_mc_process_min_samples()), " *
        "mc_process_min_mission_s=$(round(_outer_route_mc_process_min_mission_s(); digits=3))"
    )
    if outer_route_state.persist
        if outer_route_state.reset_requested
            println("Outer-route state cache reset requested; starting from empty history.")
        elseif outer_route_state.loaded_rows > 0
            println(
                "Outer-route state cache loaded rows=$(outer_route_state.loaded_rows), " *
                "signatures=$(outer_route_state.loaded_signatures) from $(outer_route_state.path)"
            )
        else
            println("Outer-route state cache path: $(outer_route_state.path) (no prior state loaded)")
        end
    else
        println("Outer-route state cache persistence: off")
    end
    println(
        "Priority thresholds: inner_sat=$(sat_threshold), inner_link=$(link_threshold), " *
        "outer_light_sat=$(outer_light_sat), outer_light_link=$(outer_light_link), " *
        "outer_light_mission_s=$(round(outer_light_mission_s; digits=3))"
    )

    planet = Earth("", SPICE_PATH)
    cases = build_cases(spec, planet)
    bench_started_ns = time_ns()
    raw_df = run_benchmarks(spec, cases, planet)
    summary_df = summarize_results(raw_df)
    density_backend_breakdown_df = summarize_density_backend_breakdown(raw_df)
    bench_elapsed_s = (time_ns() - bench_started_ns) / 1e9

    split_gate_df = nothing
    split_gate_csv_path = nothing
    split_gate_report_path = nothing
    split_gate_elapsed_s = 0.0
    if _split_rollout_enabled()
        split_gate_started_ns = time_ns()
        split_gate_result = evaluate_split_rollout_gate(spec, cases, outdir)
        split_gate_elapsed_s = (time_ns() - split_gate_started_ns) / 1e9
        split_gate_df = split_gate_result.df
        split_gate_csv_path = split_gate_result.csv_path
        split_gate_report_path = split_gate_result.report_path
    end

    multirate_gate_df = nothing
    multirate_gate_csv_path = nothing
    multirate_gate_report_path = nothing
    multirate_gate_elapsed_s = 0.0
    if _multirate_rollout_enabled()
        multirate_gate_started_ns = time_ns()
        multirate_gate_result = evaluate_multirate_rollout_gate(spec, cases, outdir)
        multirate_gate_elapsed_s = (time_ns() - multirate_gate_started_ns) / 1e9
        multirate_gate_df = multirate_gate_result.df
        multirate_gate_csv_path = multirate_gate_result.csv_path
        multirate_gate_report_path = multirate_gate_result.report_path
    end

    orbit_started_ns = time_ns()
    orbit_raw_df = run_per_orbit_for_scenarios(spec, cases, planet)
    orbit_summary_df = summarize_per_orbit_results(orbit_raw_df)
    orbit_elapsed_s = (time_ns() - orbit_started_ns) / 1e9
    entry_duration_started_ns = time_ns()
    entry_duration_raw_df = run_entry_duration_sweep(spec, cases, planet)
    entry_duration_summary_df = summarize_entry_duration_results(entry_duration_raw_df)
    entry_duration_elapsed_s = (time_ns() - entry_duration_started_ns) / 1e9
    total_elapsed_s = bench_elapsed_s + split_gate_elapsed_s + multirate_gate_elapsed_s + orbit_elapsed_s + entry_duration_elapsed_s
    outer_route_state_saved = _save_outer_route_history!(spec, outer_route_state.path, outer_route_state.persist)

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    raw_path = joinpath(outdir, "runtime_raw_$(spec.name)_$(stamp).csv")
    summary_path = joinpath(outdir, "runtime_summary_$(spec.name)_$(stamp).csv")
    orbit_raw_path = joinpath(outdir, "runtime_per_orbit_raw_$(spec.name)_$(stamp).csv")
    orbit_summary_path = joinpath(outdir, "runtime_per_orbit_summary_$(spec.name)_$(stamp).csv")
    entry_duration_raw_path = joinpath(outdir, "runtime_entry_duration_raw_$(spec.name)_$(stamp).csv")
    entry_duration_summary_path = joinpath(outdir, "runtime_entry_duration_summary_$(spec.name)_$(stamp).csv")
    stage_timing_path = joinpath(outdir, "runtime_stage_timing_$(spec.name)_$(stamp).csv")
    hardware_info_path = joinpath(outdir, "runtime_hardware_info_$(spec.name)_$(stamp).csv")
    inner_hint_layer_path = joinpath(outdir, "runtime_inner_hint_layers_$(spec.name)_$(stamp).csv")
    density_backend_breakdown_path = joinpath(outdir, "runtime_density_backend_breakdown_$(spec.name)_$(stamp).csv")
    report_path = joinpath(outdir, "runtime_report_$(spec.name)_$(stamp).md")
    hw = _runtime_hardware_snapshot()
    inner_hint_layer_df = _inner_hint_layer_report_df(spec, hw)
    plot_paths = generate_runtime_plots(
        outdir,
        spec,
        stamp,
        raw_df,
        summary_df,
        orbit_summary_df,
        entry_duration_summary_df;
        split_gate_df=split_gate_df,
        multirate_gate_df=multirate_gate_df,
        inner_hint_layer_df=inner_hint_layer_df,
        density_backend_breakdown_df=density_backend_breakdown_df
    )

    stage_names = ["run_benchmarks"]
    stage_elapsed = [bench_elapsed_s]
    if _split_rollout_enabled()
        push!(stage_names, "run_split_rollout_gate")
        push!(stage_elapsed, split_gate_elapsed_s)
    end
    if _multirate_rollout_enabled()
        push!(stage_names, "run_multirate_rollout_gate")
        push!(stage_elapsed, multirate_gate_elapsed_s)
    end
    push!(stage_names, "run_per_orbit")
    push!(stage_names, "run_entry_duration_sweep")
    push!(stage_names, "total")
    push!(stage_elapsed, orbit_elapsed_s)
    push!(stage_elapsed, entry_duration_elapsed_s)
    push!(stage_elapsed, total_elapsed_s)
    stage_timing_df = DataFrame(stage=stage_names, elapsed_s=stage_elapsed)
    hardware_info_df = DataFrame([
        (
            profile=spec.name,
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

    CSV.write(raw_path, raw_df)
    CSV.write(summary_path, summary_df)
    CSV.write(orbit_raw_path, orbit_raw_df)
    CSV.write(orbit_summary_path, orbit_summary_df)
    CSV.write(entry_duration_raw_path, entry_duration_raw_df)
    CSV.write(entry_duration_summary_path, entry_duration_summary_df)
    CSV.write(stage_timing_path, stage_timing_df)
    CSV.write(hardware_info_path, hardware_info_df)
    CSV.write(inner_hint_layer_path, inner_hint_layer_df)
    CSV.write(density_backend_breakdown_path, density_backend_breakdown_df)
    write_report(
        report_path,
        spec,
        raw_df,
        summary_df,
        orbit_summary_df,
        entry_duration_summary_df;
        plot_paths=plot_paths,
        split_gate_df=split_gate_df,
        split_gate_csv_path=split_gate_csv_path,
        split_gate_report_path=split_gate_report_path,
        multirate_gate_df=multirate_gate_df,
        multirate_gate_csv_path=multirate_gate_csv_path,
        multirate_gate_report_path=multirate_gate_report_path,
        stage_timing_df=stage_timing_df,
        inner_hint_layer_df=inner_hint_layer_df,
        inner_hint_layer_csv_path=inner_hint_layer_path,
        density_backend_breakdown_df=density_backend_breakdown_df,
        density_backend_breakdown_csv_path=density_backend_breakdown_path,
        entry_duration_summary_csv_path=entry_duration_summary_path
    )

    println("Analysis complete.")
    println("Raw results: $raw_path")
    println("Summary: $summary_path")
    println("Mission-time-sweep raw: $orbit_raw_path")
    println("Mission-time-sweep summary: $orbit_summary_path")
    println("Entry-duration-sweep raw: $entry_duration_raw_path")
    println("Entry-duration-sweep summary: $entry_duration_summary_path")
    println("Stage timing: $stage_timing_path")
    println("Hardware info: $hardware_info_path")
    println("Inner hint layers: $inner_hint_layer_path")
    println("Density backend breakdown: $density_backend_breakdown_path")
    if !(split_gate_df === nothing)
        pass_count = (:pass_all in names(split_gate_df)) ? count(Bool.(split_gate_df.pass_all)) : 0
        println("Split rollout gate: $(pass_count)/$(nrow(split_gate_df)) pass")
    end
    if !(split_gate_csv_path === nothing)
        println("Split rollout gate CSV: $(split_gate_csv_path)")
    end
    if !(split_gate_report_path === nothing)
        println("Split rollout gate report: $(split_gate_report_path)")
    end
    if !(multirate_gate_df === nothing)
        pass_count = (:pass_all in names(multirate_gate_df)) ? count(Bool.(multirate_gate_df.pass_all)) : 0
        println("Multirate rollout gate: $(pass_count)/$(nrow(multirate_gate_df)) pass")
    end
    if !(multirate_gate_csv_path === nothing)
        println("Multirate rollout gate CSV: $(multirate_gate_csv_path)")
    end
    if !(multirate_gate_report_path === nothing)
        println("Multirate rollout gate report: $(multirate_gate_report_path)")
    end
    if outer_route_state.persist
        println(
            "Outer-route state cache: $(outer_route_state_saved.path) " *
            "(rows_saved=$(outer_route_state_saved.rows), signatures=$(outer_route_state_saved.signatures))"
        )
    end
    println("Plots generated: $(length(plot_paths))")
    println("Report: $report_path")
end
