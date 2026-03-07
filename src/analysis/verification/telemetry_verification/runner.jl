function _run_simulation_dataframe(args::SimulationConfiguration, scenario_name::String, truth::AtmosphereTruthConfig, profile::Symbol)
    return mktempdir() do tmp
        cfg_run = SimulationConfiguration(
            file_paths=args.file_paths,
            simulation_settings=SimulationSettings(
                results=true,
                verbose=false,
                results_directory=tmp,
                generate_plots=false,
                generate_filenames=false,
                normalize=false,
                save_csv=true
            ),
            mission_configuration=args.mission_configuration,
            environment_model=args.environment_model,
            dynamics_model=args.dynamics_model,
            guidance_model=args.guidance_model,
            navigation_model=args.navigation_model,
            control_model=args.control_model,
            initial_time=args.initial_time,
            integration_tolerances=args.integration_tolerances
        )

        save_fields = _save_fields_for_study()
        base_maxiters = _telemetry_solver_maxiters(profile)

        function _run_once(maxiters::Int)
            solve_result = nothing
            elapsed_s = @elapsed begin
                solver_mode = _telemetry_solver_mode()
                withenv(
                    "SPACEAGORA_WARN_NORMALIZE" => "0",
                    "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
                    "SPACEAGORA_SOLVER_MODE" => solver_mode,
                    "SPACEAGORA_SOLVER_MAXITERS" => string(maxiters),
                    "SPACEAGORA_GRAM_OFFLINE_SURROGATE" => truth.gram_offline_surrogate,
                    "SPACEAGORA_GRAM_STATIC_GRID" => truth.gram_static_grid ? "on" : "off",
                    "SPACEAGORA_GRAM_TRACK_CACHE" => truth.gram_track_cache ? "on" : "off",
                    "SPACEAGORA_GRAM_GLOBAL_LOCK" => truth.gram_global_lock
                ) do
                    cd(tmp) do
                        solve_result = SimulationEngine.run_simulation(
                            cfg_run;
                            isolate_state=false,
                            save_fields=save_fields,
                            return_solution=true,
                            return_solver_metadata=true
                        )
                    end
                end
            end
            return solve_result, elapsed_s
        end

        retry_used = false
        maxiters_used = base_maxiters
        solve_result = nothing
        elapsed_s = 0.0
        try
            solve_result, elapsed_s = _run_once(base_maxiters)
        catch err
            if _is_maxiters_error(err)
                retry_used = true
                maxiters_used = _telemetry_solver_retry_maxiters(base_maxiters)
                @warn "Telemetry scenario '$scenario_name' hit MaxIters; retrying with larger SPACEAGORA_SOLVER_MAXITERS." base_maxiters=base_maxiters retry_maxiters=maxiters_used
                solve_result, elapsed_s = _run_once(maxiters_used)
            else
                rethrow(err)
            end
        end

        csv_path = joinpath(tmp, "simulation_results.csv")
        isfile(csv_path) || error("Missing simulation output for $scenario_name: $csv_path")
        results_df = CSV.read(csv_path, DataFrame)

        solver_trace = solve_result.solver_trace
        solver_sequence = isempty(solver_trace) ? "n/a" : join([meta.solver for meta in solver_trace], "->")
        solver_fallback_count = count(meta -> meta.fallback_used, solver_trace)
        fallback_triggers = [meta.trigger_retcode for meta in solver_trace if !(meta.trigger_retcode isa Missing)]
        solver_fallback_trigger = isempty(fallback_triggers) ? "n/a" : join(sort(unique(string.(fallback_triggers))), "|")
        solver_info = (
            solver_mode=String(solve_result.solver_mode),
            solver_sequence=solver_sequence,
            solver_fallback_used=solver_fallback_count > 0,
            solver_fallback_count=solver_fallback_count,
            solver_fallback_trigger=solver_fallback_trigger,
            solver_retcode=string(solve_result.solution.retcode),
            solver_maxiters=maxiters_used,
            solver_maxiters_retry_used=retry_used
        )
        return (results_df=results_df, elapsed_s=elapsed_s, solver_info=solver_info)
    end
end

function _initial_condition_from_time_aligned_telemetry(telemetry)
    sma_m = telemetry.sma_km * 1e3
    ra_m = sma_m * (1.0 + telemetry.ecc)
    rp_m = sma_m * (1.0 - telemetry.ecc)
    return InitialCondition(
        ra=ra_m,
        rp=rp_m,
        i=telemetry.inc_deg,
        ω=telemetry.aop_deg,
        Ω=telemetry.raan_deg,
        ν=telemetry.ta_deg
    )
end

function _run_single_scenario(cfg::OrbitEventsScenarioConfig, profile::Symbol)
    scenario_start = time()
    cal = cfg.calibration
    use_calibration = _calibration_active(cal, profile)

    final_is_quick = profile == :quick
    final_orbits = final_is_quick ? cfg.target_orbits_quick : cfg.target_orbits_full
    final_points = final_is_quick ? cfg.compare_points_quick : cfg.compare_points_full

    eval_profile = (use_calibration && cal.search_on_quick_subset && profile == :full) ? :quick : profile
    eval_is_quick = eval_profile == :quick
    eval_orbits = eval_is_quick ? cfg.target_orbits_quick : cfg.target_orbits_full
    eval_points = eval_is_quick ? cfg.compare_points_quick : cfg.compare_points_full

    cd_candidates = (use_calibration && cal.fit_cd_scale) ?
        _grid_values(cal.cd_scale_min, cal.cd_scale_max, cal.cd_scale_steps) : [1.0]
    cr_base = cfg.srp_cr
    cr_candidates = (use_calibration && cal.fit_cr && cfg.srp_enabled) ?
        _grid_values(cal.cr_min, cal.cr_max, cal.cr_steps) : [cr_base]

    best_cd = cd_candidates[1]
    best_cr = cr_candidates[1]
    best_score = Inf

    if use_calibration
        for cd_scale in cd_candidates, cr_value in cr_candidates
            args_eval = _make_orbit_args(cfg, eval_orbits; cd_scale=cd_scale, cr_override=cr_value)
            args_eval = _with_study_settings(args_eval; quick=eval_is_quick)
            eval_run = _run_simulation_dataframe(args_eval, cfg.name, cfg.atmosphere_truth, eval_profile)
            eval_df = eval_run.results_df
            eval_rows, eval_errors = _orbit_rows_errors(cfg, args_eval, eval_df, eval_points)
            if cal.fit_bias
                eval_bias = _estimate_event_biases(eval_errors, cal.bias_abs_max_km)
                eval_rows, eval_errors = _orbit_rows_errors(cfg, args_eval, eval_df, eval_points; bias_by_event=eval_bias)
            end
            score = _calibration_score(eval_rows, cal.objective)
            if score < best_score
                best_score = score
                best_cd = cd_scale
                best_cr = cr_value
            end
        end
    else
        best_score = NaN
    end

    args_final = _make_orbit_args(cfg, final_orbits; cd_scale=best_cd, cr_override=best_cr)
    args_final = _with_study_settings(args_final; quick=final_is_quick)
    final_run = _run_simulation_dataframe(args_final, cfg.name, cfg.atmosphere_truth, profile)
    final_df = final_run.results_df
    selected_runtime_s = final_run.elapsed_s
    solver_info = final_run.solver_info
    final_rows, final_errors = _orbit_rows_errors(cfg, args_final, final_df, final_points)
    bias_by_event = (use_calibration && cal.fit_bias) ? _estimate_event_biases(final_errors, cal.bias_abs_max_km) : Dict{String, Float64}()
    if !isempty(bias_by_event)
        final_rows, final_errors = _orbit_rows_errors(cfg, args_final, final_df, final_points; bias_by_event=bias_by_event)
    end
    final_score = use_calibration ? _calibration_score(final_rows, cal.objective) : NaN
    calibration_runtime_s = time() - scenario_start
    annotated_rows = _annotate_calibration_rows(
        final_rows,
        best_cd,
        best_cr,
        bias_by_event,
        final_score,
        selected_runtime_s,
        args_final.integration_tolerances.dt_max_orbit,
        calibration_runtime_s,
        use_calibration,
        solver_info
    )
    return annotated_rows, final_errors, calibration_runtime_s
end

function _run_single_scenario(cfg::TimeAlignedScenarioConfig, profile::Symbol)
    scenario_start = time()
    cal = cfg.calibration
    use_calibration = _calibration_active(cal, profile)

    final_is_quick = profile == :quick
    final_points = final_is_quick ? cfg.max_points_quick : cfg.max_points_full
    final_telemetry = _load_time_aligned_telemetry(cfg, final_points)
    final_ic = _initial_condition_from_time_aligned_telemetry(final_telemetry)
    final_mission_time_s = max(final_telemetry.time_s[end] - final_telemetry.time_s[1], 1.0)

    eval_profile = (use_calibration && cal.search_on_quick_subset && profile == :full) ? :quick : profile
    eval_is_quick = eval_profile == :quick
    eval_points = eval_is_quick ? cfg.max_points_quick : cfg.max_points_full
    eval_telemetry = _load_time_aligned_telemetry(cfg, eval_points)
    eval_ic = _initial_condition_from_time_aligned_telemetry(eval_telemetry)
    eval_mission_time_s = max(eval_telemetry.time_s[end] - eval_telemetry.time_s[1], 1.0)

    cd_candidates = (use_calibration && cal.fit_cd_scale) ?
        _grid_values(cal.cd_scale_min, cal.cd_scale_max, cal.cd_scale_steps) : [1.0]
    cr_base = cfg.srp_cr
    cr_candidates = (use_calibration && cal.fit_cr && cfg.srp_enabled) ?
        _grid_values(cal.cr_min, cal.cr_max, cal.cr_steps) : [cr_base]

    best_cd = cd_candidates[1]
    best_cr = cr_candidates[1]
    best_score = Inf

    if use_calibration
        for cd_scale in cd_candidates, cr_value in cr_candidates
            args_eval = _make_time_aligned_args(
                cfg,
                eval_mission_time_s,
                eval_ic;
                cd_scale=cd_scale,
                cr_override=cr_value
            )
            args_eval = _with_study_settings(args_eval; quick=eval_is_quick)
            eval_run = _run_simulation_dataframe(args_eval, cfg.name, cfg.atmosphere_truth, eval_profile)
            eval_df = eval_run.results_df
            eval_rows, eval_errors = _time_aligned_rows_errors(cfg, args_eval, eval_df, eval_telemetry)
            if cal.fit_bias
                eval_bias = _estimate_event_biases(eval_errors, cal.bias_abs_max_km)
                eval_rows, eval_errors = _time_aligned_rows_errors(cfg, args_eval, eval_df, eval_telemetry; bias_by_event=eval_bias)
            end
            score = _calibration_score(eval_rows, cal.objective)
            if score < best_score
                best_score = score
                best_cd = cd_scale
                best_cr = cr_value
            end
        end
    else
        best_score = NaN
    end

    args_final = _make_time_aligned_args(
        cfg,
        final_mission_time_s,
        final_ic;
        cd_scale=best_cd,
        cr_override=best_cr
    )
    args_final = _with_study_settings(args_final; quick=final_is_quick)
    final_run = _run_simulation_dataframe(args_final, cfg.name, cfg.atmosphere_truth, profile)
    final_df = final_run.results_df
    selected_runtime_s = final_run.elapsed_s
    solver_info = final_run.solver_info
    final_rows, final_errors = _time_aligned_rows_errors(cfg, args_final, final_df, final_telemetry)
    bias_by_event = (use_calibration && cal.fit_bias) ? _estimate_event_biases(final_errors, cal.bias_abs_max_km) : Dict{String, Float64}()
    if !isempty(bias_by_event)
        final_rows, final_errors = _time_aligned_rows_errors(cfg, args_final, final_df, final_telemetry; bias_by_event=bias_by_event)
    end
    final_score = use_calibration ? _calibration_score(final_rows, cal.objective) : NaN
    calibration_runtime_s = time() - scenario_start
    annotated_rows = _annotate_calibration_rows(
        final_rows,
        best_cd,
        best_cr,
        bias_by_event,
        final_score,
        selected_runtime_s,
        args_final.integration_tolerances.dt_max_orbit,
        calibration_runtime_s,
        use_calibration,
        solver_info
    )
    return annotated_rows, final_errors, calibration_runtime_s
end

function _run_verification(cfg::StudyConfig)::VerificationResult
    scenarios = _load_scenarios_from_manifest(cfg.manifest_path)

    summary_rows = NamedTuple[]
    error_tables = DataFrame[]
    wall_start = time()

    println("Telemetry Orbit Accuracy Study")
    println(@sprintf("profile=%s enforce=%s", String(cfg.profile), string(cfg.enforce)))
    println("manifest=$(cfg.manifest_path)")
    println("deterministic_mode=GRAM(per-scenario atmosphere_truth from manifest)")

    for sc in scenarios
        truth = sc.atmosphere_truth
        cal = sc.calibration
        println(
            "Running scenario $(sc.name) ... " *
            "(truth=$(truth.assumption_id), atmosphere=$(truth.atmosphere_model), " *
            "dataset=$(truth.atmosphere_dataset), weather=$(truth.space_weather_model), solar=$(truth.solar_flux_model), " *
            "seed=$(truth.gram_seed), " *
            "surrogate=$(truth.gram_offline_surrogate), static_grid=$(truth.gram_static_grid), " *
            "track_cache=$(truth.gram_track_cache), calibration=$(_calibration_active(cal, cfg.profile))$(_scenario_status_extra(sc)))"
        )
        row_batch, err_batch, elapsed_s = _run_single_scenario(sc, cfg.profile)
        println(@sprintf("COMPUTATIONAL TIME [%s/%s] = %.3f s", sc.name, String(cfg.profile), elapsed_s))

        for row in row_batch
            gates = _evaluate_thresholds(row, sc, cfg.profile)
            push!(summary_rows, merge(
                row,
                gates,
                (
                    profile=String(cfg.profile),
                    source_file=_source_file(sc, row.event),
                    axis_units=_axis_units(sc),
                    value_units=_value_units(sc, row.event),
                    orbit_altitude_mode=_orbit_altitude_mode(sc),
                    maneuver_count=_maneuver_count(sc),
                    simulation_runtime_s=elapsed_s,
                    timestamp_utc=string(now(UTC)),
                    atmosphere_truth_id=sc.atmosphere_truth.assumption_id,
                    atmosphere_model=sc.atmosphere_truth.atmosphere_model,
                    atmosphere_dataset=sc.atmosphere_truth.atmosphere_dataset,
                    space_weather_model=sc.atmosphere_truth.space_weather_model,
                    solar_flux_model=sc.atmosphere_truth.solar_flux_model,
                    gram_seed=sc.atmosphere_truth.gram_seed,
                    gram_perturbation_scales=join(sc.atmosphere_truth.gram_perturbation_scales, ","),
                    gram_offline_surrogate=sc.atmosphere_truth.gram_offline_surrogate,
                    gram_static_grid=sc.atmosphere_truth.gram_static_grid,
                    gram_track_cache=sc.atmosphere_truth.gram_track_cache,
                    gram_global_lock=sc.atmosphere_truth.gram_global_lock
                )
            ))
        end
        dt_by_event = Dict{String, Float64}()
        for row in row_batch
            dt_by_event[String(row.event)] = Float64(row.dt_max_orbit_s)
        end
        for err_df in err_batch
            nrow(err_df) == 0 && continue
            evt = String(err_df.event[1])
            err_df.dt_max_orbit_s = fill(get(dt_by_event, evt, NaN), nrow(err_df))
        end
        append!(error_tables, err_batch)
    end

    total_runtime_s = time() - wall_start
    println(@sprintf("TOTAL COMPUTATIONAL TIME [%s] = %.3f s", String(cfg.profile), total_runtime_s))

    summary_df = DataFrame(summary_rows)
    errors_df = isempty(error_tables) ? DataFrame() : vcat(error_tables...; cols=:union)
    if nrow(summary_df) > 0
        sort!(summary_df, [:scenario, :event])
        summary_df.total_runtime_s = fill(total_runtime_s, nrow(summary_df))
    end
    _append_display_metric_columns!(summary_df)
    _append_display_error_columns!(errors_df, summary_df)

    mkpath(dirname(cfg.out_summary))
    mkpath(dirname(cfg.out_errors))
    CSV.write(cfg.out_summary, summary_df)
    CSV.write(cfg.out_errors, errors_df)

    println("\nSummary:")
    show(summary_df, allrows=true, allcols=true)
    println("\n\nSaved:")
    println("  summary: $(cfg.out_summary)")
    println("  errors : $(cfg.out_errors)")
    plots_dir = ""
    if cfg.generate_plots
        plots_dir = _generate_plots(cfg.out_summary, cfg.out_errors, cfg.profile)
        println("  plots  : $(plots_dir)")
    else
        println("  plots  : disabled")
    end

    if cfg.enforce && nrow(summary_df) > 0
        failed = summary_df[summary_df.pass .== false, :]
        if nrow(failed) > 0
            println("\nThreshold failures:")
            show(failed, allrows=true, allcols=true)
            error("Telemetry orbit accuracy thresholds failed for $(nrow(failed)) row(s).")
        end
    end

    return VerificationResult(
        summary=summary_df,
        errors=errors_df,
        summary_path=cfg.out_summary,
        errors_path=cfg.out_errors,
        plots_dir=plots_dir,
        profile=cfg.profile,
        enforce=cfg.enforce,
        total_runtime_s=total_runtime_s
    )
end

"""
    run_verification(request::VerificationRequest) -> VerificationResult

Run the telemetry verification study using an explicit typed request.
"""
function run_verification(request::VerificationRequest)::VerificationResult
    return _run_verification(_study_config(request))
end

"""
    run_verification_cli([args=copy(ARGS)]) -> VerificationResult

CLI-oriented verification entrypoint that parses command-line arguments and runs
the telemetry verification study.
"""
function run_verification_cli(args::Vector{String}=copy(ARGS))::VerificationResult
    cfg = parse_cli(args)
    request = _request_from_study_config(cfg)
    request = VerificationRequest(
        profile=request.profile,
        out_summary=request.out_summary,
        out_errors=request.out_errors,
        manifest_path=request.manifest_path,
        enforce=cfg.enforce,
        generate_plots=cfg.generate_plots
    )
    return run_verification(request)
end

"""
    run_study()

Convenience wrapper used by script-style entrypoints. It runs the telemetry
verification CLI flow and returns the summary and error tables.
"""
function run_study()
    result = run_verification_cli(copy(ARGS))
    return (summary=result.summary, errors=result.errors)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_verification_cli(copy(ARGS))
end

