function _candidate_paths(outdir::String, stage::String, id::Int; draw::Union{Nothing, Int}=nothing)
    tag = draw === nothing ? @sprintf("cand_%04d_%s", id, stage) : @sprintf("cand_%04d_%s_draw%03d", id, stage, draw)
    manifests_dir = joinpath(outdir, "manifests")
    results_dir = joinpath(outdir, "results")
    logs_dir = joinpath(outdir, "logs")
    mkpath(manifests_dir)
    mkpath(results_dir)
    mkpath(logs_dir)
    return (
        manifest=joinpath(manifests_dir, "$(tag).toml"),
        summary=joinpath(results_dir, "$(tag)_summary.csv"),
        errors=joinpath(results_dir, "$(tag)_errors.csv"),
        log=joinpath(logs_dir, "$(tag).log")
    )
end

function _apply_uncertainty_to_scenario!(
    sc::Dict{String, Any},
    rng::AbstractRNG,
    cfg::TunerConfig
)
    if haskey(sc, "atmosphere_truth")
        at = Dict{String, Any}(sc["atmosphere_truth"])
        base_seed = Int(get(at, "gram_seed", 1001))
        at["gram_seed"] = base_seed + rand(rng, 1:1_000_000)

        base_scales = if haskey(at, "gram_perturbation_scales")
            [Float64(v) for v in at["gram_perturbation_scales"]]
        else
            [0.0, 0.0, 0.0, 0.0]
        end
        jitter = cfg.uncertainty_atm_scale
        pert = Float64[]
        for s in base_scales
            sigma = s > 0.0 ? jitter * s : jitter
            push!(pert, max(0.0, s + sigma * randn(rng)))
        end
        at["gram_perturbation_scales"] = pert
        sc["atmosphere_truth"] = at
    end

    frac = cfg.uncertainty_ic_scale
    if haskey(sc, "ra_m")
        sc["ra_m"] = _clamp_positive(Float64(sc["ra_m"]) * (1.0 + frac * randn(rng)))
    end
    if haskey(sc, "rp_altitude_m")
        sc["rp_altitude_m"] = _clamp_positive(Float64(sc["rp_altitude_m"]) + 1500.0 * frac * randn(rng))
    end
    if haskey(sc, "i_deg")
        sc["i_deg"] = clamp(Float64(sc["i_deg"]) + 0.5 * frac * randn(rng), 0.0, 180.0)
    end
    if haskey(sc, "aop_deg")
        sc["aop_deg"] = _wrap_0_360(Float64(sc["aop_deg"]) + 2.0 * frac * randn(rng))
    end
    if haskey(sc, "raan_deg")
        sc["raan_deg"] = _wrap_0_360(Float64(sc["raan_deg"]) + 2.0 * frac * randn(rng))
    end
    if haskey(sc, "ta_deg")
        sc["ta_deg"] = _wrap_0_360(Float64(sc["ta_deg"]) + 5.0 * frac * randn(rng))
    end

    return sc
end

function _apply_candidate_to_manifest(
    base_manifest::Dict{String, Any},
    cand::TuneCandidate,
    cfg::TunerConfig;
    uncertainty_draw::Union{Nothing, Int}=nothing
)::Dict{String, Any}
    doc = deepcopy(base_manifest)
    haskey(doc, "scenarios") || throw(ArgumentError("Manifest missing top-level scenarios array."))
    scenarios = doc["scenarios"]
    scenarios isa AbstractVector || throw(ArgumentError("Manifest scenarios must be an array."))

    shift_s = Int(cand.values["epoch_shift_s"])

    for i in eachindex(scenarios)
        sc = Dict{String, Any}(scenarios[i])

        if haskey(sc, "initial_time")
            init = Dict{String, Any}(sc["initial_time"])
            sec_raw = Float64(get(init, "second", 0.0))
            sec_int = floor(Int, sec_raw)
            sec_frac = sec_raw - sec_int
            sec_ms = round(Int, sec_frac * 1000)
            dt0 = DateTime(
                Int(init["year"]),
                Int(init["month"]),
                Int(init["day"]),
                Int(init["hour"]),
                Int(init["minute"]),
                sec_int,
                sec_ms
            )
            dt_shifted = dt0 + Second(shift_s)
            init["year"] = year(dt_shifted)
            init["month"] = month(dt_shifted)
            init["day"] = day(dt_shifted)
            init["hour"] = hour(dt_shifted)
            init["minute"] = minute(dt_shifted)
            init["second"] = second(dt_shifted) + millisecond(dt_shifted) / 1000
            sc["initial_time"] = init
        end

        if haskey(sc, "ra_m")
            sc["ra_m"] = _clamp_positive(Float64(sc["ra_m"]) * Float64(cand.values["ra_scale"]))
        end
        if haskey(sc, "rp_altitude_m")
            sc["rp_altitude_m"] = _clamp_positive(Float64(sc["rp_altitude_m"]) + Float64(cand.values["rp_altitude_offset_m"]))
        end
        if haskey(sc, "i_deg")
            sc["i_deg"] = clamp(Float64(sc["i_deg"]) + Float64(cand.values["i_offset_deg"]), 0.0, 180.0)
        end
        if haskey(sc, "aop_deg")
            sc["aop_deg"] = _wrap_0_360(Float64(sc["aop_deg"]) + Float64(cand.values["aop_offset_deg"]))
        end
        if haskey(sc, "raan_deg")
            sc["raan_deg"] = _wrap_0_360(Float64(sc["raan_deg"]) + Float64(cand.values["raan_offset_deg"]))
        end
        if haskey(sc, "ta_deg")
            sc["ta_deg"] = _wrap_0_360(Float64(sc["ta_deg"]) + Float64(cand.values["ta_offset_deg"]))
        end

        if haskey(sc, "spacecraft")
            spacecraft = Dict{String, Any}(sc["spacecraft"])
            haskey(spacecraft, "bus_mass_kg") && (spacecraft["bus_mass_kg"] = _clamp_positive(Float64(spacecraft["bus_mass_kg"]) * Float64(cand.values["bus_mass_scale"])))
            haskey(spacecraft, "prop_mass_kg") && (spacecraft["prop_mass_kg"] = _clamp_positive(Float64(spacecraft["prop_mass_kg"]) * Float64(cand.values["prop_mass_scale"])))
            haskey(spacecraft, "panel_mass_each_kg") && (spacecraft["panel_mass_each_kg"] = _clamp_positive(Float64(spacecraft["panel_mass_each_kg"]) * Float64(cand.values["panel_mass_scale"])))
            haskey(spacecraft, "panel_offset_y_m") && (spacecraft["panel_offset_y_m"] = _clamp_positive(Float64(spacecraft["panel_offset_y_m"]) * Float64(cand.values["panel_offset_scale"])))

            if haskey(spacecraft, "bus_dims_m")
                spacecraft["bus_dims_m"] = [Float64(v) * Float64(cand.values["bus_dims_scale"]) for v in spacecraft["bus_dims_m"]]
            end
            if haskey(spacecraft, "panel_dims_m")
                spacecraft["panel_dims_m"] = [Float64(v) * Float64(cand.values["panel_dims_scale"]) for v in spacecraft["panel_dims_m"]]
            end
            sc["spacecraft"] = spacecraft
        end

        haskey(sc, "srp_cr") && (sc["srp_cr"] = _clamp_positive(Float64(sc["srp_cr"]) * Float64(cand.values["srp_cr_scale"])))
        haskey(sc, "srp_area_m2") && (sc["srp_area_m2"] = _clamp_positive(Float64(sc["srp_area_m2"]) * Float64(cand.values["srp_area_scale"])))

        gravity_variant = String(cand.values["gravity_variant"])
        if gravity_variant == "j2_only"
            sc["gravity_model"] = "inverse_squared_j2"
            sc["gravity_harmonics_degree"] = 0
            sc["gravity_harmonics_order"] = 0
        else
            degree = if gravity_variant == "harm_deg2"
                2
            elseif gravity_variant == "harm_deg4"
                4
            elseif gravity_variant == "harm_deg8"
                8
            else
                20
            end
            sc["gravity_model"] = "inverse_squared"
            sc["gravity_harmonics_degree"] = degree
            sc["gravity_harmonics_order"] = degree
        end

        if haskey(sc, "maneuvers")
            mv = Dict{String, Any}(sc["maneuvers"])
            if haskey(mv, "orbit_numbers") && haskey(mv, "delta_v_mps")
                orbits = [Int(v) for v in mv["orbit_numbers"]]
                dvs = [Float64(v) for v in mv["delta_v_mps"]]
                tuned = similar(dvs)
                for j in eachindex(dvs)
                    scale = Float64(cand.values["dv_global_scale"])
                    if orbits[j] <= 50
                        scale *= Float64(cand.values["dv_early_scale"])
                    end
                    dv = dvs[j] * scale
                    if orbits[j] == 7
                        dv += Float64(cand.values["dv_orbit7_bias_mps"])
                    end
                    tuned[j] = dv
                end
                mv["delta_v_mps"] = tuned
                sc["maneuvers"] = mv
            end
        end

        if haskey(sc, "calibration")
            cal = Dict{String, Any}(sc["calibration"])
            cal["enabled"] = cfg.calibration_enabled
            sc["calibration"] = cal
        end

        if uncertainty_draw !== nothing
            rng = MersenneTwister(hash((cfg.seed, "uncertainty", String(get(sc, "name", "scenario")), cand.id, uncertainty_draw)))
            _apply_uncertainty_to_scenario!(sc, rng, cfg)
        end

        scenarios[i] = sc
    end

    doc["scenarios"] = scenarios
    return doc
end

@inline function _event_metric(df::DataFrame, scenario::String, event::String, col::Symbol)
    mask = (df.scenario .== scenario) .& (df.event .== event)
    idx = findfirst(mask)
    if idx === nothing
        return missing
    end
    value = df[idx, col]
    if value isa Missing
        return missing
    end
    return Float64(value)
end

@inline function _huber(x::Float64, delta::Float64)::Float64
    ax = abs(x)
    if ax <= delta
        return 0.5 * x * x
    end
    return delta * (ax - 0.5 * delta)
end

function _objective_from_summary(
    summary_df::DataFrame,
    cfg::TunerConfig;
    run_failed::Bool,
    runtime_s::Float64,
    noise_rng::Union{Nothing, AbstractRNG}=nothing
)
    if run_failed || nrow(summary_df) == 0
        rt_pen = cfg.lambda_time * max(0.0, runtime_s / cfg.runtime_budget_s - 1.0)
        score = cfg.lambda_fail + rt_pen
        return (
            objective=score,
            base_loss=0.0,
            fail_penalty=cfg.lambda_fail,
            runtime_penalty=rt_pen,
            all_pass=false,
            failed_rows=0,
            per_scenario=""
        )
    end

    base_loss = 0.0
    all_pass = true
    failed_rows = 0
    scenario_acc = Dict{String, Float64}()

    for row in eachrow(summary_df)
        scenario = String(row.scenario)
        w = get(cfg.scenario_weights, scenario, 1.0)

        rmse = Float64(row.rmse_km)
        absv = Float64(row.max_abs_km)
        lim_rmse = Float64(row.limit_max_rmse_km)
        lim_abs = Float64(row.limit_max_abs_km)

        rmse_ratio = if isfinite(lim_rmse) && lim_rmse > 0.0
            rmse / lim_rmse
        elseif isfinite(Float64(row.limit_nmae)) && Float64(row.limit_nmae) > 0.0
            Float64(row.nmae) / Float64(row.limit_nmae)
        else
            rmse
        end
        abs_ratio = absv / max(lim_abs, 1.0e-12)

        if noise_rng !== nothing
            rmse_ratio = max(0.0, rmse_ratio * (1.0 + cfg.telemetry_noise_frac * randn(noise_rng)))
            abs_ratio = max(0.0, abs_ratio * (1.0 + cfg.telemetry_noise_frac * randn(noise_rng)))
        end

        term = w * (_huber(rmse_ratio, cfg.huber_delta) + 0.5 * _huber(abs_ratio, cfg.huber_delta))
        base_loss += term
        scenario_acc[scenario] = get(scenario_acc, scenario, 0.0) + term

        pass_row = Bool(row.pass)
        all_pass &= pass_row
        pass_row || (failed_rows += 1)
    end

    fail_pen = cfg.lambda_fail * (run_failed ? 1.0 : 0.0)
    runtime_pen = cfg.lambda_time * max(0.0, runtime_s / cfg.runtime_budget_s - 1.0)
    score = base_loss + fail_pen + runtime_pen

    parts = [string(k, "=", @sprintf("%.6f", scenario_acc[k])) for k in sort(collect(keys(scenario_acc)))]
    return (
        objective=score,
        base_loss=base_loss,
        fail_penalty=fail_pen,
        runtime_penalty=runtime_pen,
        all_pass=all_pass,
        failed_rows=failed_rows,
        per_scenario=join(parts, ";")
    )
end

function _parse_run_metrics(summary_df::DataFrame, wall_s::Float64)
    if nrow(summary_df) == 0
        return (runtime_s=wall_s, scenario_count=0, event_count=0)
    end

    sim_runtime = try
        Float64(summary_df.total_runtime_s[1])
    catch
        wall_s
    end
    scenarios = length(unique(String.(summary_df.scenario)))
    events = nrow(summary_df)
    return (runtime_s=sim_runtime, scenario_count=scenarios, event_count=events)
end

function evaluate_candidate(
    cfg::TunerConfig,
    stage::String,
    profile::Symbol,
    cand::TuneCandidate,
    base_manifest::Dict{String, Any};
    uncertainty_draw::Union{Nothing, Int}=nothing,
    progress_lock::Union{Nothing, ReentrantLock}=nothing
)
    paths = _candidate_paths(cfg.outdir, stage, cand.id; draw=uncertainty_draw)
    manifest_doc = _apply_candidate_to_manifest(base_manifest, cand, cfg; uncertainty_draw=uncertainty_draw)
    open(paths.manifest, "w") do io
        TOML.print(io, manifest_doc)
    end

    cmd = `$(Base.julia_cmd()) --project=$(joinpath(REPO_ROOT, ".AGORA")) --startup-file=no $(TELEMETRY_STUDY_SCRIPT) --profile=$(String(profile)) --manifest=$(paths.manifest) --enforce=$(cfg.enforce ? "1" : "0") --out-summary=$(paths.summary) --out-errors=$(paths.errors)`

    env_pairs = Pair{String, String}[
        "SPACEAGORA_TELEMETRY_SOLVER_MODE" => String(cand.values["solver_mode"]),
        "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => string(Float64(cand.values["dt_max_orbit_s"])),
        "SPACEAGORA_TELEMETRY_DT_MAX_ATM" => string(Float64(cand.values["dt_max_atm_s"]))
    ]
    append!(env_pairs, _candidate_runtime_policy_env_pairs())
    if profile == :quick && cfg.maxiters_quick !== nothing
        push!(env_pairs, "SPACEAGORA_TELEMETRY_SOLVER_MAXITERS" => string(cfg.maxiters_quick))
    elseif profile == :full && cfg.maxiters_full !== nothing
        push!(env_pairs, "SPACEAGORA_TELEMETRY_SOLVER_MAXITERS" => string(cfg.maxiters_full))
    end

    run_ok = true
    error_text = ""
    wall_s = @elapsed begin
        open(paths.log, "w") do io
            try
                withenv(env_pairs...) do
                    proc = run(pipeline(cmd, stdout=io, stderr=io); wait=false)
                    started = time()
                    next_heartbeat = started + TUNER_HEARTBEAT_SECONDS
                    while process_running(proc)
                        sleep(2.0)
                        now_s = time()
                        if now_s >= next_heartbeat
                            elapsed_min = (now_s - started) / 60.0
                            _progress_print(
                                progress_lock,
                                @sprintf(
                                    "[heartbeat %s/%s] candidate=%d elapsed=%.1f min solver=%s",
                                    stage,
                                    String(profile),
                                    cand.id,
                                    elapsed_min,
                                    String(cand.values["solver_mode"])
                                )
                            )
                            next_heartbeat += TUNER_HEARTBEAT_SECONDS
                        end
                    end
                    wait(proc)
                    if !Base.success(proc)
                        run_ok = false
                        error_text = "nonzero_exit_code=$(proc.exitcode)"
                    end
                end
            catch err
                run_ok = false
                error_text = sprint(showerror, err)
            end
        end
    end

    summary_df = run_ok && isfile(paths.summary) ? CSV.read(paths.summary, DataFrame) : DataFrame()
    errors_df = isfile(paths.errors) ? CSV.read(paths.errors, DataFrame) : DataFrame()

    runtime_info = _parse_run_metrics(summary_df, wall_s)
    obj_rng = uncertainty_draw === nothing ? nothing : MersenneTwister(hash((cfg.seed, "telemetry_noise", cand.id, uncertainty_draw)))
    obj = _objective_from_summary(
        summary_df,
        cfg;
        run_failed=(!run_ok || nrow(summary_df) == 0),
        runtime_s=runtime_info.runtime_s,
        noise_rng=obj_rng
    )

    row = merge(
        _candidate_row_nt(cand),
        (
            stage=stage,
            profile=String(profile),
            uncertainty_draw=uncertainty_draw === nothing ? missing : uncertainty_draw,
            success=run_ok && nrow(summary_df) > 0,
            pass=obj.all_pass,
            objective=obj.objective,
            objective_base=obj.base_loss,
            penalty_fail=obj.fail_penalty,
            penalty_time=obj.runtime_penalty,
            failed_rows=obj.failed_rows,
            scenario_objective_breakdown=obj.per_scenario,
            scenario_count=runtime_info.scenario_count,
            event_count=runtime_info.event_count,
            odyssey_peri_rmse_km=_event_metric(summary_df, "odyssey", "peri", :rmse_km),
            odyssey_apo_rmse_km=_event_metric(summary_df, "odyssey", "apo", :rmse_km),
            vex_peri_rmse_km=_event_metric(summary_df, "vex", "peri", :rmse_km),
            vex_apo_rmse_km=_event_metric(summary_df, "vex", "apo", :rmse_km),
            earth_peri_rmse_km=_event_metric(summary_df, "earth_gmat", "peri", :rmse_km),
            earth_apo_rmse_km=_event_metric(summary_df, "earth_gmat", "apo", :rmse_km),
            simulation_runtime_s=runtime_info.runtime_s,
            wall_runtime_s=wall_s,
            error_message=run_ok ? "" : (isempty(error_text) ? "run_failed_or_summary_missing" : error_text),
            summary_path=paths.summary,
            errors_path=paths.errors,
            manifest_path=paths.manifest,
            log_path=paths.log,
            errors_rows=nrow(errors_df)
        )
    )
    return row
end

