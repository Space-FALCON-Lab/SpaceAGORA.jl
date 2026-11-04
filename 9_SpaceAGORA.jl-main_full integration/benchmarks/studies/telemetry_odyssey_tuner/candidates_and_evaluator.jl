function _candidate_nt(c::TuneCandidate)
    return (
        candidate_id=c.id,
        epoch_shift_s=c.epoch_shift_s,
        ra_scale=c.ra_scale,
        rp_altitude_offset_m=c.rp_altitude_offset_m,
        i_offset_deg=c.i_offset_deg,
        aop_offset_deg=c.aop_offset_deg,
        raan_offset_deg=c.raan_offset_deg,
        ta_offset_deg=c.ta_offset_deg,
        bus_mass_scale=c.bus_mass_scale,
        prop_mass_scale=c.prop_mass_scale,
        panel_mass_scale=c.panel_mass_scale,
        bus_dims_scale=c.bus_dims_scale,
        panel_dims_scale=c.panel_dims_scale,
        panel_offset_scale=c.panel_offset_scale,
        srp_cr_scale=c.srp_cr_scale,
        srp_area_scale=c.srp_area_scale,
        gravity_variant=c.gravity_variant,
        dv_global_scale=c.dv_global_scale,
        dv_early_scale=c.dv_early_scale,
        dv_orbit7_bias_mps=c.dv_orbit7_bias_mps,
        solver_mode_requested=c.solver_mode,
        dt_max_orbit_requested_s=c.dt_max_orbit_s,
        dt_max_atm_requested_s=c.dt_max_atm_s
    )
end

@inline _pick(rng::AbstractRNG, values::Vector) = values[rand(rng, 1:length(values))]

function sample_candidates(cfg::TunerConfig)::Vector{TuneCandidate}
    rng = MersenneTwister(cfg.seed)
    candidates = TuneCandidate[]
    seen = Set{String}()
    while length(candidates) < cfg.quick_candidates
        cand = TuneCandidate(
            id=length(candidates) + 1,
            epoch_shift_s=_pick(rng, EPOCH_SHIFT_VALUES_S),
            ra_scale=_pick(rng, RA_SCALE_VALUES),
            rp_altitude_offset_m=_pick(rng, RP_ALTITUDE_OFFSET_VALUES_M),
            i_offset_deg=_pick(rng, I_OFFSET_VALUES_DEG),
            aop_offset_deg=_pick(rng, AOP_OFFSET_VALUES_DEG),
            raan_offset_deg=_pick(rng, RAAN_OFFSET_VALUES_DEG),
            ta_offset_deg=_pick(rng, TA_OFFSET_VALUES_DEG),
            bus_mass_scale=_pick(rng, BUS_MASS_SCALE_VALUES),
            prop_mass_scale=_pick(rng, PROP_MASS_SCALE_VALUES),
            panel_mass_scale=_pick(rng, PANEL_MASS_SCALE_VALUES),
            bus_dims_scale=_pick(rng, BUS_DIMS_SCALE_VALUES),
            panel_dims_scale=_pick(rng, PANEL_DIMS_SCALE_VALUES),
            panel_offset_scale=_pick(rng, PANEL_OFFSET_SCALE_VALUES),
            srp_cr_scale=_pick(rng, SRP_CR_SCALE_VALUES),
            srp_area_scale=_pick(rng, SRP_AREA_SCALE_VALUES),
            gravity_variant=_pick(rng, GRAVITY_VARIANTS),
            dv_global_scale=_pick(rng, DV_GLOBAL_SCALE_VALUES),
            dv_early_scale=_pick(rng, DV_EARLY_SCALE_VALUES),
            dv_orbit7_bias_mps=_pick(rng, DV_ORBIT7_BIAS_VALUES_MPS),
            solver_mode=_pick(rng, cfg.solver_modes),
            dt_max_orbit_s=_pick(rng, cfg.dt_max_orbit_values),
            dt_max_atm_s=_pick(rng, cfg.dt_max_atm_values)
        )
        sig = _candidate_signature(cand)
        if sig in seen
            continue
        end
        push!(seen, sig)
        push!(candidates, cand)
    end
    return candidates
end

function _find_base_scenario(base_manifest::Dict{String, Any}, name::String)::Dict{String, Any}
    haskey(base_manifest, "scenarios") || throw(ArgumentError("Manifest is missing top-level 'scenarios' array"))
    scenarios = base_manifest["scenarios"]
    scenarios isa AbstractVector || throw(ArgumentError("'scenarios' must be an array in manifest"))
    for sc in scenarios
        if sc isa AbstractDict && haskey(sc, "name") && String(sc["name"]) == name
            return deepcopy(Dict{String, Any}(sc))
        end
    end
    throw(ArgumentError("Scenario '$name' not found in manifest $(DEFAULT_BASE_MANIFEST)"))
end

function _apply_candidate(base::Dict{String, Any}, cand::TuneCandidate; calibration_enabled::Bool)::Dict{String, Any}
    sc = deepcopy(base)

    init = Dict{String, Any}(sc["initial_time"])
    sec_raw = Float64(init["second"])
    sec_int = floor(Int, sec_raw)
    sec_frac = sec_raw - sec_int
    sec_ms = round(Int, sec_frac * 1000)
    dt0 = DateTime(Int(init["year"]), Int(init["month"]), Int(init["day"]), Int(init["hour"]), Int(init["minute"]), sec_int, sec_ms)
    dt_shifted = dt0 + Second(cand.epoch_shift_s)
    init["year"] = year(dt_shifted)
    init["month"] = month(dt_shifted)
    init["day"] = day(dt_shifted)
    init["hour"] = hour(dt_shifted)
    init["minute"] = minute(dt_shifted)
    init["second"] = second(dt_shifted) + millisecond(dt_shifted) / 1000
    sc["initial_time"] = init

    sc["ra_m"] = _clamp_positive(Float64(sc["ra_m"]) * cand.ra_scale)
    sc["rp_altitude_m"] = _clamp_positive(Float64(sc["rp_altitude_m"]) + cand.rp_altitude_offset_m)
    sc["i_deg"] = clamp(Float64(sc["i_deg"]) + cand.i_offset_deg, 0.0, 180.0)
    sc["aop_deg"] = _wrap_0_360(Float64(sc["aop_deg"]) + cand.aop_offset_deg)
    sc["raan_deg"] = _wrap_0_360(Float64(sc["raan_deg"]) + cand.raan_offset_deg)
    sc["ta_deg"] = _wrap_0_360(Float64(sc["ta_deg"]) + cand.ta_offset_deg)

    if haskey(sc, "spacecraft")
        spacecraft = Dict{String, Any}(sc["spacecraft"])
        spacecraft["bus_mass_kg"] = _clamp_positive(Float64(spacecraft["bus_mass_kg"]) * cand.bus_mass_scale)
        spacecraft["prop_mass_kg"] = _clamp_positive(Float64(spacecraft["prop_mass_kg"]) * cand.prop_mass_scale)
        spacecraft["panel_mass_each_kg"] = _clamp_positive(Float64(spacecraft["panel_mass_each_kg"]) * cand.panel_mass_scale)
        spacecraft["panel_offset_y_m"] = _clamp_positive(Float64(spacecraft["panel_offset_y_m"]) * cand.panel_offset_scale)
        spacecraft["bus_dims_m"] = [Float64(v) * cand.bus_dims_scale for v in spacecraft["bus_dims_m"]]
        spacecraft["panel_dims_m"] = [Float64(v) * cand.panel_dims_scale for v in spacecraft["panel_dims_m"]]
        sc["spacecraft"] = spacecraft
    end

    if haskey(sc, "srp_cr")
        sc["srp_cr"] = _clamp_positive(Float64(sc["srp_cr"]) * cand.srp_cr_scale)
    end
    if haskey(sc, "srp_area_m2")
        sc["srp_area_m2"] = _clamp_positive(Float64(sc["srp_area_m2"]) * cand.srp_area_scale)
    end

    if cand.gravity_variant == "j2_only"
        sc["gravity_model"] = "inverse_squared_j2"
        sc["gravity_harmonics_degree"] = 0
        sc["gravity_harmonics_order"] = 0
    else
        degree = if cand.gravity_variant == "harm_deg2"
            2
        elseif cand.gravity_variant == "harm_deg4"
            4
        elseif cand.gravity_variant == "harm_deg8"
            8
        else
            20
        end
        sc["gravity_model"] = "inverse_squared"
        sc["gravity_harmonics_degree"] = degree
        sc["gravity_harmonics_order"] = degree
    end

    if haskey(sc, "maneuvers")
        maneuvers = Dict{String, Any}(sc["maneuvers"])
        orbits = [Int(v) for v in maneuvers["orbit_numbers"]]
        delta_v = [Float64(v) for v in maneuvers["delta_v_mps"]]
        tuned = similar(delta_v)
        for i in eachindex(delta_v)
            scale = cand.dv_global_scale
            if orbits[i] <= 50
                scale *= cand.dv_early_scale
            end
            dv = delta_v[i] * scale
            if orbits[i] == 7
                dv += cand.dv_orbit7_bias_mps
            end
            tuned[i] = dv
        end
        maneuvers["delta_v_mps"] = tuned
        sc["maneuvers"] = maneuvers
    end

    if haskey(sc, "calibration")
        cal = Dict{String, Any}(sc["calibration"])
        cal["enabled"] = calibration_enabled
        sc["calibration"] = cal
    end
    return sc
end

@inline function _paper_ratio(event::String, rmse_km::Float64, max_abs_km::Float64)
    lim = PAPER_LIMITS[event]
    return rmse_km / lim.rmse_km, max_abs_km / lim.max_abs_km
end

@inline function _score_candidate(
    peri_rmse::Float64,
    apo_rmse::Float64,
    peri_abs::Float64,
    apo_abs::Float64
)::Float64
    peri_rmse_ratio, peri_abs_ratio = _paper_ratio("peri", peri_rmse, peri_abs)
    apo_rmse_ratio, apo_abs_ratio = _paper_ratio("apo", apo_rmse, apo_abs)
    ratios = [peri_rmse_ratio, peri_abs_ratio, apo_rmse_ratio, apo_abs_ratio]
    return mean(ratios) + 0.35 * maximum(ratios)
end

@inline function _passes_paper_limits(
    peri_rmse::Float64,
    apo_rmse::Float64,
    peri_abs::Float64,
    apo_abs::Float64
)::Bool
    peri_ok = peri_rmse <= PAPER_LIMITS["peri"].rmse_km && peri_abs <= PAPER_LIMITS["peri"].max_abs_km
    apo_ok = apo_rmse <= PAPER_LIMITS["apo"].rmse_km && apo_abs <= PAPER_LIMITS["apo"].max_abs_km
    return peri_ok && apo_ok
end

function _event_row(df::DataFrame, event::String)
    idx = findfirst((df.event .== event))
    idx === nothing && throw(ArgumentError("Missing event '$event' in summary dataframe"))
    return df[idx, :]
end

function _candidate_paths(outdir::String, stage::Symbol, id::Int)
    tag = @sprintf("cand_%04d_%s", id, String(stage))
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

function evaluate_candidate(
    cfg::TunerConfig,
    stage::Symbol,
    cand::TuneCandidate,
    tuned_scenario::Dict{String, Any};
    progress_lock::Union{Nothing, ReentrantLock}=nothing
)
    paths = _candidate_paths(cfg.outdir, stage, cand.id)
    manifest_doc = Dict{String, Any}("version" => 1, "scenarios" => Any[tuned_scenario])
    open(paths.manifest, "w") do io
        TOML.print(io, manifest_doc)
    end

    profile = String(stage)
    cmd = `$(Base.julia_cmd()) --project=$(joinpath(REPO_ROOT, ".AGORA")) --startup-file=no $(TELEMETRY_STUDY_SCRIPT) --profile=$(profile) --manifest=$(paths.manifest) --enforce=$(cfg.enforce ? "1" : "0") --out-summary=$(paths.summary) --out-errors=$(paths.errors)`

    env_pairs = Pair{String, String}[
        "SPACEAGORA_TELEMETRY_SOLVER_MODE" => cand.solver_mode,
        "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => string(cand.dt_max_orbit_s),
        "SPACEAGORA_TELEMETRY_DT_MAX_ATM" => string(cand.dt_max_atm_s)
    ]
    append!(env_pairs, _candidate_runtime_policy_env_pairs())
    if stage == :quick && cfg.maxiters_quick !== nothing
        push!(env_pairs, "SPACEAGORA_TELEMETRY_SOLVER_MAXITERS" => string(cfg.maxiters_quick))
    elseif stage == :full && cfg.maxiters_full !== nothing
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
                                    "[heartbeat %s] candidate=%d elapsed=%.1f min solver=%s dt_orbit=%.1f dt_atm=%.3f",
                                    profile,
                                    cand.id,
                                    elapsed_min,
                                    cand.solver_mode,
                                    cand.dt_max_orbit_s,
                                    cand.dt_max_atm_s
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

    if !run_ok || !isfile(paths.summary)
        return merge(
            _candidate_nt(cand),
            (
                stage=profile,
                success=false,
                pass=false,
                score=Inf,
                error_message=isempty(error_text) ? "run_failed_or_summary_missing" : error_text,
                peri_rmse_km=missing,
                apo_rmse_km=missing,
                peri_max_abs_km=missing,
                apo_max_abs_km=missing,
                peri_nmae=missing,
                apo_nmae=missing,
                solver_mode_reported=missing,
                solver_sequence=missing,
                solver_fallback_count=missing,
                solver_retcode=missing,
                coverage_peri=missing,
                coverage_apo=missing,
                simulation_runtime_s=missing,
                wall_runtime_s=wall_s,
                summary_path=paths.summary,
                errors_path=paths.errors,
                manifest_path=paths.manifest,
                log_path=paths.log
            )
        )
    end

    summary_df = CSV.read(paths.summary, DataFrame)
    scenario_df = summary_df[summary_df.scenario .== cfg.scenario_name, :]
    if nrow(scenario_df) == 0
        return merge(
            _candidate_nt(cand),
            (
                stage=profile,
                success=false,
                pass=false,
                score=Inf,
                error_message="summary_missing_scenario",
                peri_rmse_km=missing,
                apo_rmse_km=missing,
                peri_max_abs_km=missing,
                apo_max_abs_km=missing,
                peri_nmae=missing,
                apo_nmae=missing,
                solver_mode_reported=missing,
                solver_sequence=missing,
                solver_fallback_count=missing,
                solver_retcode=missing,
                coverage_peri=missing,
                coverage_apo=missing,
                simulation_runtime_s=missing,
                wall_runtime_s=wall_s,
                summary_path=paths.summary,
                errors_path=paths.errors,
                manifest_path=paths.manifest,
                log_path=paths.log
            )
        )
    end

    peri = _event_row(scenario_df, "peri")
    apo = _event_row(scenario_df, "apo")
    peri_rmse = Float64(peri.rmse_km)
    apo_rmse = Float64(apo.rmse_km)
    peri_abs = Float64(peri.max_abs_km)
    apo_abs = Float64(apo.max_abs_km)
    score = _score_candidate(peri_rmse, apo_rmse, peri_abs, apo_abs)
    pass = _passes_paper_limits(peri_rmse, apo_rmse, peri_abs, apo_abs)

    return merge(
        _candidate_nt(cand),
        (
            stage=profile,
            success=true,
            pass=pass,
            score=score,
            error_message="",
            peri_rmse_km=peri_rmse,
            apo_rmse_km=apo_rmse,
            peri_max_abs_km=peri_abs,
            apo_max_abs_km=apo_abs,
            peri_nmae=Float64(peri.nmae),
            apo_nmae=Float64(apo.nmae),
            solver_mode_reported=String(peri.solver_mode),
            solver_sequence=String(peri.solver_sequence),
            solver_fallback_count=Int(peri.solver_fallback_count),
            solver_retcode=String(peri.solver_retcode),
            coverage_peri=Float64(peri.coverage),
            coverage_apo=Float64(apo.coverage),
            simulation_runtime_s=Float64(peri.simulation_runtime_s),
            wall_runtime_s=wall_s,
            summary_path=paths.summary,
            errors_path=paths.errors,
            manifest_path=paths.manifest,
            log_path=paths.log
        )
    )
end

