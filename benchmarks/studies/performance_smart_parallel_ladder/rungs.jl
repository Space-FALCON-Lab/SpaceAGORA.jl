function _ladder_rungs(config::SmartLadderConfig)::Vector{LadderRungSpec}
    rungs = LadderRungSpec[
        LadderRungSpec(
            mode=:serial,
            label="r0_true_serial",
            description="No outer parallelism, no inner callback parallelism, no adaptive policies.",
            matrix=:outer_pinned,
            backend="none",
            inner_adaptive=false,
            outer_route_adaptive=false
        ),
        LadderRungSpec(
            mode=:outer_only_process,
            label="r1_b_outer_only_process",
            description="Outer-loop process parallel only; inner callback parallelism pinned off.",
            matrix=:outer_pinned,
            backend="process",
            inner_adaptive=false,
            outer_route_adaptive=false
        ),
        LadderRungSpec(
            mode=:outer_only,
            label="r1_a_outer_only",
            description="Outer-loop parallel only; inner callback parallelism pinned off.",
            matrix=:outer_pinned,
            backend=config.outer_only_backend,
            inner_adaptive=false,
            outer_route_adaptive=false
        ),
        LadderRungSpec(
            mode=:inner_only,
            label="r2_inner_only",
            description="Inner callback auto policies only; outer loop serial.",
            matrix=:full_auto,
            backend="none",
            inner_adaptive=false,
            outer_route_adaptive=false
        ),
        LadderRungSpec(
            mode=:outer_inner_static,
            label="r3_outer_inner_static",
            description="Outer + inner enabled with static policy behavior (adaptive off).",
            matrix=:full_auto,
            backend="auto",
            inner_adaptive=false,
            outer_route_adaptive=false
        ),
        LadderRungSpec(
            mode=:outer_inner_adaptive,
            label="r4_outer_inner_adaptive",
            description="Outer + inner enabled with adaptive policy behavior (baseline full-auto heuristic).",
            matrix=:full_auto,
            backend="auto",
            inner_adaptive=true,
            outer_route_adaptive=true
        ),
        LadderRungSpec(
            mode=:outer_inner_full_smart,
            label="r5_full_smart",
            description="Outer + inner enabled with tuned adaptive policy behavior (thermal forced on, shorter adaptive window, control tail-guard).",
            matrix=:full_auto,
            backend="auto",
            inner_adaptive=true,
            outer_route_adaptive=true
        )
    ]
    if _smart_ladder_smoke_mode()
        return LadderRungSpec[
            rungs[1],
            rungs[3],
            rungs[7],
        ]
    end
    if !config.include_layer_attribution
        return rungs
    end

    attribution_backend = _layer_attribution_backend(config)
    append!(
        rungs,
        LadderRungSpec[
            LadderRungSpec(
                mode=:attr_outer_only,
                label="la_outer_only",
                description="Layer attribution baseline: outer parallel enabled, all inner layers off.",
                matrix=:outer_pinned,
                backend=attribution_backend,
                inner_adaptive=false,
                outer_route_adaptive=false
            ),
            LadderRungSpec(
                mode=:attr_density,
                label="la_density",
                description="Outer + density inner layer only.",
                matrix=:attribution_density,
                backend=attribution_backend,
                inner_adaptive=false,
                outer_route_adaptive=false
            ),
            LadderRungSpec(
                mode=:attr_thermal,
                label="la_thermal",
                description="Outer + thermal inner layer only.",
                matrix=:attribution_thermal,
                backend=attribution_backend,
                inner_adaptive=false,
                outer_route_adaptive=false
            ),
            LadderRungSpec(
                mode=:attr_control,
                label="la_control",
                description="Outer + control inner layer only.",
                matrix=:attribution_control,
                backend=attribution_backend,
                inner_adaptive=false,
                outer_route_adaptive=false
            ),
            LadderRungSpec(
                mode=:attr_multibody,
                label="la_multibody",
                description="Outer + multibody inner layer only.",
                matrix=:attribution_multibody,
                backend=attribution_backend,
                inner_adaptive=false,
                outer_route_adaptive=false
            ),
            LadderRungSpec(
                mode=:attr_effector,
                label="la_effector",
                description="Outer + effector inner layer only.",
                matrix=:attribution_effector,
                backend=attribution_backend,
                inner_adaptive=false,
                outer_route_adaptive=false
            ),
            LadderRungSpec(
                mode=:attr_density_thermal,
                label="la_density_thermal",
                description="Outer + density + thermal pair interaction.",
                matrix=:attribution_density_thermal,
                backend=attribution_backend,
                inner_adaptive=false,
                outer_route_adaptive=false
            ),
            LadderRungSpec(
                mode=:attr_density_multibody,
                label="la_density_multibody",
                description="Outer + density + multibody pair interaction.",
                matrix=:attribution_density_multibody,
                backend=attribution_backend,
                inner_adaptive=false,
                outer_route_adaptive=false
            ),
            LadderRungSpec(
                mode=:attr_control_effector,
                label="la_control_effector",
                description="Outer + control + effector pair interaction.",
                matrix=:attribution_control_effector,
                backend=attribution_backend,
                inner_adaptive=false,
                outer_route_adaptive=false
            ),
            LadderRungSpec(
                mode=:attr_multibody_effector,
                label="la_multibody_effector",
                description="Outer + multibody + effector pair interaction.",
                matrix=:attribution_multibody_effector,
                backend=attribution_backend,
                inner_adaptive=false,
                outer_route_adaptive=false
            ),
        ]
    )
    return rungs
end

function _ladder_env_pairs(
    rung::LadderRungSpec,
    config::SmartLadderConfig;
    solver_mode::Union{Nothing, String}=nothing
)::Vector{Pair{String, Union{Nothing, String}}}
    matrix = _ladder_matrix_modes(rung.matrix)
    r5_full_smart = rung.mode == :outer_inner_full_smart
    serial_rung = rung.mode == :serial
    thermal_mode = r5_full_smart ? "on" : matrix.thermal
    inner_scheduler = r5_full_smart ? "dynamic" : "static"
    measured_reward = r5_full_smart ? "1" : "0"
    persistent_hints = r5_full_smart ? "1" : "0"
    persistent_state_persist = r5_full_smart ? "1" : "0"
    if _smart_ladder_smoke_mode()
        persistent_hints = "0"
        persistent_state_persist = "0"
    end
    pairs = Pair{String, Union{Nothing, String}}[
        "SPACEAGORA_PARALLEL_PROFILE" => (serial_rung ? "R0" : nothing),
        "SPACEAGORA_PERF_PARALLEL_BACKEND" => rung.backend,
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => (rung.inner_adaptive ? "1" : "0"),
        "SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE" => (rung.outer_route_adaptive ? "1" : "0"),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => matrix.density,
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => thermal_mode,
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => matrix.control,
        "SPACEAGORA_MULTIBODY_PARALLEL" => matrix.multibody,
        "SPACEAGORA_EFFECTOR_PARALLEL" => matrix.effector,
        "SPACEAGORA_RHS_BATCH_PARALLEL" => (serial_rung ? "off" : "auto"),
        "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => inner_scheduler,
        "SPACEAGORA_PARALLEL_POLICY_WINDOW" => (rung.mode == :outer_inner_full_smart ? "4" : "8"),
        "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD" => (rung.mode == :outer_inner_full_smart ? "1" : "0"),
        "SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD" => measured_reward,
        "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => persistent_hints,
        "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => persistent_state_persist,
        "SPACEAGORA_PERF_INCLUDE_CONTROL_STRESS_PER_ORBIT" => (config.include_control_stress_per_orbit ? "1" : "0"),
        "SPACEAGORA_PERF_CONTROL_STRESS_REPEATS_FULL" => string(config.control_stress_repeats_full),
        "SPACEAGORA_PERF_CONTROL_STRESS_WARMUP_FULL" => string(config.control_stress_warmup_full)
    ]
    if isnothing(config.process_workers)
        push!(pairs, "SPACEAGORA_PERF_PROCS" => nothing)
    else
        push!(pairs, "SPACEAGORA_PERF_PROCS" => string(config.process_workers))
    end
    if !(solver_mode === nothing)
        push!(pairs, "SPACEAGORA_PERF_SOLVER_MODE" => String(solver_mode))
    end
    return pairs
end

function _solver_variants(config::SmartLadderConfig)::Vector{NamedTuple{(:label, :solver_mode), Tuple{String, Union{Nothing, String}}}}
    if config.solver_axis == :inherit
        return [(label="inherit", solver_mode=nothing)]
    elseif config.solver_axis == :frozen
        mode = config.solver_mode
        return [(label=_safe_path_token(mode), solver_mode=mode)]
    elseif config.solver_axis == :factorial
        variants = NamedTuple{(:label, :solver_mode), Tuple{String, Union{Nothing, String}}}[]
        for mode in config.solver_factor_modes
            push!(variants, (label=_safe_path_token(mode), solver_mode=mode))
        end
        return variants
    end
    error("Unsupported solver-axis '$(config.solver_axis)'.")
end

@inline function _primary_solver_variant(config::SmartLadderConfig)::NamedTuple{(:label, :solver_mode), Tuple{String, Union{Nothing, String}}}
    variants = _solver_variants(config)
    isempty(variants) && error("No solver variants resolved for solver-axis '$(config.solver_axis)'.")
    return first(variants)
end

function _latest_artifact_path(outdir::String, prefix::String, profile_name::String, suffix::String)::String
    pattern = "$(prefix)_$(profile_name)_"
    candidates = String[]
    for entry in readdir(outdir)
        if startswith(entry, pattern) && endswith(entry, suffix)
            push!(candidates, joinpath(outdir, entry))
        end
    end
    isempty(candidates) && error("No artifact found in '$outdir' for prefix='$prefix' profile='$(profile_name)' suffix='$suffix'")
    sort!(candidates; by=path -> stat(path).mtime)
    return last(candidates)
end

function _latest_artifact_path_optional(outdir::String, prefix::String, profile_name::String, suffix::String)::Union{Nothing, String}
    pattern = "$(prefix)_$(profile_name)_"
    candidates = String[]
    for entry in readdir(outdir)
        if startswith(entry, pattern) && endswith(entry, suffix)
            push!(candidates, joinpath(outdir, entry))
        end
    end
    isempty(candidates) && return nothing
    sort!(candidates; by=path -> stat(path).mtime)
    return last(candidates)
end

@inline function _git_head_commit(repo_root::String)::String
    try
        return strip(read(`git -C $(repo_root) rev-parse HEAD`, String))
    catch
        return "unknown"
    end
end

@inline function _manifest_path_for_project(project_path::String)::Union{Nothing, String}
    primary = joinpath(project_path, "Manifest.toml")
    if isfile(primary)
        return primary
    end
    fallback = joinpath(REPO_ROOT, "Manifest.toml")
    if isfile(fallback)
        return fallback
    end
    return nothing
end

@inline function _sha256_file(path::String)::String
    return bytes2hex(sha256(read(path)))
end

@inline function _cmd_string(cmd::Cmd)::String
    return sprint(io -> print(io, cmd))
end

function _write_rung_repro_metadata(
    rung_outdir::String,
    rung::LadderRungSpec,
    config::SmartLadderConfig,
    pass_idx::Int,
    solver_label::String,
    solver_mode::Union{Nothing, String},
    cmd::Cmd,
    env_pairs::Vector{Pair{String, Union{Nothing, String}}}
)::NamedTuple
    manifest_path = _manifest_path_for_project(SMART_LADDER_PROJECT)
    manifest_sha256 = isnothing(manifest_path) ? missing : _sha256_file(manifest_path)
    metadata_path = joinpath(
        rung_outdir,
        "smart_parallel_ladder_rung_metadata_$(config.profile.name)_$(rung.label).csv"
    )
    env_path = joinpath(
        rung_outdir,
        "smart_parallel_ladder_rung_env_$(config.profile.name)_$(rung.label).csv"
    )

    metadata_df = DataFrame([
        (
            timestamp_utc=string(now(UTC)),
            profile=config.profile.name,
            pass=pass_idx,
            solver_label=solver_label,
            solver_mode=(solver_mode === nothing ? missing : solver_mode),
            rung=rung.label,
            mode=String(rung.mode),
            matrix=String(rung.matrix),
            backend=rung.backend,
            inner_adaptive=rung.inner_adaptive,
            outer_route_adaptive=rung.outer_route_adaptive,
            project_path=SMART_LADDER_PROJECT,
            command=_cmd_string(cmd),
            julia_version=string(VERSION),
            julia_threads=Threads.nthreads(),
            cpu_threads=Sys.CPU_THREADS,
            git_commit=_git_head_commit(REPO_ROOT),
            manifest_path=isnothing(manifest_path) ? missing : manifest_path,
            manifest_sha256=manifest_sha256
        )
    ])
    CSV.write(metadata_path, metadata_df)

    env_rows = NamedTuple[]
    for (key, value) in env_pairs
        push!(env_rows, (
            key=String(key),
            value=(value === nothing ? missing : String(value))
        ))
    end
    env_df = DataFrame(env_rows)
    if nrow(env_df) > 0
        sort!(env_df, :key)
    end
    CSV.write(env_path, env_df)
    return (metadata_path=metadata_path, env_path=env_path)
end

@inline function _has_column(df::DataFrame, column::Symbol)::Bool
    return column in propertynames(df)
end

function _stage_elapsed_s(
    stage_timing_path::Union{Nothing, String},
    fallback_total_s::Float64
)::Tuple{Float64, Float64, Float64, Float64, Float64}
    if isnothing(stage_timing_path)
        return fallback_total_s, 0.0, 0.0, 0.0, fallback_total_s
    end

    stage_df = CSV.read(stage_timing_path, DataFrame)
    bench_elapsed_s = 0.0
    split_gate_elapsed_s = 0.0
    orbit_elapsed_s = 0.0
    entry_duration_elapsed_s = 0.0
    total_elapsed_s = fallback_total_s

    if _has_column(stage_df, :stage) && _has_column(stage_df, :elapsed_s)
        for row in eachrow(stage_df)
            stage_name = String(row.stage)
            elapsed = Float64(row.elapsed_s)
            if stage_name == "run_benchmarks"
                bench_elapsed_s = elapsed
            elseif stage_name == "run_split_rollout_gate"
                split_gate_elapsed_s = elapsed
            elseif stage_name == "run_per_orbit"
                orbit_elapsed_s = elapsed
            elseif stage_name == "run_entry_duration_sweep"
                entry_duration_elapsed_s = elapsed
            elseif stage_name == "total"
                total_elapsed_s = elapsed
            end
        end
    end

    if total_elapsed_s <= 0.0
        total_elapsed_s = bench_elapsed_s + split_gate_elapsed_s + orbit_elapsed_s + entry_duration_elapsed_s
    end
    if total_elapsed_s <= 0.0
        total_elapsed_s = fallback_total_s
    end
    if bench_elapsed_s <= 0.0 && split_gate_elapsed_s <= 0.0 && orbit_elapsed_s <= 0.0 && entry_duration_elapsed_s <= 0.0
        bench_elapsed_s = total_elapsed_s
    end

    return bench_elapsed_s, split_gate_elapsed_s, orbit_elapsed_s, entry_duration_elapsed_s, total_elapsed_s
end

function _ladder_artifacts(
    rung::LadderRungSpec,
    config::SmartLadderConfig,
    rung_outdir::String,
    elapsed_s::Float64
)::ModeRunArtifacts
    profile_name = config.profile.name
    raw_path = _latest_artifact_path(rung_outdir, "runtime_raw", profile_name, ".csv")
    summary_path = _latest_artifact_path(rung_outdir, "runtime_summary", profile_name, ".csv")
    orbit_raw_path = _latest_artifact_path(rung_outdir, "runtime_per_orbit_raw", profile_name, ".csv")
    orbit_summary_path = _latest_artifact_path(rung_outdir, "runtime_per_orbit_summary", profile_name, ".csv")
    entry_duration_raw_path = _latest_artifact_path(rung_outdir, "runtime_entry_duration_raw", profile_name, ".csv")
    entry_duration_summary_path = _latest_artifact_path(rung_outdir, "runtime_entry_duration_summary", profile_name, ".csv")
    stage_timing_path = _latest_artifact_path_optional(rung_outdir, "runtime_stage_timing", profile_name, ".csv")
    hardware_info_path = _latest_artifact_path_optional(rung_outdir, "runtime_hardware_info", profile_name, ".csv")
    split_gate_csv_path = _latest_artifact_path_optional(rung_outdir, "split_rollout_gate", profile_name, ".csv")
    split_gate_report_path = _latest_artifact_path_optional(rung_outdir, "split_rollout_gate", profile_name, ".md")
    report_path = _latest_artifact_path(rung_outdir, "runtime_report", profile_name, ".md")

    bench_elapsed_s, split_gate_elapsed_s, orbit_elapsed_s, entry_duration_elapsed_s, total_elapsed_s = _stage_elapsed_s(stage_timing_path, elapsed_s)
    split_gate_df = isnothing(split_gate_csv_path) ? nothing : CSV.read(split_gate_csv_path, DataFrame)
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
        stage_timing_path=isnothing(stage_timing_path) ? "" : stage_timing_path,
        hardware_info_path=isnothing(hardware_info_path) ? "" : hardware_info_path,
        split_gate_elapsed_s=split_gate_elapsed_s,
        split_gate_csv_path=split_gate_csv_path,
        split_gate_report_path=split_gate_report_path,
        split_gate_df=split_gate_df,
        raw_df=CSV.read(raw_path, DataFrame),
        summary_df=CSV.read(summary_path, DataFrame),
        orbit_summary_df=CSV.read(orbit_summary_path, DataFrame)
    )
end
