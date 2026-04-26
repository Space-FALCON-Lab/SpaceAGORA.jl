function run_rung(
    rung::LadderRungSpec,
    config::SmartLadderConfig,
    pass_idx::Int,
    pass_outdir::String;
    solver_label::String="inherit",
    solver_mode::Union{Nothing, String}=nothing
)::ModeRunArtifacts
    rung_outdir = joinpath(pass_outdir, rung.label)
    mkpath(rung_outdir)

    println()
    println(
        "[smart-ladder] pass=$(pass_idx) rung=$(rung.label) backend=$(rung.backend) matrix=$(rung.matrix) " *
        "inner_adaptive=$(rung.inner_adaptive) outer_route_adaptive=$(rung.outer_route_adaptive)"
    )
    cmd = `$(Base.julia_cmd()) --project=$(SMART_LADDER_PROJECT) $(SMART_LADDER_RUNTIME_SCRIPT) --profile=$(config.profile.name) --outdir=$(rung_outdir)`
    env_pairs = _ladder_env_pairs(rung, config; solver_mode=solver_mode)
    metadata_paths = _write_rung_repro_metadata(
        rung_outdir,
        rung,
        config,
        pass_idx,
        solver_label,
        solver_mode,
        cmd,
        env_pairs
    )

    started_ns = time_ns()
    withenv(env_pairs...) do
        run(cmd)
    end
    elapsed_s = (time_ns() - started_ns) / 1e9

    artifacts = _ladder_artifacts(rung, config, rung_outdir, elapsed_s)
    println(
        "[smart-ladder] pass=$(pass_idx) rung=$(rung.label) completed total=$(round(artifacts.elapsed_s; digits=3)) s " *
        "(run_benchmarks=$(round(artifacts.bench_elapsed_s; digits=3)) s, split_gate=$(round(artifacts.split_gate_elapsed_s; digits=3)) s, " *
        "per_orbit=$(round(artifacts.orbit_elapsed_s; digits=3)) s, solver=$(solver_label), metadata=$(metadata_paths.metadata_path), env=$(metadata_paths.env_path))"
    )
    return artifacts
end

function _tag_rung_column(
    df::DataFrame,
    rung::LadderRungSpec;
    pass_idx::Union{Nothing, Int}=nothing,
    solver_label::Union{Nothing, String}=nothing,
    solver_mode::Union{Nothing, String}=nothing
)::DataFrame
    tagged = copy(df)
    tagged[!, :rung] = fill(rung.label, nrow(tagged))
    tagged[!, :rung_mode] = fill(String(rung.mode), nrow(tagged))
    tagged[!, :rung_matrix] = fill(String(rung.matrix), nrow(tagged))
    tagged[!, :rung_backend] = fill(rung.backend, nrow(tagged))
    tagged[!, :rung_inner_adaptive] = fill(rung.inner_adaptive, nrow(tagged))
    tagged[!, :rung_outer_route_adaptive] = fill(rung.outer_route_adaptive, nrow(tagged))
    if !(solver_label === nothing)
        tagged[!, :solver_label] = fill(solver_label, nrow(tagged))
    end
    if !(solver_mode === nothing)
        tagged[!, :solver_mode_factor] = fill(solver_mode, nrow(tagged))
    end
    if !isnothing(pass_idx)
        tagged[!, :pass] = fill(pass_idx, nrow(tagged))
    end
    return tagged
end

function _ordered_rungs(
    rungs::Vector{LadderRungSpec},
    rng::Random.AbstractRNG,
    randomize::Bool
)::Vector{LadderRungSpec}
    if !randomize || length(rungs) <= 1
        return copy(rungs)
    end
    perm = randperm(rng, length(rungs))
    return rungs[perm]
end

@inline function _safe_ratio(num::Real, den::Real)::Union{Missing, Float64}
    n = Float64(num)
    d = Float64(den)
    if !isfinite(n) || !isfinite(d) || d <= 0.0
        return missing
    end
    return n / d
end

@inline function _safe_ratio(num::Missing, den)::Union{Missing, Float64}
    return missing
end

@inline function _safe_ratio(num, den::Missing)::Union{Missing, Float64}
    return missing
end

@inline function _safe_ratio(num, den)::Union{Missing, Float64}
    if !(num isa Real) || !(den isa Real)
        return missing
    end
    return _safe_ratio(num::Real, den::Real)
end

@inline function _safe_cv_pct(values::Vector{Float64})::Union{Missing, Float64}
    n = length(values)
    if n == 0
        return missing
    elseif n == 1
        return 0.0
    end
    μ = mean(values)
    abs(μ) <= eps(Float64) && return missing
    return 100.0 * std(values; corrected=true) / abs(μ)
end

@inline function _safe_ci95_bounds(values::Vector{Float64})::Tuple{Union{Missing, Float64}, Union{Missing, Float64}}
    n = length(values)
    if n == 0
        return (missing, missing)
    elseif n == 1
        μ = values[1]
        return (μ, μ)
    end
    μ = mean(values)
    σ = std(values; corrected=true)
    sem = σ / sqrt(n)
    margin = 1.96 * sem
    return (μ - margin, μ + margin)
end

@inline function _pass_stats(values::Vector{Float64})::NamedTuple
    n = length(values)
    if n == 0
        return (
            samples=0,
            mean_s=missing,
            std_s=missing,
            cv_pct=missing,
            ci95_low_s=missing,
            ci95_high_s=missing
        )
    end
    ci_low, ci_high = _safe_ci95_bounds(values)
    return (
        samples=n,
        mean_s=mean(values),
        std_s=(n == 1 ? 0.0 : std(values; corrected=true)),
        cv_pct=_safe_cv_pct(values),
        ci95_low_s=ci_low,
        ci95_high_s=ci_high
    )
end

@inline function _parse_bool_env_default(name::String, default::Bool)::Bool
    raw = get(ENV, name, default ? "1" : "0")
    try
        return _parse_bool_token(raw)
    catch
        return default
    end
end

@inline function _protocol_mode_tags(config::SmartLadderConfig)::NamedTuple
    cache_cold = config.clean
    outer_state_reset = _parse_bool_env_default("SPACEAGORA_PERF_OUTER_ROUTE_STATE_RESET", false)
    inner_state_reset = _parse_bool_env_default("SPACEAGORA_PARALLEL_POLICY_STATE_RESET", false)
    start_mode =
        (cache_cold && outer_state_reset && inner_state_reset) ? "cold" :
        (!cache_cold && !outer_state_reset && !inner_state_reset) ? "warm" :
        "mixed"
    return (
        start_mode=start_mode,
        cache_mode=(cache_cold ? "cold_cache" : "warm_cache"),
        cache_cold=cache_cold,
        outer_state_reset=outer_state_reset,
        inner_state_reset=inner_state_reset
    )
end

