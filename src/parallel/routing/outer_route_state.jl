Base.@kwdef struct OuterRouteFeatures
    category::String = "deterministic"
    n_sats::Int = 1
    n_links::Int = 1
    max_links_per_sat::Int = 1
    mission_time_s::Float64 = 0.0
    has_nbody::Bool = false
    has_srp::Bool = false
    harmonics_degree::Int = 0
    has_control::Bool = false
    orientation_on::Bool = false
    density_family::String = "unknown"
    solver_mode::String = "auto"
    dt_max_orbit_s::Float64 = 0.0
    control_rate_s::Float64 = 0.0
    guidance_rate_s::Float64 = 0.0
    navigation_rate_s::Float64 = 0.0
    gram_surrogate_enabled::Bool = false
    gram_static_grid_enabled::Bool = false
    control_effector_count::Int = 0
    thermal_enabled::Bool = false
    dynamic_effector_count::Int = 0
    effector_cost_class::String = "unknown"
    montecarlo_samples::Int = 0
end

Base.@kwdef struct OuterRouteTuning
    inner_sat_threshold::Int = 8
    inner_link_threshold::Int = 12
    outer_light_sat_threshold::Int = 2
    outer_light_link_threshold::Int = 4
    outer_light_mission_threshold_s::Float64 = 14_400.0
    spice_constellation_process_enabled::Bool = true
    spice_constellation_min_sats::Int = 4
    adaptive_enabled::Bool = true
    adaptive_min_samples::Int = 2
    adaptive_exploration_c::Float64 = 1.25
    failure_penalty_s::Float64 = 120.0
    mc_process_min_samples::Int = 16
    mc_process_min_mission_s::Float64 = 3600.0
    trace::Bool = false
end

Base.@kwdef mutable struct OuterRouteStats
    samples::Int = 0
    successes::Int = 0
    failures::Int = 0
    elapsed_sum_s::Float64 = 0.0
    elapsed_sq_sum_s::Float64 = 0.0
end

Base.@kwdef mutable struct OuterRouteState
    lock::ReentrantLock = ReentrantLock()
    history::Dict{String, Dict{Symbol, OuterRouteStats}} = Dict{String, Dict{Symbol, OuterRouteStats}}()
end

function reset_outer_route_state!(state::OuterRouteState)
    lock(state.lock) do
        empty!(state.history)
    end
    return nothing
end

@inline function _route_stats_payload(stats::OuterRouteStats)::Dict{String, Any}
    return Dict{String, Any}(
        "samples" => max(0, Int(stats.samples)),
        "successes" => max(0, Int(stats.successes)),
        "failures" => max(0, Int(stats.failures)),
        "elapsed_sum_s" => max(0.0, Float64(stats.elapsed_sum_s)),
        "elapsed_sq_sum_s" => max(0.0, Float64(stats.elapsed_sq_sum_s))
    )
end

@inline function _route_payload_stats(payload)::Union{Nothing, OuterRouteStats}
    payload isa AbstractDict || return nothing
    samples = try
        max(0, Int(get(payload, "samples", 0)))
    catch
        0
    end
    successes = try
        max(0, Int(get(payload, "successes", 0)))
    catch
        0
    end
    failures = try
        max(0, Int(get(payload, "failures", 0)))
    catch
        0
    end
    elapsed_sum_s = try
        max(0.0, Float64(get(payload, "elapsed_sum_s", 0.0)))
    catch
        0.0
    end
    elapsed_sq_sum_s = try
        max(0.0, Float64(get(payload, "elapsed_sq_sum_s", NaN)))
    catch
        NaN
    end
    samples <= 0 && return nothing
    successes = min(samples, successes)
    failures = min(samples - successes, failures)
    elapsed_sum_s = max(0.0, elapsed_sum_s)
    if !isfinite(elapsed_sq_sum_s)
        # Backward-compatibility path for legacy state files without second-moment data.
        elapsed_sq_sum_s = (elapsed_sum_s^2) / samples
    end
    elapsed_sq_sum_s = max(elapsed_sq_sum_s, (elapsed_sum_s^2) / samples)
    return OuterRouteStats(
        samples=samples,
        successes=successes,
        failures=failures,
        elapsed_sum_s=elapsed_sum_s,
        elapsed_sq_sum_s=elapsed_sq_sum_s
    )
end

function save_outer_route_state(
    state::OuterRouteState,
    path::AbstractString;
    metadata::AbstractDict=Dict{String, Any}()
)::NamedTuple{(:path, :signatures, :rows), Tuple{String, Int, Int}}
    path_s = String(path)
    rows = Dict{String, Any}[]
    signatures = 0
    lock(state.lock) do
        for signature in sort!(collect(keys(state.history)))
            bucket = state.history[signature]
            isempty(bucket) && continue
            signature_rows = 0
            for route in (:none, :threads, :process)
                stats = get(bucket, route, nothing)
                stats isa OuterRouteStats || continue
                stats.samples > 0 || continue
                push!(rows, Dict{String, Any}(
                    "signature" => signature,
                    "route" => String(route),
                    "stats" => _route_stats_payload(stats)
                ))
                signature_rows += 1
            end
            signature_rows > 0 && (signatures += 1)
        end
    end

    metadata_out = Dict{String, Any}()
    for (k, v) in metadata
        key = String(k)
        metadata_out[key] = if v isa Number || v isa Bool
            v
        else
            String(v)
        end
    end

    payload = Dict{String, Any}(
        "schema_version" => 2,
        "updated_utc" => string(now(UTC)),
        "metadata" => metadata_out,
        "history" => rows
    )

    mkpath(dirname(path_s))
    open(path_s, "w") do io
        TOML.print(io, payload)
    end
    return (path=path_s, signatures=signatures, rows=length(rows))
end

function load_outer_route_state!(
    state::OuterRouteState,
    path::AbstractString;
    replace::Bool=true
)::NamedTuple{(:path, :signatures, :rows), Tuple{String, Int, Int}}
    path_s = String(path)
    isfile(path_s) || return (path=path_s, signatures=0, rows=0)
    parsed = TOML.parsefile(path_s)
    rows_in = get(parsed, "history", Any[])
    rows_in isa AbstractVector || return (path=path_s, signatures=0, rows=0)

    loaded_rows = 0
    loaded_signatures = Set{String}()
    lock(state.lock) do
        replace && empty!(state.history)
        for row in rows_in
            row isa AbstractDict || continue
            signature = strip(String(get(row, "signature", "")))
            isempty(signature) && continue
            route = Symbol(String(get(row, "route", "")))
            route in (:none, :threads, :process) || continue
            stats = _route_payload_stats(get(row, "stats", nothing))
            stats === nothing && continue

            bucket = get!(state.history, signature) do
                Dict{Symbol, OuterRouteStats}()
            end
            existing = get!(bucket, route) do
                OuterRouteStats()
            end
            existing.samples += stats.samples
            existing.successes += stats.successes
            existing.failures += stats.failures
            existing.elapsed_sum_s += stats.elapsed_sum_s
            existing.elapsed_sq_sum_s += stats.elapsed_sq_sum_s
            loaded_rows += 1
            push!(loaded_signatures, signature)
        end
    end
    return (path=path_s, signatures=length(loaded_signatures), rows=loaded_rows)
end

