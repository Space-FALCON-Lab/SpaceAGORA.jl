module Store

using Dates
using SHA
using TOML
using Arrow
using Tables

using ..Spec: CalibrationSpec, spec_to_dict
using ..ParamSpace: Candidate
using ..Backend: BackendEvaluation

export RunStore, LedgerEntry, RunState
export init_store, append_evaluation!, load_ledger_entries, load_stage_entries
export save_state!, load_state, stage_is_completed!
export spec_path, evaluations_path, state_path, best_manifest_path, report_path
export STAGE_SEQUENCE

const STAGE_SEQUENCE = (
    "prepare",
    "global_search_quick",
    "local_refine_full",
    "robustness_validation",
    "promote"
)

Base.@kwdef struct RunStore
    run_id::String
    run_dir::String
    spec_path::String
    evaluations_path::String
    state_path::String
    best_manifest_path::String
    report_path::String
end

Base.@kwdef struct LedgerEntry
    timestamp_utc::String
    stage::String
    candidate::Candidate
    evaluation::BackendEvaluation
end

Base.@kwdef mutable struct RunState
    schema_version::Int = 1
    run_id::String = ""
    current_stage::String = "prepare"
    next_candidate_id::Int = 1
    stage_status::Dict{String, String} = Dict{String, String}()
    best_candidate_id::Int = 0
    best_score::Float64 = Inf
    best_candidate_values::Dict{String, Any} = Dict{String, Any}()
    updated_utc::String = ""
end

@inline spec_path(store::RunStore)::String = store.spec_path
@inline evaluations_path(store::RunStore)::String = store.evaluations_path
@inline state_path(store::RunStore)::String = store.state_path
@inline best_manifest_path(store::RunStore)::String = store.best_manifest_path
@inline report_path(store::RunStore)::String = store.report_path

@inline function _safe_token(raw::AbstractString)::String
    token = replace(lowercase(strip(raw)), r"[^a-z0-9_\-]+" => "_")
    return isempty(token) ? "run" : token
end

@inline _sha1_hex(bytes::Vector{UInt8}) = bytes2hex(sha1(bytes))
@inline _sha1_hex(text::AbstractString) = _sha1_hex(Vector{UInt8}(codeunits(text)))

function _canonical_repr(x)::String
    if x isa AbstractDict
        keys_sorted = sort!(collect(keys(x)); by=k -> string(k))
        parts = String[]
        for k in keys_sorted
            push!(parts, string(repr(string(k)), "=>", _canonical_repr(x[k])))
        end
        return "{" * join(parts, ",") * "}"
    elseif x isa AbstractVector
        return "[" * join((_canonical_repr(v) for v in x), ",") * "]"
    elseif x isa Tuple
        return "(" * join((_canonical_repr(v) for v in x), ",") * ")"
    end
    return repr(x)
end

@inline function _hash_file(path::AbstractString)::String
    if !isfile(path)
        return "missing"
    end
    return _sha1_hex(read(path))
end

function _find_git_root(start::AbstractString)::Union{Nothing, String}
    path = abspath(start)
    while true
        if isdir(joinpath(path, ".git"))
            return path
        end
        parent = dirname(path)
        parent == path && return nothing
        path = parent
    end
end

function _git_sha(git_root::Union{Nothing, String})::String
    git_root === nothing && return "nogit"
    try
        return strip(readchomp(`git -C $(git_root) rev-parse HEAD`))
    catch
        return "nogit"
    end
end

function _dependency_lock_hash()::String
    lock_paths = String[]

    active = Base.active_project()
    if active !== nothing
        push!(lock_paths, joinpath(dirname(String(active)), "Manifest.toml"))
    end

    local_manifest = abspath(joinpath(@__DIR__, "..", "Manifest.toml"))
    push!(lock_paths, local_manifest)

    for p in lock_paths
        if isfile(p)
            return _hash_file(p)
        end
    end
    return "missing"
end

function _telemetry_hashes(spec::CalibrationSpec)::Vector{String}
    out = String[]
    for manifest_path in spec.manifest_paths
        abs_path = abspath(manifest_path)
        digest = _hash_file(abs_path)
        push!(out, string(abs_path, ":", digest))
    end
    sort!(out)
    return out
end

function _run_fingerprint(spec::CalibrationSpec)::String
    git_root = _find_git_root(abspath(joinpath(@__DIR__, "..", "..")))
    payload = Dict{String, Any}(
        "schema_version" => 1,
        "spec" => spec_to_dict(spec),
        "git_sha" => _git_sha(git_root),
        "dependency_lock_sha1" => _dependency_lock_hash(),
        "telemetry_hashes" => _telemetry_hashes(spec),
        "seeds" => Dict("seed" => spec.seed)
    )
    return _sha1_hex(_canonical_repr(payload))
end

@inline function _default_stage_status()::Dict{String, String}
    return Dict{String, String}(stage => "pending" for stage in STAGE_SEQUENCE)
end

@inline function _default_state(run_id::String)::RunState
    return RunState(
        run_id=run_id,
        current_stage=first(STAGE_SEQUENCE),
        next_candidate_id=1,
        stage_status=_default_stage_status(),
        updated_utc=string(now(UTC))
    )
end

function _state_to_dict(state::RunState)::Dict{String, Any}
    return Dict{String, Any}(
        "schema_version" => state.schema_version,
        "run_id" => state.run_id,
        "current_stage" => state.current_stage,
        "next_candidate_id" => state.next_candidate_id,
        "stage_status" => Dict(state.stage_status),
        "best_candidate_id" => state.best_candidate_id,
        "best_score" => state.best_score,
        "best_candidate_values" => Dict(state.best_candidate_values),
        "updated_utc" => state.updated_utc
    )
end

function _dict_to_state(doc::AbstractDict)::RunState
    stage_status = _default_stage_status()
    raw = get(doc, "stage_status", Dict{String, Any}())
    for (k, v) in pairs(raw)
        stage_status[String(k)] = String(v)
    end
    return RunState(
        schema_version=Int(get(doc, "schema_version", 1)),
        run_id=String(get(doc, "run_id", "")),
        current_stage=String(get(doc, "current_stage", "prepare")),
        next_candidate_id=Int(get(doc, "next_candidate_id", 1)),
        stage_status=stage_status,
        best_candidate_id=Int(get(doc, "best_candidate_id", 0)),
        best_score=Float64(get(doc, "best_score", Inf)),
        best_candidate_values=Dict{String, Any}(String(k) => v for (k, v) in pairs(get(doc, "best_candidate_values", Dict{String, Any}()))),
        updated_utc=String(get(doc, "updated_utc", ""))
    )
end

function init_store(spec::CalibrationSpec)::RunStore
    hash_hex = _run_fingerprint(spec)
    run_id = string(_safe_token(spec.id), "_", hash_hex[1:12])
    run_dir = joinpath(spec.output_root, "runs", run_id)
    mkpath(run_dir)

    store = RunStore(
        run_id=run_id,
        run_dir=run_dir,
        spec_path=joinpath(run_dir, "spec.toml"),
        evaluations_path=joinpath(run_dir, "evaluations.arrow"),
        state_path=joinpath(run_dir, "state.toml"),
        best_manifest_path=joinpath(run_dir, "best_manifest.toml"),
        report_path=joinpath(run_dir, "report.md")
    )

    frozen = spec_to_dict(spec)
    if isfile(store.spec_path)
        current = TOML.parsefile(store.spec_path)
        if _canonical_repr(current) != _canonical_repr(frozen)
            throw(ArgumentError("Existing run directory has non-matching spec.toml: $(store.spec_path)."))
        end
    else
        open(store.spec_path, "w") do io
            TOML.print(io, frozen)
        end
    end

    if !isfile(store.state_path)
        save_state!(store, _default_state(run_id))
    end

    return store
end

function save_state!(store::RunStore, state::RunState)::Nothing
    state.updated_utc = string(now(UTC))
    open(store.state_path, "w") do io
        TOML.print(io, _state_to_dict(state))
    end
    return nothing
end

function load_state(store::RunStore)::RunState
    if !isfile(store.state_path)
        return _default_state(store.run_id)
    end
    doc = TOML.parsefile(store.state_path)
    state = _dict_to_state(doc)
    isempty(state.run_id) && (state.run_id = store.run_id)
    return state
end

function stage_is_completed!(store::RunStore, state::RunState, stage::String; next_candidate_id::Int=state.next_candidate_id)
    state.stage_status[stage] = "completed"
    state.current_stage = stage
    state.next_candidate_id = next_candidate_id
    save_state!(store, state)
    return nothing
end

@inline function _dict_to_blob(d::Dict{String, Any})::String
    isempty(d) && return ""
    io = IOBuffer()
    TOML.print(io, Dict("payload" => d))
    return String(take!(io))
end

@inline function _blob_to_dict(text::String)::Dict{String, Any}
    isempty(strip(text)) && return Dict{String, Any}()
    parsed = TOML.parse(text)
    payload = get(parsed, "payload", Dict{String, Any}())
    return Dict{String, Any}(String(k) => v for (k, v) in pairs(payload))
end

@inline function _row_table(row::NamedTuple)
    return (
        timestamp_utc=[row.timestamp_utc],
        stage=[row.stage],
        candidate_id=[row.candidate_id],
        success=[row.success],
        score=[row.score],
        runtime_s=[row.runtime_s],
        error_message=[row.error_message],
        candidate_values=[row.candidate_values],
        metrics=[row.metrics],
        artifacts=[row.artifacts]
    )
end

function append_evaluation!(store::RunStore, candidate::Candidate, ev::BackendEvaluation)::Nothing
    row = (
        timestamp_utc=Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sss"),
        stage=String(ev.stage),
        candidate_id=Int64(ev.candidate_id),
        success=Bool(ev.success),
        score=Float64(ev.score),
        runtime_s=Float64(ev.runtime_s),
        error_message=String(ev.error_message),
        candidate_values=_dict_to_blob(Dict{String, Any}(candidate.values)),
        metrics=_dict_to_blob(Dict{String, Any}(ev.metrics)),
        artifacts=_dict_to_blob(Dict{String, Any}(ev.artifacts))
    )

    table = _row_table(row)
    if isfile(store.evaluations_path)
        Arrow.append(store.evaluations_path, table)
    else
        # Use stream format so we can append record batches incrementally.
        Arrow.write(store.evaluations_path, table; file=false)
    end
    return nothing
end

@inline function _metrics_dict(raw::Dict{String, Any})::Dict{String, Float64}
    out = Dict{String, Float64}()
    for (k, v) in raw
        fv = try
            Float64(v)
        catch
            continue
        end
        out[k] = fv
    end
    return out
end

@inline function _artifacts_dict(raw::Dict{String, Any})::Dict{String, String}
    return Dict{String, String}(k => String(v) for (k, v) in raw)
end

function load_ledger_entries(store::RunStore; stage::Union{Nothing, String}=nothing)::Vector{LedgerEntry}
    if !isfile(store.evaluations_path)
        return LedgerEntry[]
    end

    entries = LedgerEntry[]
    for row in Tables.rows(Arrow.Table(store.evaluations_path))
        row_stage = String(row.stage)
        if stage !== nothing && row_stage != stage
            continue
        end

        candidate = Candidate(
            id=Int(row.candidate_id),
            values=_blob_to_dict(String(row.candidate_values)),
            stage=row_stage
        )

        metrics = _metrics_dict(_blob_to_dict(String(row.metrics)))
        artifacts = _artifacts_dict(_blob_to_dict(String(row.artifacts)))

        ev = BackendEvaluation(
            candidate_id=Int(row.candidate_id),
            stage=row_stage,
            success=Bool(row.success),
            score=Float64(row.score),
            runtime_s=Float64(row.runtime_s),
            error_message=String(row.error_message),
            metrics=metrics,
            artifacts=artifacts
        )

        push!(entries, LedgerEntry(
            timestamp_utc=String(row.timestamp_utc),
            stage=row_stage,
            candidate=candidate,
            evaluation=ev
        ))
    end

    return entries
end

@inline load_stage_entries(store::RunStore, stage::String) = load_ledger_entries(store; stage=stage)

end # module Store
