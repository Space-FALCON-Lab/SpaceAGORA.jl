using ..SimulationModel.SimConfig: _with_fields
import ..SimulationModel
import ..SimulationEngine
using Arrow: Arrow
using CSV: CSV
using DataFrames: DataFrame, nrow, ncol, metadata

"""
    ConstellationEnsembleResult

What [`run_constellation_ensemble`](@ref)`(args; merge_results=true)` returns.

- `campaign::MonteCarloResult`: what the call returns without `merge_results`.
  Sample `k` is member `k`, and its `value` is what `run_simulation` returned for
  the keywords passed: `nothing` by default, the solution with
  `return_solution=true`, the `SimulationResults` with `return_results=true`.
- `members`: one row per member and saved time, a Tables.jl table (a
  `DataFrame`). It has the long layout's columns (see [`long_table`](@ref)):
  `spacecraft` (the member's position in `args.dynamics_model.spacecraft`, which
  is also its sample index and folder number), `spacecraft_id`, `plane` and
  `slot` (from its constellation slot, `missing` without one) and `time`, then
  the member's per-spacecraft fields without a `sc<i>_` prefix and its run-wide
  fields. Members with different fields give the union of the columns, in the
  order first seen, with `missing` where a member has none. Rows are in position
  order, then saved time. A failed member has no rows. When the ensemble writes
  `constellation_members.feather`, it writes it member by member and `members`
  reads it in place, with read-only columns; otherwise `members` is built in
  memory.
- `constellation::DataFrame`: one row per member, in position order:
  `spacecraft`, `spacecraft_id`, `plane`, `slot`, `member_tag` (the member's
  folder name), `success`, `error` (the error text of a failed member, `missing`
  otherwise), `elapsed_s`, `rows` (its rows in `members`) and the member's last
  row in `members` under the same names (`missing` for a member without rows).
- `files::Vector{String}`: the merged files written. With `results=true`,
  `constellation_members.feather` (with the results bundle, as the members
  write theirs), `constellation_members.csv` (with `save_csv`) and
  `constellation_summary.csv`, in `results_directory` beside the member
  folders. Empty otherwise.

Each member's rows come back from its run whole and are released once merged,
so the merge holds the members' rows about once, besides the values the members
return with `return_results=true`. For large campaigns,
`SPACEAGORA_RESULTS_THIN_STRIDE` and `SPACEAGORA_RESULTS_FINAL_ONLY` thin each
member's rows, and [`read_constellation_results`](@ref) merges member folders
already written, one member at a time.
"""
struct ConstellationEnsembleResult{M}
    campaign::MonteCarloResult
    members::M
    constellation::DataFrame
    files::Vector{String}
end

const _CONSTELLATION_IDENTITY = ("spacecraft", "spacecraft_id", "plane", "slot")

# A run-wide save field and a per-spacecraft one with the same column prefix
# would give two columns of one name in the members table, where a member's
# run-wide fields sit beside its per-spacecraft ones. Refused before any member
# runs; names that only meet after flattening are refused when the tables merge.
function _validate_merged_save_fields(run_kwargs)
    fields = get(run_kwargs, :save_fields, nothing)
    fields === nothing && return nothing
    per_spacecraft = Set{String}(field.column_prefix for field in fields if field.per_satellite)
    for field in fields
        field.per_satellite && continue
        field.column_prefix in per_spacecraft && throw(ArgumentError(
            "run_constellation_ensemble with merge_results=true puts each member's run-wide fields beside its " *
            "per-spacecraft fields, so they need different names, but the run-wide save field :$(field.name) " *
            "has the column prefix \"$(field.column_prefix)\" of a per-spacecraft field. Give it another column_prefix."))
    end
    return nothing
end

# One member's run with merge_results=true. It runs with `return_results=true`,
# so its saved rows are kept, and its final results save puts them, in the long
# layout from the long builder, in the sink. It returns what `run_simulation`
# returns for the caller's own keywords (`value`) and those tables. The
# configuration that came back with the rows is dropped here unless the caller
# asked for the results, so it never travels from a process worker.
function _constellation_member_outcome(member::SimulationConfiguration, run_kwargs)
    wants_results = get(run_kwargs, :return_results, false) === true
    wants_solution = get(run_kwargs, :return_solution, false) === true
    kwargs = Dict{Symbol, Any}(pairs(run_kwargs))
    kwargs[:return_results] = true
    sink = Ref{Any}(nothing)
    result = Base.ScopedValues.with(SimulationEngine._LONG_TABLES_SINK => sink) do
        SimulationEngine.run_simulation(member; kwargs...)
    end
    tables = sink[]
    value = wants_results ? result : wants_solution ? result.solution : nothing
    return (value=value, spacecraft=tables.spacecraft, run=tables.run,
        results_directory=result.configuration.simulation_settings.results_directory)
end

# The member run of an ensemble: `run_simulation`, or with merge_results=true the
# outcome above.
_constellation_member_run(member::SimulationConfiguration, run_kwargs, merge_results::Bool) =
    merge_results ? _constellation_member_outcome(member, run_kwargs) :
        SimulationEngine.run_simulation(member; run_kwargs...)

# Member `position`'s rows in the members table: its spacecraft table, with
# `spacecraft` set to its position in the constellation, and its run-wide
# columns beside them. A member runs one spacecraft, so its run table has one
# row per row of its spacecraft table. When the member's returned value holds
# these tables (the long layout with `return_results=true`), they are copied
# first, so that the value stays what `run_simulation` returned.
function _constellation_member_rows(position::Integer, outcome)::DataFrame
    value = hasproperty(outcome, :value) ? outcome.value : nothing
    shared = value isa SimulationEngine.SimulationResults &&
        (value.table === outcome.spacecraft || (outcome.run !== nothing && value.run_table === outcome.run))
    rows = shared ? copy(outcome.spacecraft) : outcome.spacecraft
    rows[!, "spacecraft"] = fill(Int32(position), nrow(rows))
    run = outcome.run
    set_presence(table) = SimulationModel.IOOutputs._set_long_column_presence!(table,
        Dict(name => [Int(position)] for name in names(table)
             if !(name in _CONSTELLATION_IDENTITY) && name != "time"))
    run === nothing && return set_presence(rows)
    isequal(run[!, "time"], rows[!, "time"]) ||
        error("Member $(position)'s run table and spacecraft table have different saved times.")
    for name in names(run)
        name == "time" && continue
        name in names(rows) && throw(ArgumentError(
            "run_constellation_ensemble with merge_results=true cannot put member $(position)'s run-wide column " *
            "\"$(name)\" beside its per-spacecraft column of the same name. Give the run-wide save field another " *
            "column_prefix."))
        rows[!, name] = shared ? copy(run[!, name]) : run[!, name]
    end
    return set_presence(rows)
end

# The folder name of member `position` (see `_ensemble_member_configuration`):
# the unpadded name when the member resumed from a checkpoint written there.
function _constellation_member_tag(member::SimulationConfiguration, position::Integer, count::Integer, id::Integer)
    unpadded = "sat_$(position)_id_$(id)"
    basename(normpath(member.simulation_settings.results_directory)) == unpadded && return unpadded
    return _ensemble_member_tag(position, count, id)
end

# The merged result of an ensemble run with merge_results=true, from the
# campaign whose samples hold the members' outcomes. Members are merged in
# position order, whatever order they finished in. As each member's rows are
# taken, its sample in `campaign` is replaced by one holding only the caller's
# value, so the members' tables are held once: in `member_rows`, from which
# `_write_constellation_results!` writes or appends them one at a time and
# releases each.
function _constellation_ensemble_result(args::SimulationConfiguration, member_configs, campaign::MonteCarloResult)
    spacecraft = args.dynamics_model.spacecraft
    count = length(spacecraft)
    slot_part(sc, part) = sc.constellation_slot === nothing ? missing : Int16(getproperty(sc.constellation_slot, part))
    member_rows = DataFrame[]
    presence = Dict{String, Vector{Int}}()
    last_rows = Union{Nothing, DataFrame}[]
    row_counts = Int[]
    destinations = String[]
    for (k, sample) in enumerate(campaign.samples)
        position = sample.index
        if sample.success
            outcome = sample.value
            destination = hasproperty(outcome, :results_directory) ? outcome.results_directory :
                member_configs[position].simulation_settings.results_directory
            push!(destinations, relpath(destination, member_configs[position].simulation_settings.results_directory))
            rows = _constellation_member_rows(position, outcome)
            for name in names(rows)
                (name in _CONSTELLATION_IDENTITY || name == "time") && continue
                push!(get!(presence, name, Int[]), Int(position))
            end
            push!(row_counts, nrow(rows))
            push!(last_rows, nrow(rows) == 0 ? nothing : rows[nrow(rows):nrow(rows), :])
            push!(member_rows, rows)
            # The sample now holds only the caller's value, here and in
            # `successful`, so only `member_rows` keeps this member's tables.
            released = _with_fields(sample; value=outcome.value)
            campaign.samples[k] = released
            j = findfirst(s -> s.index == position, campaign.successful)
            j === nothing || (campaign.successful[j] = released)
        else
            push!(destinations, ".")
            push!(row_counts, 0)
            push!(last_rows, nothing)
        end
    end


    positions = [sample.index for sample in campaign.samples]
    identity = DataFrame(
        spacecraft=Int32.(positions),
        spacecraft_id=Int64[spacecraft[k].id for k in positions],
        plane=Union{Missing, Int16}[slot_part(spacecraft[k], :plane) for k in positions],
        slot=Union{Missing, Int16}[slot_part(spacecraft[k], :slot) for k in positions],
        member_tag=String[_constellation_member_tag(member_configs[k], k, count, spacecraft[k].id) for k in positions],
        results_subdirectory=destinations,
    )
    constellation = _constellation_summary(identity, Bool[sample.success for sample in campaign.samples],
        Union{Missing, String}[sample.success ? missing : sprint(showerror, sample.error) for sample in campaign.samples],
        Float64[sample.elapsed_s for sample in campaign.samples], row_counts, last_rows)

    result_campaign = MonteCarloResult(campaign.samples, campaign.elapsed_s, campaign.threads;
                                       route=campaign.route, local_slots=campaign.local_slots)
    members, files = _write_constellation_results!(args, member_rows, constellation, presence)
    return ConstellationEnsembleResult(result_campaign, members, constellation, files)
end

# The constellation table: the members' identity and folder (`identity`), their
# status, their rows in the members table and their last row there, under the
# same names (`missing` for a member without rows).
function _constellation_summary(identity::DataFrame, success, error, elapsed_s, row_counts, last_rows; read_status=nothing)
    constellation = copy(identity)
    constellation[!, "success"] = success
    constellation[!, "error"] = error
    constellation[!, "elapsed_s"] = elapsed_s
    constellation[!, "rows"] = row_counts
    read_status === nothing || (constellation[!, "read_status"] = read_status)
    last_names = String[]
    for last in last_rows
        last === nothing && continue
        for name in names(last)
            name in _CONSTELLATION_IDENTITY || name in last_names || push!(last_names, name)
        end
    end
    for name in last_names
        name in names(constellation) && throw(ArgumentError(
            "merged constellation summary reserves column '$name' for campaign metadata; " *
            "choose a different SaveField column_prefix"))
        constellation[!, name] = [last === nothing || !(name in names(last)) ? missing : last[1, name]
                                  for last in last_rows]
    end
    return constellation
end

# The members table of an ensemble without rows: its identity columns.
_empty_members() = DataFrame(spacecraft=Int32[], spacecraft_id=Int64[], plane=Union{Missing, Int16}[],
    slot=Union{Missing, Int16}[], time=Float64[])

# The members table's columns from the members' tables, in the order first seen,
# each with a type that holds every member's values and `missing` where a member
# lacks the column, as `append!(...; cols=:union, promote=true)` gives them.
function _members_columns(tables)
    order = String[]
    types = Dict{String, Type}()
    present = Dict{String, Int}()
    seen = 0
    for table in tables
        seen += 1
        for name in names(table)
            T = eltype(table[!, name])
            if haskey(types, name)
                types[name] = promote_type(types[name], T)
                present[name] += 1
            else
                push!(order, name)
                types[name] = T
                present[name] = 1
            end
        end
    end
    if seen == 0
        empty = _empty_members()
        return Pair{String, Type}[name => eltype(empty[!, name]) for name in names(empty)]
    end
    return Pair{String, Type}[name => (present[name] < seen ? Union{Missing, types[name]} : types[name]) for name in order]
end

# A member's rows in the members table's `columns`.
_members_batch(rows::DataFrame, columns) = DataFrame([name => (name in names(rows) ?
    convert(Vector{T}, rows[!, name]) : Vector{T}(fill(missing, nrow(rows)))) for (name, T) in columns]; copycols=false)

# The members table as an Arrow file at `path`, written member by member, one
# record batch each, so that a member's rows are in memory only while they are
# written. `tables` gives the members' tables in position order, one at a time.
function _write_members_feather!(path::String, columns, tables, presence)
    return SimulationModel.IOSerialization._atomic_write_file(path, tmp -> begin
        schema = SimulationModel.IOOutputs._set_long_column_presence!(_empty_members(), presence)
        key = SimulationModel.IOOutputs._LONG_COLUMN_SPACECRAFT_KEY
        writer = open(Arrow.Writer, tmp; ntasks=1, metadata=Dict(key => metadata(schema, key)))
        try
            written = 0
            for rows in tables
                Arrow.write(writer, _members_batch(rows, columns))
                written += 1
            end
            written == 0 && Arrow.write(writer, _members_batch(_empty_members(), columns))
        finally
            close(writer)
        end
    end)
end

# The merged files, in `results_directory` beside the member folders, when the
# ensemble saves results: the members table as the members write theirs (the
# feather with the results bundle, written member by member, and then read in
# place; the CSV with `save_csv`), and the constellation table as
# `constellation_summary.csv`. Without the feather the members table is built in
# memory. Returns the members table and the files written.
function _write_constellation_results!(args::SimulationConfiguration, member_rows::Vector{DataFrame}, constellation::DataFrame, presence)
    settings = args.simulation_settings
    directory = settings.results_directory
    write(path, writer) = SimulationModel.IOSerialization._atomic_write_file(path, writer)
    files = String[]
    if settings.results && SimulationEngine._typed_save_bundle_enabled()
        mkpath(directory)
        path = joinpath(directory, "constellation_members.feather")
        columns = _members_columns(member_rows)
        # Each member's rows are released once written.
        pending = Iterators.map(_ -> popfirst!(member_rows), 1:length(member_rows))
        push!(files, _write_members_feather!(path, columns, pending, presence))
        members = DataFrame(Arrow.Table(path); copycols=false)
    else
        members = isempty(member_rows) ? _empty_members() : popfirst!(member_rows)
        while !isempty(member_rows)
            append!(members, popfirst!(member_rows); cols=:union, promote=true)
        end
    end
    SimulationModel.IOOutputs._set_long_column_presence!(members, presence)
    settings.results || return members, files
    mkpath(directory)
    settings.save_csv &&
        push!(files, write(joinpath(directory, "constellation_members.csv"), tmp -> CSV.write(tmp, members)))
    push!(files, write(joinpath(directory, "constellation_summary.csv"), tmp -> CSV.write(tmp, constellation)))
    return members, files
end

const _MEMBER_FOLDER = r"^sat_(\d+)_id_(-?\d+)$"

"""
    read_constellation_results(results_directory; overwrite=false) -> (members, constellation, files)

Merge the member folders a constellation ensemble wrote under
`results_directory` (`sat_<i>_id_<id>`, each with its results bundle) into the
tables of [`ConstellationEnsembleResult`](@ref), reading one member at a time.
Nothing is run again.

- `members` has the columns and rows `run_constellation_ensemble(args;
  merge_results=true)` gives: `spacecraft` (the folder's `<i>`),
  `spacecraft_id`, `plane`, `slot` and `time`, then each member's
  per-spacecraft fields without a prefix and its run-wide fields. It is written
  to `constellation_members.feather`, one record batch per member, and read in
  place: only one member's rows are in memory at a time.
- `constellation` has one row per member, as the ensemble's `constellation`
  table, and a `read_status` column: `"read"`, or why the member's rows were not
  read (`"member failed"`, `"folder missing"`, `"no results bundle"`,
  `"unfinished run"`, `"files differ from manifest"`). The members are the
  folders found and, when `constellation_summary.csv` is there, the members it
  lists, so a member whose folder is gone keeps its row. Its `success`, `error`
  and `elapsed_s` come from that summary; without one they are `missing`. It is
  written to `constellation_summary.csv`.
- `files` lists the two files written.

The two files are not replaced when they exist already (an ensemble run with
`merge_results=true` wrote them): pass `overwrite=true` to replace them.

A member written in the long layout (`results_layout=:long`) gives its tables as
they are, with its position set. A member written in the wide layout has one
spacecraft, so its per-spacecraft columns are its `sc1_` ones, renamed without
the prefix; its plane and slot come from `constellation_summary.csv`, or are
`missing` without it. A member's bundle is read only if its manifest has a known
version and layout (`"1"` wide, `"2"` long) and its files have the sizes the
manifest records; a bundle of another version, or a wide bundle of more than one
spacecraft, is refused. A member that the summary records as failed, or whose
run streamed its results and did not finish, has no rows. Members not read are
named in a warning.
"""
function read_constellation_results(results_directory::AbstractString; overwrite::Bool=false)
    directory = String(results_directory)
    isdir(directory) || throw(ArgumentError("No results directory at $(directory)."))
    members_path = joinpath(directory, "constellation_members.feather")
    summary_path = joinpath(directory, "constellation_summary.csv")
    summary = isfile(summary_path) ? CSV.read(summary_path, DataFrame) : nothing
    overwrite || !(isfile(members_path) || summary !== nothing) || throw(ArgumentError(
        "$(directory) already holds merged constellation results (constellation_members.feather or " *
        "constellation_summary.csv). Pass overwrite=true to replace them with the ones read from the member folders."))
    # The members: every folder found and every member the summary lists.
    found = Dict{String, Tuple{Int, Int}}()
    for name in readdir(directory)
        m = match(_MEMBER_FOLDER, name)
        (m === nothing || !isdir(joinpath(directory, name))) && continue
        found[name] = (parse(Int, m[1]), parse(Int, m[2]))
    end
    listed = Dict{String, Any}()
    if summary !== nothing
        for row in eachrow(summary)
            listed[String(row.member_tag)] = row
        end
    end
    tags = union(keys(found), keys(listed))
    isempty(tags) && throw(ArgumentError("No constellation member folders (sat_<i>_id_<id>) under $(directory)."))
    position_of(tag) = haskey(found, tag) ? found[tag] : (Int(listed[tag].spacecraft), Int(listed[tag].spacecraft_id))
    members_list = sort!([(position_of(tag)..., tag) for tag in tags])
    as_int16(v) = v === nothing || ismissing(v) ? missing : Int16(v)
    recorded(tag, name) = haskey(listed, tag) && Symbol(name) in propertynames(listed[tag]) ? listed[tag][name] : missing
    member_folder(tag) = SimulationModel.IOSerialization._member_results_directory(
        joinpath(directory, tag), recorded(tag, "results_subdirectory"))

    # Why a member's rows are not read, or "read".
    function status_of(member)
        _, _, tag = member
        haskey(listed, tag) && !coalesce(listed[tag].success, true) && return "member failed"
        haskey(found, tag) || return "folder missing"
        return _disk_member_status(member_folder(tag))
    end
    rows_of(member) = begin
        position, id, tag = member
        located = (plane=as_int16(recorded(tag, "plane")), slot=as_int16(recorded(tag, "slot")))
        _disk_member_rows(member_folder(tag), position, id, located)
    end
    statuses = String[status_of(member) for member in members_list]
    readable = [status == "read" for status in statuses]
    unread = [(member[3], status) for (member, status) in zip(members_list, statuses) if status ∉ ("read", "member failed")]
    isempty(unread) || @warn "read_constellation_results did not read the rows of $(length(unread)) member(s): " *
        join(("$(tag) ($(status))" for (tag, status) in unread), ", ") * "."

    # A first pass for the columns, the row counts and the last rows; the second
    # writes the members, reading each again.
    row_counts = zeros(Int, length(members_list))
    last_rows = Union{Nothing, DataFrame}[nothing for _ in members_list]
    planes = Union{Missing, Int16}[as_int16(recorded(member[3], "plane")) for member in members_list]
    slots = Union{Missing, Int16}[as_int16(recorded(member[3], "slot")) for member in members_list]
    schema_tables = DataFrame[]
    presence = Dict{String, Vector{Int}}()
    for (k, member) in enumerate(members_list)
        readable[k] || continue
        rows = rows_of(member)
        for name in names(rows)
            (name in _CONSTELLATION_IDENTITY || name == "time") && continue
            push!(get!(presence, name, Int[]), member[1])
        end
        row_counts[k] = nrow(rows)
        if nrow(rows) > 0
            last_rows[k] = DataFrame(rows[nrow(rows):nrow(rows), :])
            planes[k], slots[k] = rows[1, "plane"], rows[1, "slot"]
        end
        # Only the names and types are kept.
        push!(schema_tables, rows[1:0, :])
    end
    columns = _members_columns(schema_tables)
    identity = DataFrame(
        spacecraft=Int32[member[1] for member in members_list],
        spacecraft_id=Int64[member[2] for member in members_list],
        plane=planes,
        slot=slots,
        member_tag=String[member[3] for member in members_list],
        results_subdirectory=String[ismissing(recorded(member[3], "results_subdirectory")) ? "." :
            recorded(member[3], "results_subdirectory") for member in members_list],
    )
    constellation = _constellation_summary(identity,
        Union{Missing, Bool}[recorded(member[3], "success") for member in members_list],
        Union{Missing, String}[recorded(member[3], "error") for member in members_list],
        Union{Missing, Float64}[recorded(member[3], "elapsed_s") for member in members_list],
        row_counts, last_rows; read_status=statuses)
    # Validate the complete summary before replacing either existing merged file.
    _write_members_feather!(members_path, columns,
        (rows_of(member) for (k, member) in enumerate(members_list) if readable[k]), presence)
    members = DataFrame(Arrow.Table(members_path); copycols=false)
    SimulationModel.IOOutputs._set_long_column_presence!(members, presence)
    written = SimulationModel.IOSerialization._atomic_write_file(summary_path, tmp -> CSV.write(tmp, constellation))
    return (members=members, constellation=constellation, files=[members_path, written])
end

# Whether a member folder's results bundle can be read ("read"), or why not. A
# bundle of an unknown version or layout, or a wide bundle of more than one
# spacecraft, is refused.
function _disk_member_status(folder::String)::String
    manifest_path = joinpath(folder, "simulation_results.manifest.toml")
    isfile(manifest_path) || return "no results bundle"
    manifest = TOML.parsefile(manifest_path)
    get(manifest, "streaming", false) === true && return "unfinished run"
    version, layout = get(manifest, "schema_version", nothing), get(manifest, "layout", "wide")
    (version, layout) in (("1", "wide"), ("2", "long")) || throw(ArgumentError(
        "$(folder) holds a results bundle of schema version $(repr(version)) and layout $(repr(layout)); " *
        "read_constellation_results reads version \"1\" (wide) and \"2\" (long)."))
    count = get(manifest, "spacecraft_count", 1)
    layout == "wide" && count != 1 && throw(ArgumentError(
        "$(folder) holds a wide results bundle of $(count) spacecraft; a constellation member runs one."))
    parts = layout == "long" ? ("spacecraft_feather", "run_feather") : ("feather",)
    files = get(manifest, "files", Dict{String, Any}())
    layout == "long" && !haskey(files, "spacecraft_feather") && return "no results bundle"
    for part in parts
        haskey(files, part) || continue
        path = joinpath(folder, basename(files[part]["path"]))
        isfile(path) && filesize(path) == files[part]["size_bytes"] || return "files differ from manifest"
    end
    layout == "wide" && !isfile(joinpath(folder, "simulation_results.feather")) && return "no results bundle"
    return "read"
end

# A member folder's rows in the members table, from its results bundle, which
# `_disk_member_status` found readable. `located` gives the plane and slot of a
# member written in the wide layout.
function _disk_member_rows(folder::String, position::Int, id::Int, located)
    manifest = TOML.parsefile(joinpath(folder, "simulation_results.manifest.toml"))
    if get(manifest, "layout", "wide") == "long"
        run_path = joinpath(folder, "simulation_results_run.feather")
        spacecraft = DataFrame(Arrow.Table(joinpath(folder, "simulation_results_spacecraft.feather")); copycols=false)
        run = haskey(manifest["files"], "run_feather") ? DataFrame(Arrow.Table(run_path); copycols=false) : nothing
        return _constellation_member_rows(position, (spacecraft=spacecraft, run=run))
    end
    wide = DataFrame(Arrow.Table(joinpath(folder, "simulation_results.feather")); copycols=false)
    n = nrow(wide)
    spacecraft = DataFrame(spacecraft=fill(Int32(position), n), spacecraft_id=fill(Int64(id), n),
        plane=Vector{Union{Missing, Int16}}(fill(located.plane, n)),
        slot=Vector{Union{Missing, Int16}}(fill(located.slot, n)), time=wide[!, "time"]; copycols=false)
    run = DataFrame(time=wide[!, "time"]; copycols=false)
    SimulationModel.IOOutputs._check_long_column_names([name[5:end] for name in names(wide) if startswith(name, "sc1_")])
    for name in names(wide)
        name == "time" && continue
        startswith(name, "sc1_") ? (spacecraft[!, name[5:end]] = wide[!, name]) : (run[!, name] = wide[!, name])
    end
    return _constellation_member_rows(position, (spacecraft=spacecraft, run=ncol(run) > 1 ? run : nothing))
end
