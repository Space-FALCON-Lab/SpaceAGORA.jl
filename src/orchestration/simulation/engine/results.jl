"""
    SimulationResults

What `run_simulation(args; return_results=true)` returns.

- `table::DataFrame`: the results table, same columns and values as the results
  CSV, built in memory even when `simulation_settings.results=false`. A resumed
  run's table holds the whole run, from the rows kept beside its checkpoint.
  With `results_layout=:long` it is the spacecraft table: one row per spacecraft
  and saved time (see [`long_table`](@ref)). A long-layout run that streams its
  results at each checkpoint reads it from its Feather file in place, with
  read-only columns; `DataFrame(results.table)` copies it.
- `configuration`: the configuration that ran. Under the default
  `isolate_state=true` this is the deep copy, so mutated controller state and
  logs (for example the RPO command log) are reachable here and not on the
  caller's object. A resumed run whose checkpoint holds no usable rows writes
  to `resumed_at_<t>s`, and its configuration's `results_directory` names that
  folder.
- `files::Vector{String}`: paths of existing files the run wrote (results CSV,
  bundle, scene, checkpoint); empty when none.
- `solution`: the ODE solution when `return_solution=true` was also passed,
  otherwise `nothing`.
- `run_table`: with `results_layout=:long`, the run table (one row per saved
  time, `time` and the run-wide fields), or `nothing` when the run has none;
  `nothing` in the wide layout. Read it with [`run_table`](@ref).
"""
struct SimulationResults{C, S}
    table::DataFrame
    configuration::C
    files::Vector{String}
    solution::S
    run_table::Union{Nothing, DataFrame}
end

SimulationResults(table::DataFrame, configuration, files::Vector{String}, solution) =
    SimulationResults(table, configuration, files, solution, nothing)

# Set by `run_constellation_ensemble(merge_results=true)` around each member's
# run. The run's final results save puts its tables in the long layout there,
# built by the long builder from the saved rows, whichever layout the run writes
# and returns.
const _LONG_TABLES_SINK = Base.ScopedValues.ScopedValue{Union{Nothing, Base.RefValue{Any}}}(nothing)

"""
    run_table(results::SimulationResults) -> Union{Nothing, DataFrame}

The run table of a run with `results_layout=:long`: one row per saved time,
`time` and the run-wide fields (such as `sun_dir_1..3`). `nothing` for a wide
run, whose table holds those columns, and for a long run without run-wide
fields.
"""
run_table(results::SimulationResults) = results.run_table

"""
    long_table(results::SimulationResults) -> DataFrame

The spacecraft table in the long layout: one row per spacecraft and saved time,
with `spacecraft` (its position in the configuration), `spacecraft_id`, `plane`
and `slot` (from its constellation slot, `missing` without one) and `time` first,
then each per-spacecraft column without its `sc<i>_` prefix. Rows are
time-major. A run with `results_layout=:long` returns its table as it is. For a
wide run the table is rearranged from the `sc<i>_` columns of its spacecraft; a
column that some spacecraft lack (thruster levels of a spacecraft without
thrusters) is `missing` for them. Its run-wide columns stay in `results.table`.
"""
function long_table(results::SimulationResults)::DataFrame
    results.configuration.simulation_settings.results_layout === :long && return results.table
    wide = results.table
    spacecraft = results.configuration.dynamics_model.spacecraft
    num_sats = length(spacecraft)
    n_times = nrow(wide)
    prefix(i) = "sc$(i)_"
    # Every spacecraft's columns, in the order first seen: spacecraft can differ.
    fields = String[]
    for name in names(wide)
        m = match(r"^sc(\d+)_(.*)$", name)
        m === nothing || !(1 <= parse(Int, m[1]) <= num_sats) || m[2] in fields || push!(fields, m[2])
    end
    SimulationModel.IOOutputs._check_long_column_names(fields)
    slot_part(sc, part) = sc.constellation_slot === nothing ? missing : Int16(getproperty(sc.constellation_slot, part))
    long = DataFrame(
        spacecraft=repeat(Int32.(1:num_sats), n_times),
        spacecraft_id=repeat(Int64[sc.id for sc in spacecraft], n_times),
        plane=repeat(Union{Missing, Int16}[slot_part(sc, :plane) for sc in spacecraft], n_times),
        slot=repeat(Union{Missing, Int16}[slot_part(sc, :slot) for sc in spacecraft], n_times),
        time=repeat(Vector{Float64}(wide.time); inner=num_sats),
    )
    presence = Dict{String, Vector{Int}}()
    for field in fields
        presence[field] = [i for i in 1:num_sats if prefix(i) * field in names(wide)]
        columns = [prefix(i) * field in names(wide) ? wide[!, prefix(i) * field] : nothing for i in 1:num_sats]
        long[!, field] = [columns[sat] === nothing ? missing : columns[sat][row] for row in 1:n_times for sat in 1:num_sats]
    end
    return SimulationModel.IOOutputs._set_long_column_presence!(long, presence)
end

# The link schedule's columns before any term columns, in order, with their types.
const _LINK_SCHEDULE_COLUMNS = (:window => String, :link => Int, :kind => String, :directed => Bool,
    :a_spacecraft => Int, :a_terminal => Int, :b_spacecraft => Int, :b_terminal => Int,
    :start_s => Float64, :end_s => Float64, :duration_s => Float64, :steps => Int,
    :score_start => Float64, :score_end => Float64, :score_min => Float64, :score_max => Float64)

"""
    link_schedule(model::InterlinkModel; policy=nothing) -> DataFrame
    link_schedule(results::SimulationResults) -> DataFrame

The link schedule of a run, from its interlink history (`model.history`): one
row for each span of accepted steps over which a link was selected
(`window == "active"`) and, when the history records the available links
(`InterlinkModel(; record_available)`), one for each span over which a link was
available (`window == "contact"`). Rows are sorted by start time, the active
rows first, then by link key. A run that writes results writes the table beside
them, as `simulation_results_link_schedule.csv` and `.feather`.

| Column | Unit | Meaning |
| --- | --- | --- |
| `window` | | `"active"` or `"contact"` |
| `link` | | the link's number, its position in registration order |
| `kind` | | the link kind's type name with its type parameters, without modules, for example `DataLink{ShannonBudget}` |
| `directed` | | whether the link is directed |
| `a_spacecraft`, `a_terminal` | | the first endpoint of the link's key, the transmitter of a directed link; a ground station end is `-s` with its antenna number |
| `b_spacecraft`, `b_terminal` | | the second endpoint |
| `a_name`, `b_name` | | only for station runs: the station name or `spacecraft N` |
| `start_s`, `end_s`, `duration_s` | s | the span: its first sample, the next sample without the link or the history's last sample, and their difference |
| `steps` | | the history samples in the span |
| `score_start`, `score_end`, `score_min`, `score_max` | the policy's | the score that selected the link at the span's first and last sample, and its extremes; `NaN` for a contact row and for samples without scores |
| `term1_start`, `term2_start`, ... | each term's | with `InterlinkModel(; record_terms=true)`: each term's score at the span's first sample |

The scores are in the policy's unit: m/s for `SemimajorAxisRatePolicy`, 1/s for
`EccentricityRatePolicy`, W for `PowerDeliveryPolicy` and
`EnergyBalancingPolicy`, bit/s for `DataThroughputPolicy`, and 1 for
`WeightedPolicy`. The spans are those of `contact_windows`, at the scheduler's
resolution, the run's accepted steps.

The table's metadata:
- `policy` and `score_unit`, from `policy` (the run's scheduler or its policy);
- `history_start_s` and `history_end_s`, the times of the history's first and
  last samples (NaN without samples): the span the schedule covers. A resumed
  run restores usable history from beside its checkpoint. Without usable history,
  a continuing run records only the resumed interval; a resume that solves nothing
  returns an empty schedule;
- `record_available`, whether the history records the available links, so that
  no contact row means none was available rather than none was recorded.

With term columns, each `term<k>_start` column's metadata gives its term's
`policy`, `score_unit`, `weight` and `scale`.

The `results` form reads the configuration that ran. A run without interlinks
gives an empty table with the same columns. It reads the interlink model as it
is when called: with `isolate_state=false`, a later run on the same
configuration refills that model's history, and so changes this table; and a
resume that solves nothing returns any validated checkpoint history without advancing it.
"""
function link_schedule(model::SimulationModel.InterlinkModel; policy=nothing)
    history = model.history
    term_count = maximum(sample -> size(sample.terms, 1), history; init=0)
    rows = [_link_schedule_row(window_kind, window, history, model.graph, term_count)
            for (window_kind, active) in (("active", true), ("contact", false))
            for window in SimulationModel.InterLinkModels._link_windows(history; active)]
    sort!(rows; by=row -> (row.start_s, row.window == "active" ? 0 : 1,
        (row.a_spacecraft, row.a_terminal, row.b_spacecraft, row.b_terminal)))
    table = DataFrame([name => T[getproperty(row, name) for row in rows] for (name, T) in _LINK_SCHEDULE_COLUMNS])
    for k in 1:term_count
        table[!, "term$(k)_start"] = Float64[row.terms[k] for row in rows]
    end
    # Preserve spacecraft-only table schemas; station runs also expose readable endpoints.
    if !isempty(model.ground_stations)
        endpoint_name(node) = node < 0 ? model.ground_stations[-node].name : "spacecraft $(node)"
        table[!, :a_name] = endpoint_name.(table.a_spacecraft)
        table[!, :b_name] = endpoint_name.(table.b_spacecraft)
    end
    description = _link_schedule_policy(policy)
    metadata!(table, "policy", description.policy; style=:note)
    metadata!(table, "score_unit", description.score_unit; style=:note)
    metadata!(table, "history_start_s", isempty(history) ? NaN : first(history).time; style=:note)
    metadata!(table, "history_end_s", isempty(history) ? NaN : last(history).time; style=:note)
    metadata!(table, "record_available", SimulationModel.InterLinkModels.records_available(model); style=:note)
    for (k, term) in enumerate(description.terms)
        k <= term_count || break
        for (key, value) in pairs(term)
            colmetadata!(table, "term$(k)_start", String(key), value; style=:note)
        end
    end
    return table
end

function link_schedule(results::SimulationResults)
    configuration = results.configuration
    model = hasproperty(configuration, :interlink_model) ? configuration.interlink_model : nothing
    model === nothing && return link_schedule(SimulationModel.InterlinkModel())
    return link_schedule(model; policy=configuration.scheduling_policy_model)
end

function _link_schedule_row(window_kind::String, window, history, graph, term_count::Int)
    key = window.key
    edge = get(graph.edge_index_lookup, key, 0)
    link = edge == 0 ? nothing : graph.edge_properties[edge]
    (a_spacecraft, a_terminal), (b_spacecraft, b_terminal) = key
    active = window_kind == "active"
    scores = active ? Float64[_sample_score(history[index], key) for index in window.first:window.last] : [NaN]
    terms = ntuple(k -> active ? _sample_term(history[window.first], key, k) : NaN, term_count)
    return (window=window_kind, link=edge, kind=link === nothing ? "" : _link_kind_name(typeof(link)),
        directed=link !== nothing && SimulationModel.InterLinkModels.link_directed(link),
        a_spacecraft, a_terminal, b_spacecraft, b_terminal,
        start_s=window.start_s, end_s=window.end_s, duration_s=window.end_s - window.start_s,
        steps=window.last - window.first + 1,
        score_start=first(scores), score_end=last(scores), score_min=minimum(scores), score_max=maximum(scores),
        terms)
end

# The score that selected `key` in a history sample, NaN when the sample has none.
function _sample_score(sample, key)
    isempty(sample.scores) && return NaN
    return sample.scores[findfirst(==(key), sample.active)]
end

# Term `k`'s score for `key` in a history sample, NaN when the sample has no such term.
function _sample_term(sample, key, k::Int)
    size(sample.terms, 1) >= k || return NaN
    return sample.terms[k, findfirst(==(key), sample.active)]
end

# A link kind's name with its parameters, without its module.
function _link_kind_name(kind::DataType)
    isempty(kind.parameters) && return string(nameof(kind))
    parameter_name(parameter) = parameter isa DataType ? _link_kind_name(parameter) : string(parameter)
    return string(nameof(kind), "{", join(map(parameter_name, kind.parameters), ", "), "}")
end

# The policy's name and score unit, and for a weighted policy its terms: a run's
# scheduler or a policy, or nothing.
function _link_schedule_policy(policy)
    SP = SimulationModel.SchedulingPolicies
    objective = policy isa SP.SchedulingPolicyModel ? policy.objective : policy
    objective === nothing && return (policy="", score_unit="", terms=NamedTuple[])
    terms = objective isa SP.WeightedPolicy ?
        NamedTuple[(policy=string(nameof(typeof(term.policy))), score_unit=_score_unit(term.policy),
            weight=Float64(term.weight), scale=Float64(term.scale)) for term in objective.terms] : NamedTuple[]
    return (policy=string(nameof(typeof(objective))), score_unit=_score_unit(objective), terms=terms)
end

function _score_unit(objective)
    SP = SimulationModel.SchedulingPolicies
    objective isa SP.SemimajorAxisRatePolicy && return "m/s"
    objective isa SP.EccentricityRatePolicy && return "1/s"
    objective isa Union{SP.PowerDeliveryPolicy, SP.EnergyBalancingPolicy} && return "W"
    objective isa Union{SP.DataThroughputPolicy, SP.DataBacklogPolicy} && return "bit/s"
    objective isa SP.WeightedPolicy && return "1"
    return ""
end

# The link schedule's entry in the results manifest: the policy and its unit, the
# number of windows, the recording flags, the span the history covers (absent
# without samples) and, for a weighted mix, its terms.
function _link_schedule_manifest(table::DataFrame, policy, model)
    description = _link_schedule_policy(policy)
    entry = Dict{String, Any}("policy" => description.policy, "score_unit" => description.score_unit,
        "rows" => nrow(table), "record_available" => SimulationModel.InterLinkModels.records_available(model),
        "record_terms" => model.record_terms)
    history = model.history
    isempty(history) || (entry["history_start_s"] = first(history).time; entry["history_end_s"] = last(history).time)
    isempty(description.terms) || (entry["terms"] = [Dict{String, Any}(String(name) => value
        for (name, value) in pairs(term)) for term in description.terms])
    return entry
end
