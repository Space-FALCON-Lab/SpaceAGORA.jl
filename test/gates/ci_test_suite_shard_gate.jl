# The default entrypoint is sharded across CI jobs. That is only safe while the
# shards named in the workflow cover every suite exactly once: a suite in no
# shard simply stops running, and every shard stays green while it does. A suite
# in two shards only wastes time, but is still worth naming.
#
# This gate reads the shard definitions out of .github/workflows/julia-ci.yml and
# checks them against the suite list the runner actually knows about, so adding a
# suite without adding it to a shard fails here instead of going untested.

const SHARD_GATE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _shard_gate_known_suites()::Vector{String}
    runner = read(joinpath(SHARD_GATE_ROOT, "test", "integration", "runtests.jl"), String)
    block = match(r"const _ALL_SUITES = \[(.*?)\]"s, runner)
    block === nothing && error("Could not find _ALL_SUITES in test/integration/runtests.jl")
    return [m.captures[1] for m in eachmatch(r"(\d\d)_[a-z0-9_]+\.jl", block.captures[1])]
end

function _shard_gate_workflow()
    workflow = read(joinpath(SHARD_GATE_ROOT, ".github", "workflows", "julia-ci.yml"), String)
    # The shard entries live in the tests-matrix `include:` list, one per
    # `- shard:` item; `suites:` is a sibling key on the following line, so this
    # matches the key wherever it is indented rather than assuming it leads the
    # list item.
    suites = [m.captures[1] for m in eachmatch(r"^\s+suites: \"([0-9,]+)\""m, workflow)]
    isempty(suites) && error("No `suites: \"...\"` shard entries found in julia-ci.yml")
    probes = [m.captures[1] for m in eachmatch(r"^\s+probe_shard: \"([0-9]+/[0-9]+)\""m, workflow)]
    skips = length(collect(eachmatch(r"^\s+skip_unit_driver: \"1\""m, workflow)))
    return suites, probes, skips
end

function ci_test_suite_shard_gate()
    known = _shard_gate_known_suites()
    isempty(known) && error("Parsed no suites from _ALL_SUITES")
    shard_suites, probe_shards, skip_unit_count = _shard_gate_workflow()

    covered = String[]
    for shard in shard_suites, token in split(shard, ",")
        tok = strip(token)
        isempty(tok) || push!(covered, lpad(tok, 2, '0'))
    end

    absent = setdiff(known, covered)
    isempty(absent) || error("Suites in no CI shard, so they would never run: " * join(absent, ", "))

    unknown = setdiff(covered, known)
    isempty(unknown) || error("CI shards name suites that do not exist: " * join(unknown, ", "))

    # Suite 09 is deliberately in more than one shard: its probe list is split by
    # SPACEAGORA_PROBE_SHARD. Any other repeat is a mistake.
    repeated = [s for s in known if count(==(s), covered) > 1 && s != "09"]
    isempty(repeated) || error("Suites in more than one CI shard: " * join(repeated, ", "))

    if !isempty(probe_shards)
        counts = unique([parse(Int, split(p, "/")[2]) for p in probe_shards])
        length(counts) == 1 || error("Probe shards disagree on the shard count: " * join(probe_shards, ", "))
        n = counts[1]
        indices = sort([parse(Int, split(p, "/")[1]) for p in probe_shards])
        indices == collect(1:n) ||
            error("Probe shards must be 1..$(n) exactly once, got " * join(probe_shards, ", "))
        # The unit tree is one indivisible subprocess, so exactly one probe shard
        # runs it and the rest must opt out.
        skip_unit_count == n - 1 ||
            error("With $(n) probe shards, exactly $(n - 1) must set skip_unit_driver: \"1\"; found $(skip_unit_count)")
    end

    println("ci_test_suite_shard_gate_ok (suites=$(length(known)), shards=$(length(shard_suites)))")
    return nothing
end

ci_test_suite_shard_gate()
