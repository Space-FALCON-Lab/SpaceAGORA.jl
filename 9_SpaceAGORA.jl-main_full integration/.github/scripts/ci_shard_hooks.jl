# CI sharding hooks for the default test entrypoint and the example smoke.
#
# Julia CI splits the test suite and the example smoke across parallel jobs.
# `.github/scripts/ci_shard_plan.py` assigns every test item to one shard and
# passes the shard its list in SPACEAGORA_CI_SHARD_ITEMS. This module runs
# exactly those items and writes a report of what it ran, which
# `.github/scripts/ci_shard_verify.py` checks against the whole universe so a
# test that stops running fails CI instead of going quiet.
#
# The items are:
#   suite:NN     a numbered test/suites file, raw-included as usual
#   probe:FILE   one entry of the probe list in test/suites/09_probe_drivers.jl
#   unit:PATH    one include of test/unit/runtests.jl
#
# Suite 09 and the unit runner are not edited to support this. Suite 09 is
# included with a `mapexpr` that narrows its `probe_files` list to this shard's
# probes and points its unit driver at `ci_unit_shard.jl`, which includes only
# this shard's unit files. The rewrite is strict: if either file's structure
# drifts from what it expects, it errors instead of guessing.
#
# Without SPACEAGORA_CI_SHARD_ITEMS none of this is loaded and the entrypoint
# runs every suite in one process, as it always has.
module CIShard

using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SUITES_DIR = joinpath(REPO_ROOT, "test", "suites")
const PROBE_SUITE = "09_probe_drivers.jl"
const PROBE_TESTSET = "Probe Drivers"
const UNIT_TESTSET = "Standalone Unit Suite Driver"
const UNIT_RUNNER = joinpath(REPO_ROOT, "test", "unit", "runtests.jl")
const UNIT_SHARD_RUNNER = joinpath(@__DIR__, "ci_unit_shard.jl")
const EXAMPLE_SMOKE = joinpath(REPO_ROOT, "test", "smoke", "ci_examples_suite_smoke.jl")
# Items the workflow runs as steps of its own; the Julia side only accepts them.
const WORKFLOW_ITEMS = ("grid", "depwarn")

const TIMINGS = Dict{String, Float64}()
const ATTEMPTED = String[]
const FAILED = String[]

function requested_items()
    raw = get(ENV, "SPACEAGORA_CI_SHARD_ITEMS", "")
    items = String.(split(raw, [' ', ',', '\n', '\t']; keepempty=false))
    length(unique(items)) == length(items) || error("SPACEAGORA_CI_SHARD_ITEMS repeats an item: $(raw)")
    return items
end

function report_dir()
    dir = get(ENV, "SPACEAGORA_CI_SHARD_REPORT_DIR", "")
    isempty(dir) && error("SPACEAGORA_CI_SHARD_REPORT_DIR must be set when SPACEAGORA_CI_SHARD_ITEMS is.")
    mkpath(dir)
    return dir
end

function record!(item::String, seconds::Float64; failed::Bool=false)
    TIMINGS[item] = round(seconds; digits=1)
    push!(ATTEMPTED, item)
    failed && push!(FAILED, item)
    println("ci_shard_item ", item, " seconds=", TIMINGS[item], failed ? " FAILED" : "")
    flush(stdout)
    return nothing
end

function write_report(name::String, universe::Vector{String})
    path = joinpath(report_dir(), name * ".toml")
    open(path, "w") do io
        TOML.print(io, Dict(
            "universe" => sort(universe),
            "requested" => sort(requested_items()),
            "attempted" => sort(ATTEMPTED),
            "failed" => sort(FAILED),
            "timings" => TIMINGS,
        ))
    end
    println("ci_shard_report ", path)
    return path
end

# --- AST helpers -------------------------------------------------------------

parse_file(path) = Meta.parseall(read(path, String); filename=path)

function toplevel_exprs(path)
    return [ex for ex in parse_file(path).args if !(ex isa LineNumberNode)]
end

function string_leaves(ex, out=String[])
    if ex isa String
        push!(out, ex)
    elseif ex isa Expr
        for arg in ex.args
            string_leaves(arg, out)
        end
    end
    return out
end

function testset_name(ex)
    (ex isa Expr && ex.head == :macrocall && ex.args[1] == Symbol("@testset")) || return nothing
    for arg in ex.args[2:end]
        arg isa String && return arg
    end
    return nothing
end

is_include_call(ex) = ex isa Expr && ex.head == :call && ex.args[1] == :include && length(ex.args) == 2

# `include(joinpath(@__DIR__, "parallel", "x_tests.jl"))` -> "parallel/x_tests.jl".
# The planner derives the same label from the source text.
function include_label(ex)
    label = join(string_leaves(ex.args[2]), "/")
    isempty(label) && error("Cannot label unit include without a path literal: $(ex)")
    return label
end

# Every `name = rhs` assignment anywhere inside `ex`.
function find_assignments(ex, name::Symbol, out=Expr[])
    if ex isa Expr
        if ex.head == :(=) && ex.args[1] == name
            push!(out, ex)
        end
        for arg in ex.args
            find_assignments(arg, name, out)
        end
    end
    return out
end

function replace_single_assignment(ex, name::Symbol, rhs, where::String)
    ex = deepcopy(ex)
    hits = find_assignments(ex, name)
    length(hits) == 1 || error("Expected exactly one `$(name) = ...` in $(where), found $(length(hits)). " *
                               "Update .github/scripts/ci_shard_hooks.jl to match the file.")
    hits[1].args[2] = rhs
    return ex
end

# --- Universe ------------------------------------------------------------------

function probe_universe()
    path = joinpath(SUITES_DIR, PROBE_SUITE)
    hits = Expr[]
    for ex in toplevel_exprs(path)
        append!(hits, find_assignments(ex, :probe_files))
    end
    length(hits) == 1 || error("Expected one `probe_files = [...]` in $(PROBE_SUITE), found $(length(hits)).")
    rhs = hits[1].args[2]
    (rhs isa Expr && rhs.head == :vect && all(a -> a isa String, rhs.args)) ||
        error("`probe_files` in $(PROBE_SUITE) is no longer a literal list of file names.")
    probes = String.(rhs.args)
    length(unique(probes)) == length(probes) || error("`probe_files` in $(PROBE_SUITE) repeats a probe.")
    return probes
end

function unit_universe()
    labels = [include_label(ex) for ex in toplevel_exprs(UNIT_RUNNER) if is_include_call(ex)]
    length(unique(labels)) == length(labels) || error("test/unit/runtests.jl includes a file twice.")
    return labels
end

function check_suite_list(suites::Vector{String})
    on_disk = sort([f for f in readdir(SUITES_DIR) if occursin(r"^\d\d_.*\.jl$", f)])
    sort(suites) == on_disk || error(
        "test/integration/runtests.jl lists suites $(sort(suites)) but test/suites holds $(on_disk). " *
        "A suite missing from the list would never run.")
    PROBE_SUITE in suites || error("$(PROBE_SUITE) is missing from the suite list.")
    return nothing
end

suite_item(file) = "suite:" * file[1:2]

function test_universe(suites::Vector{String})
    return vcat([suite_item(s) for s in suites if s != PROBE_SUITE],
                "probe:" .* probe_universe(),
                "unit:" .* unit_universe())
end

function check_requested(universe::Vector{String}, requested::Vector{String})
    unknown = [i for i in requested if !(i in universe) && !(i in WORKFLOW_ITEMS)]
    isempty(unknown) || error("This shard was asked for items that do not exist: $(join(unknown, ", ")). " *
                              "Known items: $(join(universe, ", "))")
    return nothing
end

# --- Suite 09 rewrite -----------------------------------------------------------

# Iterates like the probe list it replaces and times each probe: a probe's time
# runs from the loop handing it out to the loop asking for the next one, and it
# only counts as attempted once the loop has come back for more.
mutable struct TimedProbeList
    probes::Vector{String}
    last::Float64
end

Base.length(t::TimedProbeList) = length(t.probes)
Base.eltype(::Type{TimedProbeList}) = String
function Base.iterate(t::TimedProbeList, i::Int=1)
    now = time()
    i > 1 && record!("probe:" * t.probes[i - 1], now - t.last)
    i > length(t.probes) && return nothing
    t.last = now
    return (t.probes[i], i + 1)
end

function suite09_mapexpr(probes::Vector{String}, run_units::Bool, seen::Set{String})
    return function (ex)
        name = testset_name(ex)
        if name == PROBE_TESTSET
            push!(seen, name)
            isempty(probes) && return nothing
            return replace_single_assignment(ex, :probe_files, TimedProbeList(probes, time()), PROBE_SUITE)
        elseif name == UNIT_TESTSET
            push!(seen, name)
            run_units || return nothing
            return replace_single_assignment(ex, :unit_script, UNIT_SHARD_RUNNER, PROBE_SUITE)
        end
        error("$(PROBE_SUITE) has a top-level expression the CI shard hook does not know how to split: " *
              "$(ex isa Expr ? ex.head : typeof(ex)) $(something(name, "")). " *
              "Update .github/scripts/ci_shard_hooks.jl so it runs in exactly one shard.")
    end
end

# --- Entry points ---------------------------------------------------------------

function run_integration(suites::Vector{String})
    check_suite_list(suites)
    requested = requested_items()
    universe = test_universe(suites)
    check_requested(universe, requested)
    probes = [p for p in probe_universe() if ("probe:" * p) in requested]
    units = [u for u in unit_universe() if ("unit:" * u) in requested]

    start = tryparse(Float64, get(ENV, "SPACEAGORA_CI_STEP_START", ""))
    start === nothing || (TIMINGS["preamble"] = round(time() - start; digits=1))
    println("ci_shard suites=", join([s for s in suites if suite_item(s) in requested], ","),
            " probes=", length(probes), " units=", length(units))

    for suite in suites
        path = joinpath(SUITES_DIR, suite)
        if suite == PROBE_SUITE
            isempty(probes) && isempty(units) && continue
            seen = Set{String}()
            t0 = time()
            try
                Base.include(suite09_mapexpr(probes, !isempty(units), seen), Main, path)
                seen == Set([PROBE_TESTSET, UNIT_TESTSET]) ||
                    error("$(PROBE_SUITE) no longer defines the testsets \"$(PROBE_TESTSET)\" and \"$(UNIT_TESTSET)\".")
            catch err
                push!(FAILED, "suite:09")
                @error "Suite $(suite) failed in this shard" exception=(err, catch_backtrace())
            end
            TIMINGS["suite:09"] = round(time() - t0; digits=1)
        else
            item = suite_item(suite)
            item in requested || continue
            t0 = time()
            failed = false
            try
                Base.include(Main, path)
            catch err
                failed = true
                @error "Suite $(suite) failed" exception=(err, catch_backtrace())
            end
            record!(item, time() - t0; failed)
        end
    end

    write_report("integration", universe)
    isempty(FAILED) || error("CI shard failures: $(join(FAILED, ", "))")
    return nothing
end

# Runs in the unit-suite subprocess that suite 09 starts.
function run_unit_shard()
    requested = requested_items()
    universe = "unit:" .* unit_universe()
    wanted = Set(i for i in requested if startswith(i, "unit:"))
    isempty(setdiff(wanted, universe)) || error("Unknown unit items: $(join(setdiff(wanted, universe), ", "))")

    mapexpr = function (ex)
        is_include_call(ex) || return ex
        item = "unit:" * include_label(ex)
        item in wanted || return nothing
        return :($(timed_include)($item, () -> $ex))
    end
    Base.include(mapexpr, Main, UNIT_RUNNER)

    write_report("unit", universe)
    isempty(FAILED) || error("Unit shard failures: $(join(FAILED, ", "))")
    return nothing
end

function timed_include(item::String, f)
    t0 = time()
    failed = false
    try
        f()
    catch err
        failed = true
        @error "Unit file $(item) failed" exception=(err, catch_backtrace())
    end
    record!(item, time() - t0; failed)
    return nothing
end

# Runs the example smoke on this shard's examples only.
function run_example_shard()
    requested = requested_items()
    universe = String[]
    replaced = Ref(0)
    select = function (files)
        append!(universe, basename.(files))
        unknown = setdiff(requested, universe)
        isempty(unknown) || error("Unknown examples requested: $(join(unknown, ", "))")
        chosen = [f for f in files if basename(f) in requested]
        append!(ATTEMPTED, basename.(chosen))
        return chosen
    end
    mapexpr = function (ex)
        if ex isa Expr && ex.head == :(=) && ex.args[1] == :examples &&
           ex.args[2] == :(list_examples())
            replaced[] += 1
            return :(examples = $(select)(list_examples()))
        end
        return ex
    end
    try
        Base.include(mapexpr, Main, EXAMPLE_SMOKE)
    catch
        append!(FAILED, ATTEMPTED)
        rethrow()
    finally
        replaced[] == 1 || error("Expected one `examples = list_examples()` in the example smoke, found $(replaced[]).")
        write_report("examples", universe)
    end
    return nothing
end

end # module CIShard
