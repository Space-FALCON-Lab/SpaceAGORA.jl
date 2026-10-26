# Reports: plain values every check returns, readable without Test.

"""
    ConformanceCheck(id, status, since, detail)

One check of a [`ConformanceReport`](@ref).

- `id::Symbol`: the check, for example `:repeatable_run`.
- `status::Symbol`: `:pass`, `:fail`, `:skip` (the check could not apply, `detail`
  says why) or `:advisory_fail` (a check added in the current minor release that
  did not pass; it becomes required in the next one).
- `since::VersionNumber`: the SpaceAGORA release that added the check.
- `detail::String`: what differed, or why the check was skipped; empty on a pass.

A passing check certifies interface behaviour only. Scientific validation of the
model stays with its author.
"""
struct ConformanceCheck
    id::Symbol
    status::Symbol
    since::VersionNumber
    detail::String
    function ConformanceCheck(id::Symbol, status::Symbol, since::VersionNumber, detail::AbstractString)
        status in CHECK_STATUSES || throw(ArgumentError(
            "ConformanceCheck status must be one of $(join(repr.(CHECK_STATUSES), ", ")); got $(repr(status))."))
        return new(id, status, since, String(detail))
    end
end

const CHECK_STATUSES = (:pass, :fail, :skip, :advisory_fail)

"""
    ConformanceReport(target, kind, checks)

The result of one `check_*` call: `target` names the object checked, `kind` the
contract it was checked against (for example `:force_torque_model`), and `checks`
holds one [`ConformanceCheck`](@ref) per check id, in a fixed order.
[`passed`](@ref) tells whether every required check passed; `show` prints one line
per check.
"""
struct ConformanceReport
    target::String
    kind::Symbol
    checks::Vector{ConformanceCheck}
end

"""
    passed(report::ConformanceReport) -> Bool

`true` when no check of `report` has status `:fail`. Skipped checks and advisory
failures do not fail a report.
"""
passed(report::ConformanceReport) = !any(check -> check.status === :fail, report.checks)

"""
    passed(check::ConformanceCheck) -> Bool

`true` unless the check's status is `:fail`.
"""
passed(check::ConformanceCheck) = check.status !== :fail

# The first release that ships the suite: its checks are required from it.
const FIRST_RELEASE = v"0.3.0"

# K-4: a check added after the first release is advisory while the package is in
# the minor release that added it (or an earlier one), and required afterwards.
# Outside a package (the standalone contract tests) the version is the first release.
function _failure_status(since::VersionNumber,
        current::VersionNumber=something(pkgversion(parentmodule(@__MODULE__)), FIRST_RELEASE))
    since > FIRST_RELEASE && (since.major, since.minor) >= (current.major, current.minor) && return :advisory_fail
    return :fail
end

_check(id::Symbol, since::VersionNumber, ok::Bool, detail::AbstractString="") =
    ConformanceCheck(id, ok ? :pass : _failure_status(since), since, ok ? "" : detail)
_skip(id::Symbol, since::VersionNumber, why::AbstractString) = ConformanceCheck(id, :skip, since, why)

# An error's message for a check's detail: at most six lines.
function _summary(err)
    err isa TaskFailedException && (err = err.task.exception)
    text = sprint(showerror, err)
    lines = split(text, '\n')
    return length(lines) > 6 ? join(lines[1:6], '\n') * "\n…" : text
end

# `f()`, with the global random number generator restored afterwards: a defect a
# check provokes (an algorithm or family that calls `rand()`) does not leak into the
# caller's later draws.
function _preserving_global_rng(f)
    rng = Random.default_rng()
    saved = copy(rng)
    try
        return f()
    finally
        copy!(rng, saved)
    end
end

# `f()` with the global random number generator seeded `seed`, restored afterwards.
function _with_global_seed(f, seed)
    rng = Random.default_rng()
    saved = copy(rng)
    try
        Random.seed!(seed)
        return f()
    finally
        copy!(rng, saved)
    end
end

function Base.show(io::IO, ::MIME"text/plain", report::ConformanceReport)
    verdict = passed(report) ? "passed" : "FAILED"
    println(io, "Conformance of ", report.target, " (", report.kind, "): ", verdict)
    for check in report.checks
        print(io, "  ", rpad(String(check.id), 30), " ", check.status)
        isempty(check.detail) || print(io, ": ", replace(check.detail, "\n" => "\n      "))
        println(io)
    end
    print(io, "Interface conformance only; scientific validation stays with the model's author.")
end

Base.show(io::IO, report::ConformanceReport) =
    print(io, "ConformanceReport(", repr(report.target), ", :", report.kind, ", ",
        count(passed, report.checks), "/", length(report.checks), " checks not failed)")

# test_report also writes each report into this directory when it is set; the
# contributor CI workflow reads the files for its summary (design section 5, step 4).
const REPORT_DIRECTORY_VARIABLE = "SPACEAGORA_CONFORMANCE_REPORT_DIR"
const _REPORTS_WRITTEN = Threads.Atomic{Int}(0)

# One TOML file per report, named so that a directory listing sorts in the order the
# reports were written: `target`, `kind`, and one `[[checks]]` table per check with
# `id`, `status`, `since` and `detail`.
function _write_report_file(report::ConformanceReport, directory::AbstractString=get(ENV, REPORT_DIRECTORY_VARIABLE, ""))
    isempty(directory) && return nothing
    mkpath(directory)
    n = Threads.atomic_add!(_REPORTS_WRITTEN, 1) + 1
    name = replace(report.target, r"[^A-Za-z0-9_.-]" => "_")
    path = joinpath(directory, string(lpad(n, 4, '0'), "-", getpid(), "-", name, "-", report.kind, ".toml"))
    checks = [Dict("id" => String(check.id), "status" => String(check.status), "since" => string(check.since),
        "detail" => check.detail) for check in report.checks]
    data = Dict("target" => report.target, "kind" => String(report.kind), "checks" => checks)
    open(io -> TOML.print(io, data; sorted=true), path, "w")
    return path
end

"""
    test_report(report::ConformanceReport)

Record `report` in the current `Test` testset: a testset named after the target,
with one testset per check. A failed check fails its test and prints its detail; a
skipped check is recorded as skipped and an advisory failure as broken. Returns
`report`.

When the environment variable `SPACEAGORA_CONFORMANCE_REPORT_DIR` names a directory,
the report is also written there as a TOML file (`target`, `kind`, and one
`[[checks]]` table per check with `id`, `status`, `since` and `detail`), which the
contributor CI workflow summarizes.

Needs the `Test` standard library: `using Test` loads the method.
"""
function test_report end
