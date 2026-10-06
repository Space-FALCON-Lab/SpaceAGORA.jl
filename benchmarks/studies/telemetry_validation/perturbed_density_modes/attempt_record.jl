# One run tag's outputs belong to one attempt (run_mode.jl).
#
# begin_attempt! removes every output an earlier attempt may have left in the tag
# directory before the new attempt runs, so a failed rerun cannot leave older
# trajectory files beside its own summary. finish_attempt! records each output
# file's SHA-256 in the summary and publishes the summary last, by rename, so a
# summary on disk always names exactly the files of the attempt it describes.
# analyze.py compares runs only when both summaries say complete and the
# recorded identities match the files.
#
# No SpaceAGORA dependency, so the rules are testable without a run.
module PerturbedModeAttempt

using Dates
using SHA
using TOML

export ATTEMPT_OUTPUTS, begin_attempt!, attempt_status, finish_attempt!, file_sha256

const ATTEMPT_OUTPUTS = ("run_summary.toml", "extrema.csv", "simulation_results.csv", "perturbation_log.csv")
const SOLVER_OK = ("Success", "Terminated")

"Remove the tag's previous outputs and return a new attempt identity."
function begin_attempt!(dir::AbstractString)::String
    mkpath(dir)
    for f in ATTEMPT_OUTPUTS
        rm(joinpath(dir, f); force=true)
    end
    return string(Dates.format(now(UTC), "yyyymmddTHHMMSS.sss"), "Z-pid", getpid(), "-", string(time_ns(); base=16))
end

file_sha256(path::AbstractString)::String = open(io -> bytes2hex(sha256(io)), path)

"""
    attempt_status(; error, retcode, have_trajectory, n_apo, n_peri, orbits,
                   completed_orbits, termination_cause) -> (status, reason)

`complete`, `incomplete` (a successful solve short of its orbit count) or `failed`.
The runtime orbit counter must reach the request, with a consistent termination
cause. Saved extrema alone cannot prove completion: the final apoapsis may be
absent from saved samples, and an early stop may have the same extrema counts.
Once completion is established, saved coverage needs `orbits - 1` of each apsis.
"""
function attempt_status(; error::AbstractString, retcode::AbstractString, have_trajectory::Bool,
                        n_apo::Integer, n_peri::Integer, orbits::Integer,
                        completed_orbits=nothing, termination_cause=nothing)
    isempty(strip(error)) || return ("failed", "the solve threw")
    retcode in SOLVER_OK || return ("failed", "retcode $(retcode)")
    have_trajectory || return ("failed", "no trajectory was written")
    (!(orbits isa Bool) && orbits >= 1) || return ("failed", "invalid requested orbit count")
    (completed_orbits isa Integer && !(completed_orbits isa Bool) && completed_orbits >= 0) ||
        return ("failed", "missing or invalid completed orbit count")
    expected_cause = retcode == "Success" ? "end_of_time_span" :
                     completed_orbits >= orbits ? "orbit_count" : "terminated_before_orbit_count"
    (termination_cause isa AbstractString && termination_cause == expected_cause) ||
        return ("failed", "missing or inconsistent termination cause")
    completed_orbits >= orbits ||
        return ("incomplete", "$(completed_orbits) completed orbit events for $(orbits) requested")
    (n_apo >= orbits - 1 && n_peri >= orbits - 1) ||
        return ("incomplete", "$(n_apo) apoapses and $(n_peri) periapses for $(orbits) orbits")
    return ("complete", "")
end

"Record the identity of every output present and publish the summary."
function finish_attempt!(dir::AbstractString, summary::AbstractDict)
    hashes = Dict{String, Any}()
    for f in ATTEMPT_OUTPUTS
        f == "run_summary.toml" && continue
        path = joinpath(dir, f)
        isfile(path) && (hashes[f] = file_sha256(path))
    end
    summary["output_sha256"] = hashes
    partial = joinpath(dir, ".run_summary.toml.partial")
    open(io -> TOML.print(io, summary), partial, "w")
    mv(partial, joinpath(dir, "run_summary.toml"); force=true)
    return summary
end

end
