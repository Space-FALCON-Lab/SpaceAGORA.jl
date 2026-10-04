# Outcome of one dispersed-campaign member, from its summary row (run_campaign.jl).
#
#   dispatch_failure   the dispatcher caught an exception and holds no member value
#   solver_failure     the member's own solve threw (odyssey_sample.jl catches it and
#                      returns retcode "ERROR" with the error text) or returned a
#                      retcode other than Success or Terminated
#   early_termination  a successful solve that stopped before the analysis horizon:
#                      fewer than horizon_passes + 1 apoapses or horizon_passes passes
#                      (Terminated covers both reaching the orbit count and impact)
#   complete           everything else; only these enter fixed-horizon statistics
#
# The classification is by coverage of the horizon, not by termination cause: a
# member that impacts after its horizon's last apoapsis still has every pass of
# the horizon. The cause itself is recorded per member (termination_cause,
# odyssey_sample.jl).
#
# No SpaceAGORA dependency, so the rules are testable without a campaign.
module DispersedMemberStatus

export member_status, status_counts, member_execution_totals, MEMBER_STATUSES

const MEMBER_STATUSES = ("complete", "early_termination", "solver_failure", "dispatch_failure")
const SOLVER_OK = ("Success", "Terminated")

_field(row, k::Symbol, default) = hasproperty(row, k) ? getproperty(row, k) : default

function member_status(row; horizon_passes::Int)::String
    horizon_passes >= 1 || throw(ArgumentError("horizon_passes must be >= 1, got $horizon_passes"))
    _field(row, :dispatch_success, false) === true || return "dispatch_failure"
    err = _field(row, :error, "")
    (err isa AbstractString && !isempty(strip(err))) && return "solver_failure"
    string(_field(row, :retcode, "ERROR")) in SOLVER_OK || return "solver_failure"
    n_apo = _field(row, :n_apo, 0); n_peri = _field(row, :n_peri, 0)
    (n_apo >= horizon_passes + 1 && n_peri >= horizon_passes) || return "early_termination"
    return "complete"
end

"Counts per status, every status present (zero when absent)."
function status_counts(rows; horizon_passes::Int)::Dict{String, Int}
    counts = Dict(s => 0 for s in MEMBER_STATUSES)
    for r in rows
        counts[member_status(r; horizon_passes=horizon_passes)] += 1
    end
    return counts
end

"Execution totals from returned members; entirely failed dispatches contribute zero."
function member_execution_totals(rows)
    solve_s = 0.0
    pids = Set{Int}()
    for row in rows
        _field(row, :dispatch_success, false) === true || continue
        elapsed = _field(row, :solve_s, missing)
        ismissing(elapsed) || (solve_s += elapsed)
        pid = _field(row, :pid, missing)
        ismissing(pid) || push!(pids, pid)
    end
    return (sum_member_solve_s=solve_s, distinct_member_pids=length(pids))
end

end
