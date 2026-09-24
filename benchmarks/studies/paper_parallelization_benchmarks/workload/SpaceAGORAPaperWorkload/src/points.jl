# The precompile workload's point list, derived from the paper harness's own
# definitions rather than copied from them.
#
# Needs, already defined in the including scope: PAPER_BENCHMARK_PHASES
# (paper_parallelization_benchmarks/cli.jl), ppc_case_catalog and
# PPC_GRAM_LIVE_CASES (parallelization_performance/cases.jl), and
# ppc_mode_specs (parallelization_performance/modes.jl). The workload package
# includes those files before this one; the coverage gate
# (test/gates/ci_ppb_workload_coverage_gate.jl) includes them into a sandbox
# module, so the list can be checked without building the image.
#
# A case added to a P phase is picked up here automatically. The only way a
# P-phase case can be absent from the list is the documented exclusion below,
# and the gate fails if anything else goes missing.

"""
    PPBWorkloadPoint

One workload entry: a harness case under one mode, run the way a harness worker
would run it (`kind = :perf` for a timed point, `:parity` for a trajectory-parity
point). `samples` is the number of Monte Carlo samples the workload itself runs:
two for a campaign case (enough to reach the campaign runner and the batch
dispatch), one otherwise. The sample count does not enter any type, so the
specializations are the ones the real point needs.
"""
struct PPBWorkloadPoint
    phase::String
    case::String
    mode::String
    kind::Symbol
    samples::Int
end

# The paper figures' phases: P1-P6 and P6p.
ppb_workload_phase_ids() = [p.id for p in PAPER_BENCHMARK_PHASES if occursin(r"^P[0-9]", p.id)]

ppb_workload_phases() = [p for p in PAPER_BENCHMARK_PHASES if p.id in ppb_workload_phase_ids()]

"""
    ppb_workload_excluded(case) -> Bool

True for cases that build a live native GRAM atmosphere, which the workload
leaves out; the set is `PPC_GRAM_LIVE_CASES`, the same list the harness uses to
decide when to load GRAMSuite.

Warming them at precompile time would need GRAMSuite as a dependency of the
workload package, which makes the image unloadable on a machine without the GRAM
build. It would also run native GRAM inside the precompile process: GRAMSuite
`include`s its native wrapper module on first construction, which is an
evaluation into a closed module during incremental compilation, and any model
object reachable from the workload module's globals would be serialized into the
image with its native handles. Harness points for these cases are launched
without the workload.
"""
ppb_workload_excluded(case::AbstractString) = any(c -> occursin(c, case), PPC_GRAM_LIVE_CASES)

function ppb_workload_points(phases=ppb_workload_phases(); catalog=ppc_case_catalog())
    points = PPBWorkloadPoint[]
    seen = Set{Tuple{String, String, Symbol}}()
    add!(phase, case, mode, kind) = begin
        key = (case, mode, kind)
        key in seen && return nothing
        push!(seen, key)
        spec = catalog[case]
        push!(points, PPBWorkloadPoint(phase, case, mode, kind, spec.montecarlo ? 2 : 1))
        return nothing
    end
    for phase in phases
        for case in phase.cases, mode in phase.modes
            ppb_workload_excluded(case) || add!(phase.id, case, mode, :perf)
        end
        # The controller runs parity for every non-serial mode of the phase.
        for case in phase.parity_cases, mode in phase.modes
            mode == "serial" && continue
            ppb_workload_excluded(case) || add!(phase.id, case, mode, :parity)
        end
    end
    return points
end

# (phase, case) pairs the P phases name but the workload leaves out.
function ppb_workload_skipped(phases=ppb_workload_phases())
    skipped = Tuple{String, String}[]
    for phase in phases, case in unique(vcat(phase.cases, phase.parity_cases))
        ppb_workload_excluded(case) && push!(skipped, (phase.id, case))
    end
    return skipped
end
