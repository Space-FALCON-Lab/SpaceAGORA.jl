# What the measured per-spacecraft footprint changes in outer-route selection.
#
# Prints, for two machine profiles and three constellation sizes, the memory
# charge per process worker, `effective_process_workers`, whether the process
# route is offered at all (`>= 2` affordable workers), and the mixed-dispatch
# local slots beside those workers -- once with the shipped-before constant
# (90 MB/spacecraft, reproduced through `SPACEAGORA_GRAM_SAT_MEMORY_MB=90`) and
# once with the measured defaults.
#
#   julia --project=. --threads=8 \
#       benchmarks/studies/gram_memory_footprint/route_fit_table.jl
#
# Machine facts are supplied through the existing override hooks, so this runs
# on any host and answers for the profile, not for the host it runs on:
#
#   SPACEAGORA_CORE_BUDGET          usable cores
#   SPACEAGORA_MEMORY_BUDGET_GB     memory a run may spend across its processes
#   SPACEAGORA_MEMORY_AVAILABLE_GB  headroom the kernel would hand out
#   SPACEAGORA_PERF_WORKER_MEMORY_GB  the package + SPICE + GRAM image per worker
#
# The coordinator's own thread count is whatever this process was started with
# and is deliberately NOT part of the comparison: the local-slot term
# `min(nthreads - 1, usable_cores - W)` is unchanged by this work, so both the
# before and the after column see the same value.

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
import Pkg
if something(Base.active_project(), "") != joinpath(REPO_ROOT, "Project.toml")
    Pkg.activate(REPO_ROOT; io=devnull)
end

using Printf
using SpaceAGORA
const PP = SpaceAGORA.ParallelProfiles

# Machine facts. `total_gb` and `cores` name the machine; `budget_gb` is that
# total less the 10% OS reserve `machine_topology` applies, and `headroom_gb` is
# the budget less the coordinator's own resident set, which is what
# `memory_worker_cap` has left to spend. Pinning the headroom (rather than
# reading this host's `MemAvailable`) is what makes the table a statement about
# the profile instead of about the box the script happens to run on.
struct Profile
    name::String
    cores::Int
    total_gb::Float64
    worker_image_gb::Float64
end

budget_gb(p::Profile) = p.total_gb * 0.90
headroom_gb(p::Profile) = budget_gb(p) - p.worker_image_gb

const PROFILES = [
    # 60 GB / 12 physical cores: this workstation (space-falcon-1).
    Profile("workstation-60GB-12c", 12, 60.0, 1.5),
    # 250 GB / 64 physical cores: the TRX50-class profile.
    Profile("trx50-250GB-64c", 64, 250.0, 1.5),
]

const SIZES = [256, 1024, 4096]

function profile_env(p::Profile, gram_mb::Union{Nothing, Real})
    pairs = Pair{String, Union{String, Nothing}}[
        "SPACEAGORA_PARALLEL_POLICY_V2" => "1",
        "SPACEAGORA_CORE_BUDGET" => string(p.cores),
        "SPACEAGORA_MEMORY_BUDGET_GB" => string(round(budget_gb(p); digits=3)),
        "SPACEAGORA_MEMORY_AVAILABLE_GB" => string(round(headroom_gb(p); digits=3)),
        "SPACEAGORA_PERF_WORKER_MEMORY_GB" => string(p.worker_image_gb),
        "SPACEAGORA_GRAM_SAT_MEMORY_MB" => gram_mb === nothing ? nothing : string(gram_mb),
    ]
    return pairs
end

function row(p::Profile, n::Int, gram_mb::Union{Nothing, Real})
    withenv(profile_env(p, gram_mb)...) do
        PP.refresh_machine_topology!()
        f = PP.OuterRouteFeatures(
            category="montecarlo", n_sats=n, density_family="gram_point",
            mission_time_s=600.0, montecarlo_samples=64
        )
        t = PP.OuterRouteTuning(memory_aware=true, mixed_dispatch=true)
        extra = PP.native_gram_worker_extra_bytes(n)
        workers = PP.effective_process_workers(f, t)
        slots = PP.mixed_local_slots(f, t, workers)
        return (extra_gb=extra / (1 << 30), workers=workers,
                fits=workers >= 2, slots=slots, capacity=workers + slots)
    end
end

function main()
    println("coordinator threads = $(Threads.nthreads())")
    println()
    @printf("%-22s %6s  %-6s %10s %8s %6s %7s %9s\n",
            "profile", "N", "charge", "extra_GB", "workers", "fits", "slots", "capacity")
    for p in PROFILES, n in SIZES
        before = row(p, n, 90)       # the shipped-before constant
        after = row(p, n, nothing)   # the measured defaults
        for (label, r) in (("before", before), ("after", after))
            @printf("%-22s %6d  %-6s %10.2f %8d %6s %7d %9d\n",
                    p.name, n, label, r.extra_gb, r.workers, r.fits ? "yes" : "NO",
                    r.slots, r.capacity)
        end
    end
    # Leave the cached snapshot describing the real host again.
    PP.refresh_machine_topology!()
    return nothing
end

main()
