using Test
using SpaceAGORA

# The precompile workload must leave nothing in the package image that
# describes the host it ran on. Its campaign warmups reach machine_topology()
# through the route planners, and a const cache filled during precompilation is
# served, unchanged, to every process that loads the image: the build host's
# cores and memory, and SPACEAGORA_CORE_BUDGET / SPACEAGORA_PHYSICAL_CORES read
# at build time instead of at run time.
#
# Checked in a fresh process, which is the only place the image's state can be
# seen: every module-level cache right after loading, then the whole workload
# run again and reset the way the precompile block resets it, then the same
# caches again. They must match.

const PS_PROJECT = dirname(Base.active_project())

const PS_SNAPSHOT_SCRIPT = raw"""
using SpaceAGORA

# Values that differ between two snapshots by design, not by leftover state:
# the native-lock window's start time (reset to "now"), the per-acquisition
# scratch the tracked lock overwrites on every entry before reading it, and
# write/generation counters that only ever count up.
const PS_VOLATILE = Set(["_NATIVE_LOCK_RESET_NS", "_native_lock_entry_ns",
                         "_native_lock_entry_wait_ns", "_ATOMIC_WRITE_COUNTER",
                         "_MACHINE_CONSTANTS_GENERATION"])

function ps_describe(v)
    v isa Base.RefValue && return isassigned(v) ? "Ref " * repr(v[]; context = :limit => true) : "Ref #undef"
    v isa Threads.Atomic && return "Atomic " * repr(v[])
    (v isa AbstractDict || v isa AbstractSet) && return "$(nameof(typeof(v))) len=$(length(v))"
    return nothing
end

function ps_snapshot()
    out = Dict{String, String}()
    seen = Set{Module}()
    function walk(m::Module)
        m in seen && return
        push!(seen, m)
        for n in names(m; all = true, imported = false)
            isdefined(m, n) || continue
            v = try getfield(m, n) catch; continue end
            if v isa Module
                parentmodule(v) === m && v !== m && walk(v)
                continue
            end
            (isconst(m, n) && !(String(n) in PS_VOLATILE)) || continue
            d = ps_describe(v)
            d === nothing || (out["$(m).$(n)"] = first(d, 400))
        end
    end
    walk(SpaceAGORA)
    return out
end

loaded = ps_snapshot()
topology_at_load = SpaceAGORA.ParallelProfiles._TOPOLOGY_CACHE[]

# The workload, exactly as the precompile block runs it, then its reset.
SpaceAGORA._run_spaceagora_precompile_workload()
SpaceAGORA._warm_predictive_campaign()
SpaceAGORA._warm_mixed_dispatch_campaign()
SpaceAGORA.SimulationCampaigns._warm_campaign_dispatchers()
filled = SpaceAGORA.ParallelProfiles._TOPOLOGY_CACHE[] !== nothing
SpaceAGORA.SimulationModel.Planets._reset_furnished_kernels!()
SpaceAGORA._reset_process_local_state!()
after = ps_snapshot()

println("TOPOLOGY_AT_LOAD=", topology_at_load === nothing ? "nothing" : "filled")
println("WORKLOAD_FILLED_TOPOLOGY=", filled)
for k in sort!(collect(union(keys(loaded), keys(after))))
    a = get(loaded, k, "<absent>")
    b = get(after, k, "<absent>")
    a == b || println("LEFTOVER ", k, " :: ", a, " => ", b)
end
t = SpaceAGORA.ParallelProfiles.machine_topology()
println("USABLE_CORES=", t.usable_cores, " SOURCE=", t.source)
"""

function ps_run(env::Pair...)
    script = tempname() * ".jl"
    write(script, PS_SNAPSHOT_SCRIPT)
    cmd = `$(Base.julia_cmd()) --startup-file=no --threads=2 --project=$(PS_PROJECT) $(script)`
    out = IOBuffer()
    ok = success(pipeline(setenv(addenv(cmd, env...); dir = mktempdir()); stdout = out, stderr = devnull))
    return ok, String(take!(out))
end

@testset "the precompile workload leaves no process or host state behind" begin
    ok, out = ps_run("SPACEAGORA_CORE_BUDGET" => "3",
                     "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0",
                     "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0")
    @test ok
    ok || println(out)
    # A fresh process starts with no topology reading ...
    @test occursin("TOPOLOGY_AT_LOAD=nothing", out)
    # ... the workload does fill it, which is why the reset exists ...
    @test occursin("WORKLOAD_FILLED_TOPOLOGY=true", out)
    # ... and after the reset every cache is as it was at load.
    leftovers = filter(l -> startswith(l, "LEFTOVER"), split(out, '\n'))
    @test isempty(leftovers)
    isempty(leftovers) || foreach(println, leftovers)
    # The run-time override is the one the topology honors.
    @test occursin("USABLE_CORES=3 SOURCE=override", out)
end

@testset "the reset clears the caches a run fills" begin
    PP = SpaceAGORA.ParallelProfiles
    SE = SpaceAGORA.SimulationEngine
    PPol = SpaceAGORA.SimulationModel.ParallelPolicy
    PP.machine_topology()
    SE._calib_machine_label()
    PPol.hint_overhead_ns()
    @test PP._TOPOLOGY_CACHE[] !== nothing
    @test !isempty(SE._CALIB_MACHINE_LABEL[])
    @test PPol._HINT_OVERHEAD_NS[] >= 0.0
    SpaceAGORA._reset_process_local_state!()
    @test PP._TOPOLOGY_CACHE[] === nothing
    @test isempty(SE._CALIB_MACHINE_LABEL[])
    @test PPol._HINT_OVERHEAD_NS[] == -1.0
    # A changed override is read on the next use.
    withenv("SPACEAGORA_CORE_BUDGET" => "5") do
        @test PP.usable_core_budget() == 5
        SpaceAGORA._reset_process_local_state!()
    end
    haskey(ENV, "SPACEAGORA_CORE_BUDGET") || @test PP.machine_topology().source !== :override
end
