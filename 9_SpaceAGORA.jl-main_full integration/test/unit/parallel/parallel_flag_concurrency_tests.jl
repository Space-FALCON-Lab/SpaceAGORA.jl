using Test
using SpaceAGORA

# SolverConfig(parallel=true) scopes its profile on the process-global ENV.
# Two flagged runs started from different tasks do not nest, so a per-call
# save/restore let the second run's "previous" values (the first run's profile)
# outlive both. These tests drive the overlay with runs that finish in the
# opposite order they started, with and without an exception, and through the
# SimulationEngineConfig path, and check that ENV ends exactly as it began.

const CF_SE = SpaceAGORA.SimulationEngine
const CF_PR = SpaceAGORA.ParallelProfiles
const CF_PC = SpaceAGORA.SimulationModel.ParallelCost

const CF_TMP = mktempdir()
const CF_CONSTANTS = joinpath(CF_TMP, "cost_constants.toml")
const CF_FLAG_PAIRS = CF_PR.parallel_flag_env_pairs()
const CF_ENGINE_CONFIG = CF_SE.SimulationEngineConfig()
const CF_KEYS = sort!(unique!(vcat(
    first.(CF_FLAG_PAIRS),
    collect(keys(CF_SE._engine_env_overrides(CF_ENGINE_CONFIG; parallel_flag = true))),
)))

cf_snapshot() = Dict(k => get(ENV, k, nothing) for k in CF_KEYS)
cf_flag_view() = Dict(k => get(ENV, k, nothing) for k in first.(CF_FLAG_PAIRS))
const CF_FLAG_ENV = Dict{String, Union{Nothing, String}}(k => v for (k, v) in CF_FLAG_PAIRS)

# A flagged scope that signals when it is inside and waits to be released.
function cf_start(entry::Function; fail::Bool = false)
    inside = Channel{Nothing}(1)
    release = Channel{Nothing}(1)
    task = Threads.@spawn entry() do
        put!(inside, nothing)
        take!(release)
        fail && error("concurrent flagged run failed")
        return cf_flag_view()
    end
    take!(inside)
    return (task = task, release = release)
end

cf_flag_entry(f) = CF_SE._with_parallel_flag(f, true)
cf_engine_entry(f) = CF_SE._with_engine_env_overrides(CF_ENGINE_CONFIG, f; parallel_flag = true)

withenv("SPACEAGORA_COST_CONSTANTS_PATH" => CF_CONSTANTS,
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
        "SPACEAGORA_INNER_THREAD_BUDGET" => nothing,
        "SPACEAGORA_CAMPAIGN_PLANNER" => nothing,
        "SPACEAGORA_PARALLEL_POLICY_V2" => nothing,
        "SPACEAGORA_PARALLEL_PROFILE" => nothing) do
    # Settle the constants path without measuring anything: a failed
    # calibration is logged once and not retried in this process.
    @test_logs (:warn,) (:info,) match_mode = :any CF_PC.ensure_machine_constants!(
        path = CF_CONSTANTS, calibrate = () -> error("no calibration in this test"))

    @testset "overlapping flagged runs that finish out of order leave ENV unchanged" begin
        before = cf_snapshot()
        a = cf_start(cf_flag_entry)
        b = cf_start(cf_flag_entry)
        @test CF_SE._parallel_flag_env_depth() == 2
        put!(a.release, nothing)
        @test fetch(a.task) == CF_FLAG_ENV
        # B is still running: its profile must still be in place.
        @test cf_flag_view() == CF_FLAG_ENV
        put!(b.release, nothing)
        @test fetch(b.task) == CF_FLAG_ENV
        @test CF_SE._parallel_flag_env_depth() == 0
        @test cf_snapshot() == before
    end

    @testset "a flagged run that throws while another is active" begin
        before = cf_snapshot()
        a = cf_start(cf_flag_entry; fail = true)
        b = cf_start(cf_flag_entry)
        put!(a.release, nothing)
        @test_throws TaskFailedException fetch(a.task)
        @test cf_flag_view() == CF_FLAG_ENV
        put!(b.release, nothing)
        fetch(b.task)
        @test CF_SE._parallel_flag_env_depth() == 0
        @test cf_snapshot() == before

        # The last one out throwing still restores.
        a = cf_start(cf_flag_entry)
        b = cf_start(cf_flag_entry; fail = true)
        put!(a.release, nothing)
        fetch(a.task)
        put!(b.release, nothing)
        @test_throws TaskFailedException fetch(b.task)
        @test CF_SE._parallel_flag_env_depth() == 0
        @test cf_snapshot() == before
    end

    @testset "the engine-config path and the plain path overlap" begin
        before = cf_snapshot()
        a = cf_start(cf_engine_entry)
        b = cf_start(cf_flag_entry)
        put!(a.release, nothing)
        @test fetch(a.task) == CF_FLAG_ENV
        @test cf_flag_view() == CF_FLAG_ENV
        put!(b.release, nothing)
        fetch(b.task)
        @test cf_snapshot() == before
        @test CF_SE._engine_active_overrides_ref[] === nothing

        a = cf_start(cf_flag_entry)
        b = cf_start(cf_engine_entry)
        put!(a.release, nothing)
        fetch(a.task)
        @test cf_flag_view() == CF_FLAG_ENV
        put!(b.release, nothing)
        @test fetch(b.task) == CF_FLAG_ENV
        @test CF_SE._parallel_flag_env_depth() == 0
        @test cf_snapshot() == before
        # Neither order leaves an engine-config override set active.
        @test CF_SE._engine_active_overrides_ref[] === nothing
        @test CF_SE._engine_active_config_ref[] === nothing
    end

    @testset "an env_overrides entry that differs from the flag still wins, and is restored" begin
        before = cf_snapshot()
        cfg = CF_SE.SimulationEngineConfig(env_overrides = Dict("SPACEAGORA_CAMPAIGN_PLANNER" => "bandit"))
        seen = CF_SE._with_engine_env_overrides(cfg, () -> get(ENV, "SPACEAGORA_CAMPAIGN_PLANNER", nothing);
                                                parallel_flag = true)
        @test seen == "bandit"
        @test cf_snapshot() == before
    end
end
