using Test
using SpaceAGORA

# `allow_inner_with_outer` is permission, not instruction.
#
# Nested under an active outer split, the shipped rule threaded the inner work
# whenever the flag was set and the item count cleared the threshold, whatever
# the per-sample budget was. On mcgrid_16sat_8mc (8 samples of 16 spacecraft,
# 32 threads) that means 16 satellites split across `fld(32, 8) = 4` threads,
# and the split costs more than it returns: outer_inner_static, which is that
# flag with no adaptive routing at all, loses to outer_threads by nearly the
# same margin R6 does. Under R6 the nested decision is now cost-aware -- it
# declines a split whose per-thread item count it cannot amortize -- while every
# other profile, outer_inner_static included, keeps the shipped behavior so the
# pinned baselines it is scored against do not move.

const PP = SpaceAGORA.SimulationModel.ParallelPolicy

_nit_env(; v2::Bool, budget::Union{Nothing, Int}) = withenv(
    "SPACEAGORA_PARALLEL_POLICY_V2" => v2 ? "1" : nothing,
    "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
    "SPACEAGORA_INNER_THREAD_BUDGET" => budget === nothing ? nothing : string(budget),
) do
    PP.snapshot_policy_decision_env()
end

_nit_decide(items, env; outer_active, allow_with_outer = true, source = :nit_probe) =
    PP.with_policy_context() do
        PP.thread_policy_decision(items; mode = :auto, threshold = 1,
                                  outer_active = outer_active,
                                  allow_with_outer = allow_with_outer,
                                  source = source, env = env)
    end

@testset "nested inner threading is declined when the budget cannot amortize" begin
    if Base.Threads.nthreads() >= 4
        env = _nit_env(v2 = true, budget = 4)
        @test env.policy_v2
        @test env.inner_thread_budget == 4
        # The P5 point: 16 spacecraft on a 4-thread per-sample budget, under an
        # active outer split that has already claimed the machine.
        d = _nit_decide(16, env; outer_active = true)
        @test !d.use_threads
        @test d.allotment == 1
        # The same shape with no outer split above it keeps threading: this
        # gate is about nesting, not about 16 items on 4 threads as such.
        @test _nit_decide(16, env; outer_active = false).use_threads
        # Permission still withdrawn explicitly means no, as before.
        @test !_nit_decide(4, env; outer_active = true, allow_with_outer = false).use_threads
    else
        @test_skip "needs >= 4 Julia threads"
    end
end

@testset "nested inner threading is still allowed when it can amortize" begin
    if Base.Threads.nthreads() >= 4
        budget = Base.Threads.nthreads()
        env = _nit_env(v2 = true, budget = budget)
        @test env.inner_thread_budget == budget
        per_thread = PP.NESTED_INNER_MAX_ITEMS_PER_THREAD
        # At the bound the split is allowed; one item past it, it is not.
        allowed = _nit_decide(per_thread * budget, env; outer_active = true)
        declined = _nit_decide(per_thread * budget + 1, env; outer_active = true)
        @test allowed.use_threads
        @test allowed.allotment == budget
        @test !declined.use_threads
    else
        @test_skip "needs >= 4 Julia threads"
    end
end

@testset "the gate is R6-only: the pinned hybrid baseline does not move" begin
    if Base.Threads.nthreads() >= 4
        # outer_inner_static is R3 with allow_inner_with_outer on, and it is one
        # of the routes R6 is scored against. Its answer must be the shipped one.
        env_r5 = _nit_env(v2 = false, budget = 4)
        @test !env_r5.policy_v2
        @test _nit_decide(16, env_r5; outer_active = true).use_threads
        # And with no snapshot at all (the live-ENV path), likewise unchanged.
        d = withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "4") do
            _nit_decide(16, nothing; outer_active = true)
        end
        @test d.use_threads
    else
        @test_skip "needs >= 4 Julia threads"
    end
end
