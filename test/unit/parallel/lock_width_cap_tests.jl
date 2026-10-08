using Test
using SpaceAGORA

const RSC = SpaceAGORA.RuntimeServices
const SEC = SpaceAGORA.SimulationEngine

# Evaluate the actual runtime definitions with a controlled clock and snapshot.
# The loaded package keeps its real clock, counters and lock-accounting methods.
module LockWidthClockFixture
const CLOCK = Ref{UInt64}(0)
const _NATIVE_LOCK_RESET_NS = Ref{UInt64}(0)
const SNAPSHOT = Ref((hold_ns = Int64(0), wait_hold_ratio = 0.0, acquisitions = Int64(0)))
time_ns() = CLOCK[]
native_lock_stats_snapshot() = SNAPSHOT[]

const SOURCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "simulation", "runtime_services.jl"))
const SOURCE_AST = Meta.parseall(read(SOURCE_PATH, String); filename = SOURCE_PATH)
const RUNTIME_MODULE = only(filter(SOURCE_AST.args) do node
    node isa Expr && node.head === :module && node.args[2] === :RuntimeServices
end)

function definition_name(node)
    if node isa Expr && node.head === :macrocall && node.args[1] === GlobalRef(Core, Symbol("@doc"))
        node = node.args[end]
    end
    node isa Expr && node.head === :function || return nothing
    signature = node.args[1]
    if signature isa Expr && signature.head === :(::)
        signature = signature.args[1]
    end
    signature isa Expr && signature.head === :call || return nothing
    return (name = signature.args[1], definition = node)
end

for name in (:native_lock_occupancy, :lock_width_ceiling)
    definitions = filter(!isnothing, definition_name.(RUNTIME_MODULE.args[3].args))
    selected = only(filter(def -> def.name === name, definitions))
    Core.eval(@__MODULE__, selected.definition)
end

function sample!(hold_ns, window_ns; wait_hold = 0.25, acquisitions = 7)
    # A nonzero origin also checks that occupancy uses time since the reset.
    _NATIVE_LOCK_RESET_NS[] = 100
    CLOCK[] = 100 + window_ns
    SNAPSHOT[] = (hold_ns = Int64(hold_ns), wait_hold_ratio = wait_hold, acquisitions = Int64(acquisitions))
    return nothing
end
end

const LWC = LockWidthClockFixture

@testset "native-lock width cap" begin
    @testset "ceiling from occupancy" begin
        RSC.reset_native_lock_stats!()
        # No lock activity means no constraint -- not a ceiling of one.
        @test RSC.lock_width_ceiling() == typemax(Int)
        o = RSC.native_lock_occupancy()
        @test o.rho == 0.0
        @test o.window_s >= 0.0
    end

    @testset "controlled occupancy and advancing window" begin
        LWC.sample!(300, 800)
        occ = LWC.native_lock_occupancy()
        @test occ == (rho = 0.375, wait_hold = 0.25, acquisitions = 7, window_s = 800 / 1.0e9)
        @test LWC.lock_width_ceiling() == 3
        # A later call samples a later window. The old test compared the later
        # ceiling to the earlier occupancy and raced with coverage overhead.
        LWC.CLOCK[] += 200
        @test LWC.native_lock_occupancy().rho == 0.3
        @test LWC.lock_width_ceiling() == 4
        @test occ.rho == 0.375

        LWC.sample!(0, 1000; wait_hold = 0.0, acquisitions = 0)
        @test LWC.native_lock_occupancy().rho == 0.0
        @test LWC.lock_width_ceiling() == typemax(Int)
        LWC.sample!(300, 0)
        @test LWC.native_lock_occupancy().rho == 0.0
        @test LWC.lock_width_ceiling() == typemax(Int)
        LWC.sample!(1200, 1000)
        @test LWC.native_lock_occupancy().rho == 1.0
        @test LWC.lock_width_ceiling() == 1
    end

    @testset "worker scaling and minimum denominator" begin
        LWC.sample!(300, 800)
        @test LWC.native_lock_occupancy(workers = 2).rho == 0.1875
        @test LWC.lock_width_ceiling(workers = 2) == 6
        for workers in (0, -2)
            @test LWC.native_lock_occupancy(workers = workers).rho == 0.375
            @test LWC.lock_width_ceiling(workers = workers) == 3
        end
    end

    @testset "floor_rho refuses to build a ceiling out of noise" begin
        LWC.sample!(20, 1000)
        @test LWC.lock_width_ceiling() == 50  # Inclusive default floor.
        LWC.sample!(19, 1000)
        @test LWC.lock_width_ceiling() == typemax(Int)
        LWC.sample!(300, 800)
        @test LWC.lock_width_ceiling(floor_rho = 0.375) == 3
        @test LWC.lock_width_ceiling(floor_rho = 0.4) == typemax(Int)
    end

    @testset "the clamp only removes width, never adds it" begin
        params = (; shared_buffers = (; rhs_width_ceiling = Ref{Int}(3)))
        wide = SEC._make_calib_flat_plan(12, :static)
        @test SEC._clamp_plan_to_lock_ceiling(wide, params).allotment == 3
        # Already under the ceiling: untouched, and the same object's other
        # fields survive.
        narrow = SEC._make_calib_flat_plan(2, :static)
        clamped = SEC._clamp_plan_to_lock_ceiling(narrow, params)
        @test clamped.allotment == 2
        @test clamped.scheduler === narrow.scheduler
        @test clamped.mode === narrow.mode
        # Unconstrained ceiling changes nothing.
        open_params = (; shared_buffers = (; rhs_width_ceiling = Ref{Int}(typemax(Int))))
        @test SEC._clamp_plan_to_lock_ceiling(wide, open_params).allotment == 12
        # No params, no clamp -- callers outside a solve are unaffected.
        @test SEC._clamp_plan_to_lock_ceiling(wide, nothing).allotment == 12
    end

    @testset "satellite_batch is exempt" begin
        # It takes its width from Polyester's own pool and honours neither
        # `allotment` nor the inner thread budget, so clamping the field would
        # change what the plan reports without changing what it runs.
        params = (; shared_buffers = (; rhs_width_ceiling = Ref{Int}(1)))
        batch = SEC._make_calib_satellite_batch_plan()
        @test SEC._clamp_plan_to_lock_ceiling(batch, params) === batch
    end
end
