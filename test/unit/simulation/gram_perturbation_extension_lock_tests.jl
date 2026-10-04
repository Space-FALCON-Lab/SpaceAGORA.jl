# Native-free extension regression. This evaluates the exact four
# selected method bodies from the extension against stand-ins. It does
# not load SpaceAGORA, GRAMSuite, SPICE or any native library, and does not test
# native CSPICE safety or the callback's first_update flag lifecycle.
module GramPerturbationExtensionLockTests
using Test

const SOURCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "ext", "SpaceAGORAGRAMSuiteExt.jl"))
const EVENTS = Any[]
const LOCK_STACK = Any[]
const GLOBAL_MUTEX = ReentrantLock()

struct TrackedLock
    site::Symbol
    mutex::ReentrantLock
end

module RS
using ..GramPerturbationExtensionLockTests: TrackedLock, GLOBAL_MUTEX
tracked_lock(site::Symbol) = TrackedLock(site, GLOBAL_MUTEX)
end

module EM
const LOCK_SCOPE = Ref(:model)
_gram_lock_scope() = LOCK_SCOPE[]
struct GRAMAtmosphereModel
    core::Any
    instance_lock::ReentrantLock
end
function _gram_walk_sample end
function _gram_walk_reseed! end
end

lock_label(x) = x isa TrackedLock ? x.site : :model
current_lock() = isempty(LOCK_STACK) ? :unlocked : lock_label(last(LOCK_STACK))

mutable struct StubAtmosphere
    seed::Int
    native_calls::Int
    draws_since_seed::Int
    state::Any
end

function stub_density_state(a::StubAtmosphere)
    @assert current_lock() !== :unlocked
    push!(EVENTS, (kind=:read, lock=current_lock()))
    return a.state
end

function stub_reseed!(a::StubAtmosphere, seed::Int)
    @assert current_lock() !== :unlocked
    push!(EVENTS, (kind=:reseed, lock=current_lock(), seed=seed))
    a.seed = seed
    a.draws_since_seed = 0
    return nothing
end

module GRAMSuite
using ..GramPerturbationExtensionLockTests: TrackedLock, EVENTS, LOCK_STACK, current_lock, lock_label
function _with_gram_lock(f, lock_obj)
    mutex = lock_obj isa TrackedLock ? lock_obj.mutex : lock_obj
    lock(mutex) do
        push!(EVENTS, (kind=:enter, lock=lock_label(lock_obj)))
        push!(LOCK_STACK, lock_obj)
        try
            return f()
        finally
            pop!(LOCK_STACK)
            push!(EVENTS, (kind=:exit, lock=lock_label(lock_obj)))
        end
    end
end
function _gram_density_state_native(core, h, lat, lon, t, wind)
    @assert current_lock() !== :unlocked
    a = core.gram_atmosphere
    a.native_calls += 1
    a.draws_since_seed += 1
    push!(EVENTS, (kind=:query, lock=current_lock(), h=h, lat=lat, lon=lon, t=t, wind=wind))
    a.state = (perturbedDensity=2.0 + 0.01 * a.seed + 0.1 * a.draws_since_seed,
               density=2.0, densityStandardDeviation=0.05, relativeStepSize=0.25)
    return nothing
end
end

function extract_between(source::String, first_marker::String, next_marker::String)
    first_range = findfirst(first_marker, source)
    first_range === nothing && error("Missing source marker: $first_marker")
    next_range = findnext(next_marker, source, last(first_range) + 1)
    next_range === nothing && error("Missing source marker: $next_marker")
    return source[first(first_range):prevind(source, first(next_range))]
end

const SOURCE = read(SOURCE_PATH, String)
const EXTRACTED = [
    extract_between(SOURCE, "@inline function _gram_call_lock(", "function __init__()"),
    extract_between(SOURCE, "@inline function _gram_density_state_tuple(", "function EM._gram_walk_clone("),
    extract_between(SOURCE, "function EM._gram_walk_sample(", "function EM._gram_walk_reseed!("),
    extract_between(SOURCE, "function EM._gram_walk_reseed!(", "EM._gram_recipe_seed("),
]
_tl(site::Symbol) = RS.tracked_lock(site)
for method_source in EXTRACTED
    Base.include_string(@__MODULE__, method_source, SOURCE_PATH * " [extracted test boundary]")
end

function new_model(seed=11)
    atmos = StubAtmosphere(seed, 0, 0, nothing)
    gram = NamedTuple{(:get_density_state, Symbol("set_seed!"))}((stub_density_state, stub_reseed!))
    return EM.GRAMAtmosphereModel((gram=gram, gram_atmosphere=atmos), ReentrantLock())
end

function checked_sample(model; first_update=nothing)
    before_calls = model.core.gram_atmosphere.native_calls
    before_events = length(EVENTS)
    value = if first_update === nothing
        EM._gram_walk_sample(model, -100.0, 0.2, -0.3, 12.0)
    else
        EM._gram_walk_sample(model, -100.0, 0.2, -0.3, 12.0; first_update=first_update)
    end
    events = EVENTS[(before_events + 1):end]
    @test model.core.gram_atmosphere.native_calls == before_calls + 1
    @test getproperty.(events, :kind) == [:enter, :query, :read, :exit]
    @test all(e -> e.lock == events[1].lock, events)
    @test events[2].h == -30.0
    @test (events[2].lat, events[2].lon, events[2].t, events[2].wind) == (0.2, -0.3, 12.0, false)
    @test isempty(LOCK_STACK)
    @test value isa NTuple{4, Float64}
    return value, events[1].lock
end

@testset "Extracted production GRAM walk method: lock and query boundary" begin
    println("source=", SOURCE_PATH)
    println("scope=exact extracted methods; stand-in locking/native calls; no package or native load")
    for scope in (:model, :global)
        EM.LOCK_SCOPE[] = scope
        empty!(EVENTS)
        m = new_model()
        normal_lock = scope === :model ? :model : :gram_density
        _, first_lock = checked_sample(m; first_update=true)
        @test first_lock === :gram_setup
        _, later_lock = checked_sample(m; first_update=false)
        @test later_lock === normal_lock
        calls_before = m.core.gram_atmosphere.native_calls
        EM._gram_walk_reseed!(m, 29)
        @test m.core.gram_atmosphere.native_calls == calls_before
        @test EVENTS[end - 1] == (kind=:reseed, lock=normal_lock, seed=29)
        reseeded, reseed_lock = checked_sample(m; first_update=true)
        @test reseed_lock === :gram_setup
        @test m.core.gram_atmosphere.native_calls == 3
        @test count(e -> e.kind === :query, EVENTS) == 3
        @test count(e -> e.kind === :read, EVENTS) == 3
        @test count(e -> e.kind === :reseed, EVENTS) == 1
        _, default_lock = checked_sample(m)
        @test default_lock === normal_lock

        # Changing the lock flag cannot add a draw or change the returned sample.
        a, b = new_model(29), new_model(29)
        value_global, _ = checked_sample(a; first_update=true)
        value_normal, _ = checked_sample(b; first_update=false)
        @test value_global == value_normal == reseeded
        @test a.core.gram_atmosphere.native_calls == b.core.gram_atmosphere.native_calls == 1
        println("lock_scope=", scope, " first=", first_lock, " ordinary=", later_lock,
                " reseeded=", reseed_lock, " default=", default_lock,
                " sequence_queries=3 sequence_reseeds=1 draw_parity=true")
    end
end
end
