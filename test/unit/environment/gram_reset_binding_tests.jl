# Native-free regression of the actual canonical reset and seed-selection bodies.
# A binding is installed inside an already-running caller, matching first native
# wrapper loading. This does not load GRAM, SPICE or SpaceAGORA, or establish
# native/scientific acceptance.
module GramResetBindingTests
using Test

const SOURCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "ext", "SpaceAGORAGRAMSuiteExt.jl"))
const GLOBAL_LOCK = ReentrantLock()
const EVENTS = Any[]

module EM
const LOCK_SCOPE = Ref(:global)
_gram_lock_scope() = LOCK_SCOPE[]
struct GRAMAtmosphereModel
    core::Any
    instance_lock::ReentrantLock
    constructor_kwargs::Union{Nothing, Dict{Symbol, Any}}
end
function reset_density_model_history! end
end

_tl(::Symbol) = GLOBAL_LOCK
function _ensure_native_static_tables!(core)
    push!(EVENTS, (:warm, islocked(GLOBAL_LOCK), islocked(core.gram_atmosphere.model_lock)))
    return nothing
end

function extract_between(source::String, first_marker::String, next_marker::String)
    first_range = findfirst(first_marker, source)
    first_range === nothing && error("Missing source marker: $first_marker")
    next_range = findnext(next_marker, source, last(first_range) + 1)
    next_range === nothing && error("Missing source marker: $next_marker")
    return source[first(first_range):prevind(source, first(next_range))]
end
const SOURCE = read(SOURCE_PATH, String)
for method_source in (
    extract_between(SOURCE, "@inline function _gram_call_lock(", "function __init__()"),
    extract_between(SOURCE, "@inline function _gram_core_seed(", "function EM._rebuild_gram_epoch_model("),
)
    Base.include_string(@__MODULE__, method_source, SOURCE_PATH * " [extracted reset boundary]")
end

mutable struct StubAtmosphere
    seed::Int
    step::Int
    model_lock::ReentrantLock
    expected_lock::ReentrantLock
    other_lock::ReentrantLock
    events::Vector{Any}
end

function model_with(gram; recipe=Dict{Symbol, Any}(:seed => 29), core_recipe=nothing)
    model_lock = ReentrantLock()
    expected_lock, other_lock = EM.LOCK_SCOPE[] === :global ?
        (GLOBAL_LOCK, model_lock) : (model_lock, GLOBAL_LOCK)
    atmosphere = StubAtmosphere(71, 4, model_lock, expected_lock, other_lock, EVENTS)
    core = core_recipe === nothing ? (gram=gram, gram_atmosphere=atmosphere) :
        (gram=gram, gram_atmosphere=atmosphere, _constructor_kwargs=core_recipe)
    return EM.GRAMAtmosphereModel(core, model_lock, recipe)
end

function install_setter!(gram)
    Core.eval(gram, quote
        function set_seed!(atmosphere, seed)
            push!(atmosphere.events,
                (:seed, seed, islocked(atmosphere.expected_lock), islocked(atmosphere.other_lock)))
            atmosphere.seed = seed
            atmosphere.step = 0
            return nothing
        end
    end)
    return nothing
end

# A deterministic stand-in walk distinguishes reseeding from merely returning
# true, and detects accidental draws during reset without a native dependency.
function sequence!(atmosphere)
    return [begin
        atmosphere.step += 1
        (atmosphere.seed, atmosphere.step)
    end for _ in 1:3]
end

function cold_reset_case(scope, recipe, core_recipe, expected_seed)
    EM.LOCK_SCOPE[] = scope
    empty!(EVENTS)
    gram = Module(gensym(:LateGRAM))
    model = model_with(gram; recipe, core_recipe)
    initial_recipe = deepcopy(recipe)
    initial_core_recipe = deepcopy(core_recipe)
    atmosphere = model.core.gram_atmosphere
    advanced = sequence!(atmosphere)
    fresh = StubAtmosphere(expected_seed, 0, ReentrantLock(), ReentrantLock(), ReentrantLock(), Any[])
    reference = sequence!(fresh)
    @test advanced != reference

    # Do not enter invokelatest around this caller or reset: the actual method
    # must observe a binding created after this function's execution began.
    @test !Base.invokelatest(isdefined, gram, :set_seed!)
    install_setter!(gram)
    @test Base.invokelatest(isdefined, gram, :set_seed!)
    @test EM.reset_density_model_history!(model)
    @test atmosphere.seed == expected_seed
    @test atmosphere.step == 0
    @test EVENTS == [(:warm, false, false), (:seed, expected_seed, true, false)]
    @test !islocked(GLOBAL_LOCK) && !islocked(model.instance_lock)
    @test sequence!(atmosphere) == reference
    @test isequal(model.constructor_kwargs, initial_recipe)
    @test core_recipe === nothing || isequal(model.core._constructor_kwargs, initial_core_recipe)

    # A second ordinary reset uses the same semantics after the binding exists.
    empty!(EVENTS)
    @test EM.reset_density_model_history!(model)
    @test atmosphere.step == 0
    @test sequence!(atmosphere) == reference
    @test EVENTS == [(:warm, false, false), (:seed, expected_seed, true, false)]
end

function checked_cold_reset_case(args...)
    # Julia 1.12 can permit stale binding access with a runtime warning under
    # ordinary settings, but reject it with --depwarn=error. Check stderr as
    # well as the result so default CI catches the faulty lookup too.
    mktemp() do path, io
        redirect_stderr(io) do
            cold_reset_case(args...)
        end
        flush(io)
        @test read(path, String) == ""
    end
end

@testset "Canonical GRAM reset sees cold native bindings" begin
    for scope in (:global, :model)
        checked_cold_reset_case(scope, Dict{Symbol, Any}(:seed => Int32(29)), nothing, 29)
        checked_cold_reset_case(scope, Dict{Symbol, Any}(), nothing, 1001)
        checked_cold_reset_case(scope, nothing, Dict{Symbol, Any}(:seed => 43), 43)
        checked_cold_reset_case(scope, nothing, Dict{Symbol, Any}(), 1001)
        # An explicit wrapper recipe takes precedence over the raw-core recipe.
        checked_cold_reset_case(scope, Dict{Symbol, Any}(:seed => 29), Dict{Symbol, Any}(:seed => 43), 29)
    end
end

@testset "Canonical GRAM reset retains unsupported cases and errors" begin
    EM.LOCK_SCOPE[] = :global
    gram = Module(gensym(:NoSeedBinding))
    for model in (
        model_with(gram),
        model_with(gram; recipe=nothing),
        EM.GRAMAtmosphereModel((gram=gram,), ReentrantLock(), Dict{Symbol, Any}(:seed => 7)),
        EM.GRAMAtmosphereModel((gram_atmosphere=nothing,), ReentrantLock(), Dict{Symbol, Any}(:seed => 7)),
        EM.GRAMAtmosphereModel((gram=gram, gram_atmosphere=nothing, _constructor_kwargs=nothing),
            ReentrantLock(), nothing),
    )
        empty!(EVENTS)
        @test !EM.reset_density_model_history!(model)
        @test isempty(EVENTS)
    end

    # A present setter that fails must still propagate its error and release
    # the selected mutex. This is separate from an absent binding returning false.
    for scope in (:global, :model)
        EM.LOCK_SCOPE[] = scope
        empty!(EVENTS)
        throwing = Module(gensym(:ThrowingSeedBinding))
        model = model_with(throwing)
        Core.eval(throwing, quote
            function set_seed!(atmosphere, seed)
                push!(atmosphere.events,
                    (:seed_error, seed, islocked(atmosphere.expected_lock), islocked(atmosphere.other_lock)))
                error("fixture seed failure")
            end
        end)
        @test_throws ErrorException EM.reset_density_model_history!(model)
        @test EVENTS == [(:warm, false, false), (:seed_error, 29, true, false)]
        @test !islocked(GLOBAL_LOCK) && !islocked(model.instance_lock)
        @test model.core.gram_atmosphere.seed == 71
        @test model.core.gram_atmosphere.step == 4
    end
end
end
