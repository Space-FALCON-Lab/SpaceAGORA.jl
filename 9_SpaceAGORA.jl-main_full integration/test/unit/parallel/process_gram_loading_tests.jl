# Exercise the production bootstrap and shared path helper without importing
# SpaceAGORA/GRAMSuite, spawning processes, installing packages, or using SPICE.
module ProcessGramLoadingTests
using Test
using Logging

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const WORKER_SOURCE = joinpath(ROOT, "src", "orchestration", "parallel", "process", "worker_pool.jl")
const SERVICES_SOURCE = joinpath(ROOT, "src", "orchestration", "simulation", "runtime_services.jl")
const PROJECT = joinpath("synthetic", "repo", "Project.toml")
const VENDOR = joinpath(dirname(PROJECT), "data", "GRAMSuite.jl")
const INITIAL_PATH = ["@", "synthetic-fallback", "@stdlib"]
const HOST_PATH = copy(LOAD_PATH)
const HOST_PROJECT = Base.active_project()
const HOST_PACKAGES = Set(keys(Base.loaded_modules))

Base.@kwdef mutable struct Fixture
    discoverable::Bool = false
    vendor_present::Bool = true
    load_path::Vector{String} = copy(INITIAL_PATH)
    package_failure::Any = nothing
    preload_failure::Any = nothing
    discovery_failure::Any = nothing
    import_failure::Any = nothing
    events::Vector{Any} = Any[]
    worker_main::Module = Module(gensym(:WorkerMain))
    package::Module = Module(gensym(:WorkerSpaceAGORA))
    services::Module = Module(gensym(:WorkerRuntimeServices))
    gram::Module = Module(gensym(:WorkerGRAMSuite))
end

function find_package(fixture, name)
    @assert name == "GRAMSuite"
    push!(fixture.events, (:find_package, name))
    fixture.discovery_failure === nothing || throw(fixture.discovery_failure)
    return fixture.discoverable ? "synthetic-installed-GRAMSuite" : nothing
end
function vendor_isdir(fixture, path)
    @assert path == VENDOR
    push!(fixture.events, (:isdir, path))
    return fixture.vendor_present
end

function remote_eval(fixture, target, workers, expr)
    @assert target === fixture.worker_main
    @assert workers == [77]
    if expr == :(using SpaceAGORA)
        push!(fixture.events, (:spaceagora,))
        fixture.package_failure === nothing || throw(fixture.package_failure)
        if !isdefined(target, :SpaceAGORA)
            Core.eval(target, :(const SpaceAGORA = $(fixture.package)))
        end
    elseif expr == :(import GRAMSuite)
        push!(fixture.events, (:import, copy(fixture.load_path)))
        fixture.import_failure === nothing || throw(fixture.import_failure)
        if !isdefined(target, :GRAMSuite)
            Core.eval(target, :(const GRAMSuite = $(fixture.gram)))
        end
    else
        error("Unexpected remote expression: $expr")
    end
    return nothing
end
function remote_wait(body, fixture, worker, path)
    # A real remote request enters the receiving worker in its current world,
    # including lookup of the binding installed by the preceding remote eval.
    return Base.invokelatest() do
        @assert worker == 77 && path == VENDOR
        @assert isdefined(fixture.worker_main, :SpaceAGORA)
        push!(fixture.events, (:remote_path, path))
        body(path)
    end
end
function preload(fixture, worker)
    @assert worker == 77
    push!(fixture.events, (:preload,))
    fixture.preload_failure === nothing || throw(fixture.preload_failure)
    return nothing
end
function furnish(fixture, worker)
    @assert worker == 77
    push!(fixture.events, (:kernels,))
    return nothing
end

struct FixtureLogger <: AbstractLogger
    fixture::Fixture
end
Logging.min_enabled_level(::FixtureLogger) = Logging.Warn
Logging.shouldlog(::FixtureLogger, level, mod, group, id) = true
Logging.catch_exceptions(::FixtureLogger) = false
function Logging.handle_message(logger::FixtureLogger, level, message, mod, group, id, file, line; kwargs...)
    @assert level == Logging.Warn
    push!(logger.fixture.events, (:warning, String(message), kwargs[:exception]))
end

function definition_name(expr)
    expr isa Expr && expr.head === :function || return nothing
    sig = expr.args[1]
    sig isa Expr && sig.head === :(::) && (sig = sig.args[1])
    return sig isa Expr && sig.head === :call ? sig.args[1] : nothing
end
function definition(path, name; owner=nothing)
    expressions = Meta.parseall(read(path, String); filename=path).args
    if owner !== nothing
        modules = filter(expr -> expr isa Expr && expr.head === :module && expr.args[2] === owner, expressions)
        length(modules) == 1 || error("Expected one $owner module in $path")
        expressions = only(modules).args[end].args
    end
    matches = filter(expr -> definition_name(expr) === name, expressions)
    length(matches) == 1 || error("Expected one $name in $path")
    return only(matches)
end

const SEAMS = Dict(
    :(Base.find_package) => :find_package,
    :isdir => :vendor_isdir,
    :(Distributed.remotecall_eval) => :remote_eval,
    :remotecall_wait => :remote_wait,
    :_preload_worker_packages! => :preload,
    :_furnish_default_spice_kernels! => :furnish,
)
const HELPER_CALL = :(Main.SpaceAGORA.RuntimeServices._prepare_gramsuite_load_path!)
const SAFE_CALLS = Set{Any}((
    :_bootstrap_process_worker!, :_prepare_gramsuite_load_path!,
    :dirname, :joinpath, :catch_backtrace, :(===), :pushfirst!, :push!,
))
function redirect(expr, counts)
    # Main in the transported closure denotes this fixture's worker namespace;
    # keep the full SpaceAGORA.RuntimeServices lookup and real helper call.
    expr === :Main && return :_WORKER_MAIN
    expr isa Expr || return expr
    if expr.head === :quote
        expr.args[1] in (:(using SpaceAGORA), :(import GRAMSuite)) ||
            error("Unexpected quoted bootstrap expression")
        return expr
    elseif expr.head in (:using, :import)
        error("Unredirected package loading")
    elseif expr.head === :macrocall
        expr.args[1] === Symbol("@warn") || error("Unexpected bootstrap macro")
    elseif expr.head === :call
        callee = expr.args[1]
        if haskey(SEAMS, callee)
            seam = SEAMS[callee]
            counts[seam] = get(counts, seam, 0) + 1
            return Expr(:call, GlobalRef(@__MODULE__, seam), :_FIXTURE,
                        (redirect(arg, counts) for arg in expr.args[2:end])...)
        elseif callee == HELPER_CALL
            counts[:helper] = get(counts, :helper, 0) + 1
        else
            callee in SAFE_CALLS || error("Unexpected bootstrap call: $callee")
        end
    end
    return Expr(expr.head, (redirect(arg, counts) for arg in expr.args)...)
end

function caller(fixture=Fixture())
    target = Module(gensym(:BootstrapCaller))
    for mod in (target, fixture.services)
        Core.eval(mod, :(const _FIXTURE = $fixture))
        Core.eval(mod, :(const _WORKER_MAIN = $(fixture.worker_main)))
    end
    Core.eval(fixture.services, :(const LOAD_PATH = $(fixture.load_path)))
    Core.eval(fixture.package, :(const RuntimeServices = $(fixture.services)))
    counts = Dict{Symbol,Int}()
    helper = redirect(definition(SERVICES_SOURCE, :_prepare_gramsuite_load_path!; owner=:RuntimeServices), counts)
    counts == Dict(:find_package => 1, :vendor_isdir => 1) || error("Shared helper boundary changed: $counts")
    Core.eval(fixture.services, helper)
    empty!(counts)
    bootstrap = redirect(definition(WORKER_SOURCE, :_bootstrap_process_worker!), counts)
    counts == Dict(:remote_eval => 2, :remote_wait => 1, :preload => 1, :furnish => 1, :helper => 1) ||
        error("Worker bootstrap boundary changed: $counts")
    Core.eval(target, bootstrap)
    return target, fixture
end
bootstrap(target, fixture) = with_logger(FixtureLogger(fixture)) do
    Base.invokelatest(() -> getfield(target, :_bootstrap_process_worker!)(77, PROJECT))
end
kinds(fixture) = first.(fixture.events)
function capture(fn)
    try
        fn()
        return nothing
    catch err
        return err
    end
end

@testset "worker GRAM bootstrap uses shared bookkeeping" begin
    for discoverable in (false, true)
        target, fixture = caller(Fixture(discoverable=discoverable, vendor_present=!discoverable))
        @test bootstrap(target, fixture) === nothing
        expected_path = discoverable ? INITIAL_PATH : [VENDOR; INITIAL_PATH]
        @test fixture.load_path == expected_path
        expected_events = discoverable ?
            [:spaceagora, :preload, :remote_path, :find_package, :import, :kernels] :
            [:spaceagora, :preload, :remote_path, :find_package, :isdir, :import, :kernels]
        @test kinds(fixture) == expected_events
        @test only(filter(event -> event[1] === :import, fixture.events))[2] == expected_path
        @test Base.invokelatest(getfield, fixture.worker_main, :GRAMSuite) === fixture.gram
        @test !isdefined(target, :GRAMSuite)
        @test !isdefined(fixture.package, :GRAMSuite)
        @test !isdefined(fixture.services, :GRAMSuite)
    end

    # A successful import does not make the synthetic environment discoverable:
    # subsequent attempts must retain the old duplicate-path behavior.
    target, fixture = caller()
    @test bootstrap(target, fixture) === nothing
    @test bootstrap(target, fixture) === nothing
    @test fixture.load_path == [VENDOR; VENDOR; INITIAL_PATH]
    @test count(==(:import), kinds(fixture)) == 2
end

@testset "optional worker GRAM failures warn then furnish kernels" begin
    for vendor_present in (false, true)
        failure = ArgumentError("synthetic missing GRAM dependency")
        target, fixture = caller(Fixture(vendor_present=vendor_present, import_failure=failure))
        @test bootstrap(target, fixture) === nothing
        @test kinds(fixture) == [:spaceagora, :preload, :remote_path, :find_package, :isdir, :import, :warning, :kernels]
        @test fixture.load_path == (vendor_present ? [VENDOR; INITIAL_PATH] : INITIAL_PATH)
        warning = only(filter(event -> event[1] === :warning, fixture.events))
        @test warning[2] == "Process worker 77 could not load GRAMSuite; campaigns using a GRAM density model will fail on this worker."
        @test warning[3][1] === failure
        @test !isempty(warning[3][2])
        @test !isdefined(fixture.worker_main, :GRAMSuite)
        @test count(==(:import), kinds(fixture)) == 1
    end
    failure = ErrorException("synthetic discovery failure")
    target, fixture = caller(Fixture(discovery_failure=failure))
    @test bootstrap(target, fixture) === nothing
    @test kinds(fixture) == [:spaceagora, :preload, :remote_path, :find_package, :warning, :kernels]
    @test fixture.load_path == INITIAL_PATH
    @test only(filter(event -> event[1] === :warning, fixture.events))[3][1] === failure
end

@testset "mandatory worker setup failures still propagate" begin
    for stage in (:package_failure, :preload_failure)
        fixture = Fixture()
        failure = ErrorException("synthetic mandatory setup failure")
        setfield!(fixture, stage, failure)
        target, fixture = caller(fixture)
        @test capture(() -> bootstrap(target, fixture)) === failure
        @test kinds(fixture) == (stage === :package_failure ? [:spaceagora] : [:spaceagora, :preload])
        @test fixture.load_path == INITIAL_PATH
        @test !isdefined(fixture.worker_main, :GRAMSuite)
    end
end

@testset "worker fixture stays inside its native-free boundary" begin
    @test_throws ErrorException redirect(:(run(`false`)), Dict{Symbol,Int}())
    @test_throws ErrorException redirect(:(import GRAMSuite), Dict{Symbol,Int}())
    @test_throws ErrorException redirect(:(@eval import GRAMSuite), Dict{Symbol,Int}())
    @test_throws ErrorException redirect(:(Main.SpaceAGORA.RuntimeServices.unexpected!()), Dict{Symbol,Int}())
    @test LOAD_PATH == HOST_PATH
    @test Base.active_project() == HOST_PROJECT
    forbidden = Set(("SpaceAGORA", "GRAMSuite", "SPICE", "Pkg", "Distributed"))
    @test isempty(filter(id -> id.name in forbidden, setdiff(Set(keys(Base.loaded_modules)), HOST_PACKAGES)))
end
end # module ProcessGramLoadingTests
