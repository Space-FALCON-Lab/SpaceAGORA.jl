# Execute selected production loader definitions against synthetic package and
# environment boundaries. Never include the full example/benchmark entrypoints.
module GramLoadingContractTests
using Test

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const EXAMPLES = joinpath(ROOT, "examples", "common.jl")
const RUNTIME = joinpath(ROOT, "benchmarks", "studies", "performance_runtime_analysis", "main.jl")
const PACKAGE = joinpath(ROOT, "benchmarks", "studies", "parallelization_performance", "cases.jl")
const REAL_PATH = copy(Base.LOAD_PATH)
const REAL_PROJECT = Base.active_project()
const REAL_DIRECTORY = pwd()
const REAL_ARGS = copy(ARGS)
const REAL_PACKAGES = Set(keys(Base.loaded_modules))
const REPO_SENTINEL = joinpath("synthetic", "SpaceAGORA")
const VENDOR_SENTINEL = joinpath(REPO_SENTINEL, "data", "GRAMSuite.jl")
const INITIAL_PATH = ["@", joinpath("synthetic", "other-project"), "@stdlib"]

Base.@kwdef mutable struct Environment
    discoverable::Bool = false
    vendor_present::Bool = true
    available_after_import::Bool = true
    load_path::Vector{String} = copy(INITIAL_PATH)
    project::Union{Nothing, String} = joinpath("synthetic", "caller", "Project.toml")
    outcomes::Vector{Any} = Any[nothing]
    instantiate_error::Any = nothing
    events::Vector{Any} = Any[]
    gram::Module = Module(gensym(:SyntheticGRAMSuite))
end

function find_package(env, name)
    @assert name == "GRAMSuite"
    push!(env.events, (:find_package, name))
    return env.discoverable ? joinpath("synthetic", "available", "GRAMSuite.jl") : nothing
end
function vendor_isdir(env, path)
    @assert path == VENDOR_SENTINEL
    push!(env.events, (:isdir, path))
    return env.vendor_present
end
function active_project(env)
    push!(env.events, (:active_project, env.project))
    return env.project
end
function require_pkg(env, id)
    @assert id == Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg")
    push!(env.events, (:require_pkg, id.name))
    return (
        activate=function(path; io)
            @assert io === devnull
            push!(env.events, (:activate, path))
            env.project = joinpath(path, "Project.toml")
            nothing
        end,
        instantiate=function(; io)
            @assert io === devnull
            push!(env.events, (:instantiate, env.project))
            env.instantiate_error === nothing || throw(env.instantiate_error)
            nothing
        end,
    )
end
function import_gram!(env, caller)
    push!(env.events, (:import, caller, copy(env.load_path)))
    outcome = isempty(env.outcomes) ? nothing : popfirst!(env.outcomes)
    outcome === nothing || throw(outcome)
    if isdefined(caller, :GRAMSuite)
        @assert getfield(caller, :GRAMSuite) === env.gram
    else
        Core.eval(caller, :(const GRAMSuite = $(env.gram)))
    end
    env.discoverable |= env.available_after_import
    return nothing
end

# Restrict redirection to known external operations. The control flow and real
# path/list/error/alias operations remain in the parsed production definitions.
const SEAMS = Dict(
    :(Base.find_package) => :find_package,
    :(Base.active_project) => :active_project,
    :(Base.require) => :require_pkg,
    :isdir => :vendor_isdir,
)
const SAFE_CALLS = Set{Any}((
    :something, :!, :dirname, :isempty, :sprint, :occursin, :ErrorException,
    :pushfirst!, :push!, :throw, :*, :(===), :isdefined, :joinpath,
    :(Base.UUID), :(Base.PkgId), :(Pkg.activate), :(Pkg.instantiate), :(Core.eval),
    :_instantiate_vendored_gramsuite!, :_missing_dependency_error,
    :ensure_gramsuite_loaded!, :setup_gram_example!,
    :_ensure_perf_gramsuite_loaded!, :ppc_ensure_gramsuite_loaded!,
))
const IMPORT_EXPR = :(import GRAMSuite)
const EXPECTED_SEAMS = Dict(
    :_instantiate_vendored_gramsuite! => Dict(:require_pkg => 1, :active_project => 1),
    :_missing_dependency_error => Dict{Symbol, Int}(),
    :ensure_gramsuite_loaded! => Dict(:find_package => 1, :vendor_isdir => 2, :import_gram! => 2),
    :setup_gram_example! => Dict{Symbol, Int}(),
    :_ensure_perf_gramsuite_loaded! => Dict(:find_package => 1, :vendor_isdir => 1, :import_gram! => 1),
    :ppc_ensure_gramsuite_loaded! => Dict(:find_package => 1, :vendor_isdir => 1, :import_gram! => 1),
)
function redirect(ex, counts)
    ex isa Expr || return ex
    if ex.head === :macrocall
        if ex.args[1] === Symbol("@eval") && ex.args[end] == IMPORT_EXPR
            counts[:import_gram!] = get(counts, :import_gram!, 0) + 1
            return Expr(:call, GlobalRef(@__MODULE__, :import_gram!), :_GRAM_FIXTURE, :(@__MODULE__))
        end
        ex.args[1] === Symbol("@__MODULE__") || error("Unexpected loader macro: $(ex.args[1])")
    elseif ex.head in (:import, :using)
        error("Unredirected package loading in fixture")
    elseif ex.head === :call
        callee = ex.args[1]
        if haskey(SEAMS, callee)
            seam = SEAMS[callee]
            counts[seam] = get(counts, seam, 0) + 1
            return Expr(:call, GlobalRef(@__MODULE__, seam), :_GRAM_FIXTURE,
                        (redirect(x, counts) for x in ex.args[2:end])...)
        end
        callee in SAFE_CALLS || error("Unexpected loader call: $callee")
    end
    return Expr(ex.head, (redirect(x, counts) for x in ex.args)...)
end
function definition_name(ex)
    ex isa Expr && ex.head === :function || return nothing
    sig = ex.args[1]
    sig isa Expr && sig.head === :(::) && (sig = sig.args[1])
    return sig isa Expr && sig.head === :call ? sig.args[1] : nothing
end
function source_definition(path, name)
    tree = Meta.parseall(read(path, String); filename=path)
    matches = filter(ex -> definition_name(ex) === name, tree.args)
    length(matches) == 1 || error("Expected one $name definition in $path")
    return only(matches)
end
function load_definition!(caller, path, name)
    counts = Dict{Symbol, Int}()
    ex = redirect(source_definition(path, name), counts)
    counts == EXPECTED_SEAMS[name] || error("Loader boundary changed for $name: $counts")
    Core.eval(caller, ex)
end

struct InitialTimeSentinel end
struct AtmosphereSentinel end
function new_caller(route; guarded=false, env=Environment())
    caller = Module(gensym(:GramLoadingCaller))
    package = Module(gensym(:SyntheticSpaceAGORA))
    guarded && Core.eval(package, :(const GRAMSuite = $(env.gram)))
    sm = (InitialTime=InitialTimeSentinel, GRAMAtmosphereModel=AtmosphereSentinel)
    for (name, value) in (:SpaceAGORA => package, :SM => sm, :REPO_ROOT => REPO_SENTINEL,
                          :PPC_REPO_ROOT => REPO_SENTINEL, :LOAD_PATH => env.load_path,
                          :_GRAM_FIXTURE => env)
        Core.eval(caller, :(const $name = $value))
    end
    if route === :examples
        for name in (:_instantiate_vendored_gramsuite!, :_missing_dependency_error,
                     :ensure_gramsuite_loaded!, :setup_gram_example!)
            load_definition!(caller, EXAMPLES, name)
        end
    elseif route === :runtime
        load_definition!(caller, RUNTIME, :_ensure_perf_gramsuite_loaded!)
    elseif route === :package
        load_definition!(caller, PACKAGE, :ppc_ensure_gramsuite_loaded!)
    else
        error("Unknown route $route")
    end
    return caller, env
end
const ROUTE_FN = Dict(:examples => :ensure_gramsuite_loaded!, :runtime => :_ensure_perf_gramsuite_loaded!,
                      :package => :ppc_ensure_gramsuite_loaded!)
binding(caller, name) = Base.invokelatest(getfield, caller, name)
invoke(caller, name, args...) = Base.invokelatest(() -> getfield(caller, name)(args...))
load(caller, route) = invoke(caller, ROUTE_FN[route])
capture(f) = try f(); nothing catch err; err end
kinds(env) = first.(env.events)
activations(env) = [event[2] for event in env.events if event[1] === :activate]
imports(env) = filter(event -> event[1] === :import, env.events)

@testset verbose=true "GRAM loading: production control flow, synthetic environment" begin
    @testset "discovery and distinct load-path policies" begin
        for route in (:examples, :runtime, :package)
            caller, env = new_caller(route)
            @test load(caller, route) === nothing
            expected = route === :runtime ? [INITIAL_PATH; VENDOR_SENTINEL] : [VENDOR_SENTINEL; INITIAL_PATH]
            @test env.load_path == expected
            @test kinds(env) == [:find_package, :isdir, :import]
            @test only(imports(env))[3] == expected
            @test binding(caller, :GRAMSuite) === env.gram
            @test !isdefined(binding(caller, :SpaceAGORA), :GRAMSuite)

            caller, env = new_caller(route; env=Environment(discoverable=true))
            @test load(caller, route) === nothing
            @test env.load_path == INITIAL_PATH
            @test kinds(env) == [:find_package, :import]

            failure = ArgumentError("synthetic package unavailable")
            caller, env = new_caller(route; env=Environment(vendor_present=false, outcomes=Any[failure]))
            err = capture(() -> load(caller, route))
            @test env.load_path == INITIAL_PATH
            @test length(imports(env)) == 1
            @test !isdefined(caller, :GRAMSuite)
            @test isempty(activations(env))
            if route === :examples
                @test err isa ErrorException
                @test occursin(sprint(showerror, failure), sprint(showerror, err))
            else
                @test err === failure
            end
        end
    end

    @testset "caller binding does not satisfy the package-binding guard" begin
        for route in (:examples, :package)
            caller, env = new_caller(route)
            @test load(caller, route) === nothing
            first_path = copy(env.load_path)
            @test load(caller, route) === nothing
            @test length(imports(env)) == 2
            @test binding(caller, :GRAMSuite) === env.gram
            @test !isdefined(binding(caller, :SpaceAGORA), :GRAMSuite)
            @test env.load_path == first_path
            @test kinds(env) == [:find_package, :isdir, :import, :find_package, :import]

            caller, env = new_caller(route; guarded=true)
            @test load(caller, route) === nothing
            @test isempty(env.events)
            @test env.load_path == INITIAL_PATH
            @test !isdefined(caller, :GRAMSuite)
        end
    end

    @testset "repeated imports do not deduplicate an undiscoverable vendor" begin
        for route in (:examples, :runtime, :package)
            caller, env = new_caller(route; env=Environment(available_after_import=false))
            @test load(caller, route) === nothing
            @test load(caller, route) === nothing
            expected = route === :runtime ? [INITIAL_PATH; VENDOR_SENTINEL; VENDOR_SENTINEL] :
                                            [VENDOR_SENTINEL; VENDOR_SENTINEL; INITIAL_PATH]
            @test env.load_path == expected
            @test length(imports(env)) == 2
            @test kinds(env) == repeat([:find_package, :isdir, :import], 2)
            @test binding(caller, :GRAMSuite) === env.gram
            @test !isdefined(binding(caller, :SpaceAGORA), :GRAMSuite)
        end
    end

    @testset "example retry preserves cause and restores project" begin
        for message in ("Package Fixture is required but does not seem to be installed",
                        "Run `Pkg.instantiate()` to install all recorded dependencies")
            caller, env = new_caller(:examples; env=Environment(outcomes=Any[ErrorException(message), nothing]))
            prior = env.project
            @test load(caller, :examples) === nothing
            @test length(imports(env)) == 2
            @test activations(env) == [VENDOR_SENTINEL, dirname(prior)]
            @test env.project == prior
            @test binding(caller, :GRAMSuite) === env.gram
            @test kinds(env) == [:find_package, :isdir, :import, :isdir, :require_pkg,
                                 :active_project, :activate, :instantiate, :activate, :import]
        end
        first_error = ErrorException("Fixture is required but does not seem to be installed")
        retry_error = ArgumentError("synthetic retry failure")
        caller, env = new_caller(:examples; env=Environment(outcomes=Any[first_error, retry_error]))
        prior = env.project
        err = capture(() -> load(caller, :examples))
        @test err isa ErrorException
        @test occursin("after instantiating its vendored project", sprint(showerror, err))
        @test occursin(sprint(showerror, first_error), sprint(showerror, err))
        @test occursin(sprint(showerror, retry_error), sprint(showerror, err))
        @test length(imports(env)) == 2
        @test env.project == prior
        @test !isdefined(caller, :GRAMSuite)

        for (failure, vendor) in ((ArgumentError("unrelated import failure"), true), (first_error, false))
            caller, env = new_caller(:examples; env=Environment(vendor_present=vendor, outcomes=Any[failure]))
            err = capture(() -> load(caller, :examples))
            @test err isa ErrorException
            @test occursin("GRAM-backed examples require loading", sprint(showerror, err))
            @test occursin(sprint(showerror, failure), sprint(showerror, err))
            @test length(imports(env)) == 1
            @test isempty(activations(env))
        end
        for prior in (nothing, joinpath("synthetic", "prior", "Project.toml")), fails in (false, true)
            failure = ErrorException("synthetic instantiate failure")
            caller, env = new_caller(:examples; env=Environment(project=prior,
                outcomes=Any[first_error], instantiate_error=fails ? failure : nothing))
            err = capture(() -> load(caller, :examples))
            restore = prior === nothing ? REPO_SENTINEL : dirname(prior)
            @test err === (fails ? failure : nothing)
            @test activations(env) == [VENDOR_SENTINEL, restore]
            @test env.project == joinpath(restore, "Project.toml")
            @test length(imports(env)) == (fails ? 1 : 2)
            @test count(==(:instantiate), kinds(env)) == 1
        end
    end

    @testset "benchmark errors never install or retry" begin
        for route in (:runtime, :package)
            failure = ErrorException("Fixture is required but does not seem to be installed")
            caller, env = new_caller(route; env=Environment(outcomes=Any[failure]))
            prior = env.project
            @test capture(() -> load(caller, route)) === failure
            @test length(imports(env)) == 1
            @test kinds(env) == [:find_package, :isdir, :import]
            @test env.project == prior
            @test !isdefined(caller, :GRAMSuite)
        end
    end

    @testset "setup aliases use the requested caller and preserve existing bindings" begin
        caller, env = new_caller(:examples)
        @test invoke(caller, :setup_gram_example!) === nothing
        @test binding(caller, :InitialTime) === InitialTimeSentinel
        @test binding(caller, :GRAMAtmosphereModel) === AtmosphereSentinel
        target = Module(gensym(:AliasTarget))
        @test invoke(caller, :setup_gram_example!, target) === nothing
        @test binding(target, :InitialTime) === InitialTimeSentinel
        @test binding(target, :GRAMAtmosphereModel) === AtmosphereSentinel
        @test !isdefined(target, :GRAMSuite)
        target = Module(gensym(:ExistingAliasTarget))
        Core.eval(target, :(const InitialTime = :existing_epoch))
        Core.eval(target, :(const GRAMAtmosphereModel = :existing_atmosphere))
        @test invoke(caller, :setup_gram_example!, target) === nothing
        @test binding(target, :InitialTime) === :existing_epoch
        @test binding(target, :GRAMAtmosphereModel) === :existing_atmosphere
        @test length(imports(env)) == 3

        caller, env = new_caller(:examples; env=Environment(outcomes=Any[ErrorException("synthetic failure")]))
        target = Module(gensym(:FailedAliasTarget))
        @test capture(() -> invoke(caller, :setup_gram_example!, target)) isa ErrorException
        @test !isdefined(target, :InitialTime)
        @test !isdefined(target, :GRAMAtmosphereModel)
    end

    @testset "fixture boundary and host isolation" begin
        @test_throws ErrorException redirect(:(run(`false`)), Dict{Symbol, Int}())
        @test_throws ErrorException redirect(:(@eval import UnrelatedPackage), Dict{Symbol, Int}())
        @test_throws ErrorException redirect(:(import GRAMSuite), Dict{Symbol, Int}())
        @test Base.LOAD_PATH == REAL_PATH
        @test Base.active_project() == REAL_PROJECT
        @test pwd() == REAL_DIRECTORY
        @test ARGS == REAL_ARGS
        forbidden = Set(("SpaceAGORA", "GRAMSuite", "SPICE", "Pkg", "Distributed"))
        @test isempty(filter(id -> id.name in forbidden, setdiff(Set(keys(Base.loaded_modules)), REAL_PACKAGES)))
    end
end
end # module GramLoadingContractTests
