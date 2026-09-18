module DensityServiceRecipeTests

using Test
using Distributed
using Logging
using StaticArrays
using SpaceAGORA

const PP = SpaceAGORA.ParallelProcess
const CB = SpaceAGORA.SimulationModel.SimulationCallbacks
const EM = SpaceAGORA.SimulationModel.EnvironmentModels

# Stand-ins make construction settings observable without GRAM or native assets.
# These tests cover configuration isolation, not native random-stream equality.
recipe(; seed=17) = Dict{Symbol, Any}(
    :planet_name => "earth", :seed => seed,
    :initial_time => (year=2026, month=1, day=2, hour=3, minute=4, second=0.0),
    :gram_root_directory => "/fixture/gram", :gram_data_directory => "/fixture/data",
    :spice_directory => "/fixture/spice", :gram_library_path => "/fixture/libgram",
    :gram_min_relative_step_size => 0.02, :gram_perturbation_scales => (0.25, 0.5, 0.75, 0.0),
)

# Trigger an intervening request only when the coordinator copies the query
# vector, after ensure_density_workers! and before density_service_dispatch.
# This owns its array methods; no production method or hook is replaced.
struct RecipeSwitchVector{F} <: AbstractVector{Float64}
    values::Vector{Float64}
    on_first_read::F
    triggered::Base.RefValue{Bool}
end
Base.size(v::RecipeSwitchVector) = size(v.values)
Base.IndexStyle(::Type{<:RecipeSwitchVector}) = IndexLinear()
function Base.getindex(v::RecipeSwitchVector, i::Int)
    @boundscheck checkbounds(v.values, i)
    if !v.triggered[]
        v.triggered[] = true
        v.on_first_read()
    end
    return v.values[i]
end

function with_fixture(f)
    saved = (PP._DENSITY_SERVICE_BUILD_MODEL_FN[], PP._DENSITY_SERVICE_EVAL_FN[],
             PP._WORKER_DENSITY_MODEL[], PP._WORKER_DENSITY_RECIPE[])
    builds = Dict{Symbol, Any}[]
    active, peak = Ref(0), Ref(0)
    events = Tuple{Int, Float64}[]
    before_eval = Ref{Function}((_model, _h) -> nothing)
    warmup_fail = Ref(false) # fail the warm-up of ANY recipe while set: a setup fault, not a recipe key
    builder = function (kwargs)
        owned = deepcopy(kwargs)
        push!(builds, owned) # Model a native constructor mutating global state first.
        get(kwargs, :fail_build, false) && error("fixture late construction failure")
        # A constructor may consume/mutate its keyword container. Cache identity
        # must remain independent of both that container and the caller's one.
        kwargs[:builder_annotation] = true
        return (recipe=owned, id=length(builds))
    end
    evaluator = function (model, h, lat, lon, elapsed, wind, vacuum_temperature)
        active[] += 1
        peak[] = max(peak[], active[])
        try
            model.id == length(builds) || error("fixture stale native state")
            (get(model.recipe, :fail_warmup, false) || warmup_fail[]) && h == 1.0e5 && error("fixture warmup failure")
            before_eval[](model, h)
            yield() # Exercise task interleaving even on a one-thread test runner.
            get(model.recipe, :fail_query, false) && h == 12.0 && error("fixture query failure")
            h == 1.0e5 || push!(events, (model.id, h))
            return (Float64(get(model.recipe, :seed, 0)), vacuum_temperature + h,
                    (Float64(lat), Float64(lon), Float64(elapsed)))
        finally
            active[] -= 1
        end
    end
    lock(PP._WORKER_DENSITY_LOCK) do
        PP._WORKER_DENSITY_MODEL[] = nothing
        PP._WORKER_DENSITY_RECIPE[] = nothing
        PP._DENSITY_SERVICE_BUILD_MODEL_FN[] = builder
        PP._DENSITY_SERVICE_EVAL_FN[] = evaluator
    end
    PP.clear_density_service_failures!()
    try
        f((; builds, peak, events, before_eval, warmup_fail))
    finally
        PP.clear_density_service_failures!()
        lock(PP._WORKER_DENSITY_LOCK) do
            PP._DENSITY_SERVICE_BUILD_MODEL_FN[] = saved[1]
            PP._DENSITY_SERVICE_EVAL_FN[] = saved[2]
            PP._WORKER_DENSITY_MODEL[] = saved[3]
            PP._WORKER_DENSITY_RECIPE[] = saved[4]
        end
    end
end

batch(kwargs; hs=[11.0, 12.0, 13.0]) = PP.density_batch_remote(
    hs, fill(2.0, length(hs)), fill(3.0, length(hs)), fill(4.0, length(hs)), true, 200.0, kwargs)

@testset "density service construction recipes" begin
    @testset "owned recipe, reuse and invalidation" begin
        with_fixture() do ctx
            original = recipe()
            @test PP.install_worker_density_model!(original)
            first_model = PP._WORKER_DENSITY_MODEL[]
            @test length(ctx.builds) == 1
            @test isequal(PP._WORKER_DENSITY_RECIPE[], original)
            @test !haskey(PP._WORKER_DENSITY_RECIPE[], :builder_annotation)
            @test PP._WORKER_DENSITY_RECIPE[] !== original
            @test PP._WORKER_DENSITY_RECIPE[] !== first_model.recipe
            @test PP.install_worker_density_model!(deepcopy(original))
            @test PP._WORKER_DENSITY_MODEL[] === first_model
            @test length(ctx.builds) == 1
            original[:initial_time] = merge(original[:initial_time], (second=12.0,))
            original[:gram_perturbation_scales] = (0.8, 0.5, 0.75, 0.0)
            @test PP._WORKER_DENSITY_RECIPE[][:initial_time].second == 0
            @test first_model.recipe[:gram_perturbation_scales][1] == 0.25
            @test PP.install_worker_density_model!(original)
            @test PP._WORKER_DENSITY_MODEL[] !== first_model
            @test isequal(PP._WORKER_DENSITY_RECIPE[], original)
            for (key, value) in (
                :planet_name => "mars", :seed => 29,
                :initial_time => (year=2027, month=1, day=2, hour=3, minute=4, second=0.0),
                :gram_root_directory => "/other/gram", :gram_data_directory => "/other/data",
                :spice_directory => "/other/spice", :gram_library_path => "/other/libgram",
                :gram_min_relative_step_size => 0.04, :gram_perturbation_scales => (1.0, 0.5, 0.75, 0.0),
            )
                PP.install_worker_density_model!(original)
                changed = deepcopy(original)
                changed[key] = value
                previous = PP._WORKER_DENSITY_MODEL[]
                count = length(ctx.builds)
                @test PP.install_worker_density_model!(changed)
                @test PP._WORKER_DENSITY_MODEL[] !== previous
                @test length(ctx.builds) == count + 1
                @test isequal(PP._WORKER_DENSITY_MODEL[].recipe, changed)
                @test PP.install_worker_density_model!(deepcopy(changed))
                @test length(ctx.builds) == count + 1
            end
        end
    end

    @testset "construction and warmup failures do not publish" begin
        with_fixture() do ctx
            good = recipe()
            for flag in (:fail_build, :fail_warmup)
                @test_throws ErrorException PP.install_worker_density_model!(merge(good, Dict(flag => true)))
                @test PP._WORKER_DENSITY_MODEL[] === nothing
                @test PP._WORKER_DENSITY_RECIPE[] === nothing
            end
            @test PP.install_worker_density_model!(good)
            for flag in (:fail_build, :fail_warmup)
                previous = PP._WORKER_DENSITY_MODEL[]
                bad = merge(good, Dict(flag => true))
                @test_throws ErrorException PP.install_worker_density_model!(bad)
                @test PP._WORKER_DENSITY_MODEL[] === nothing
                @test PP._WORKER_DENSITY_RECIPE[] === nothing
                @test_throws ErrorException batch(bad)
                @test PP._WORKER_DENSITY_MODEL[] === nothing
                @test PP._WORKER_DENSITY_RECIPE[] === nothing
                @test all(==(17.0), fetch(@async batch(good))[1])
                @test PP._WORKER_DENSITY_MODEL[] !== previous
            end
            @test_throws ErrorException batch(merge(good, Dict(:fail_query => true)))
            @test all(==(17.0), fetch(@async batch(good))[1])
            @test !islocked(PP._WORKER_DENSITY_LOCK)
        end
    end

    @testset "legacy planet API and alternating requests" begin
        with_fixture() do ctx
            @test PP.install_worker_density_model!("earth")
            @test PP._WORKER_DENSITY_RECIPE[] == Dict(:planet_name => "earth")
            @test PP.install_worker_density_model!(Dict(:planet_name => "earth"))
            @test length(ctx.builds) == 1
            @test PP.density_batch_remote([11.0], [0.0], [0.0], [0.0], false, 200.0)[1] == [0.0]
            a, b = recipe(seed=17), recipe(seed=29)
            @test PP.install_worker_density_model!(a)
            @test PP.install_worker_density_model!(b) # Another request after ensure(A).
            result = batch(a) # Dispatch(A) must select A again, not evaluate B.
            @test result == (fill(17.0, 3), [211.0, 212.0, 213.0], fill(2.0, 3), fill(3.0, 3), fill(4.0, 3))
            @test all(==(29.0), batch(b)[1])
            @test all(==(17.0), batch(a)[1])
        end
    end

    @testset "yielding requests hold one model for the entire batch" begin
        with_fixture() do ctx
            entered, release, started = Channel{Nothing}(1), Channel{Nothing}(1), Channel{Nothing}(1)
            ctx.before_eval[] = function (model, h)
                if get(model.recipe, :seed, 0) == 17 && h == 11.0
                    put!(entered, nothing)
                    take!(release)
                end
            end
            a = @async try
                batch(recipe(seed=17))
            finally
                # If construction fails before the evaluator signals, let the
                # test report that failure instead of blocking on take! forever.
                isready(entered) || put!(entered, nothing)
            end
            b = nothing
            try
                take!(entered)
                istaskdone(a) && fetch(a)
                b = @async begin
                    put!(started, nothing)
                    batch(recipe(seed=29))
                end
                take!(started)
                yield()
                @test !istaskdone(b)
                @test length(ctx.builds) == 1
            finally
                put!(release, nothing)
                # Join both tasks even if the first failed. fetch below (or the
                # active exception above) reports failure after state is safe.
                for task in (a, b)
                    task === nothing && continue
                    try
                        wait(task)
                    catch
                    end
                end
            end
            @test all(==(17.0), fetch(a)[1])
            @test all(==(29.0), fetch(b)[1])
            @test ctx.peak[] == 1
            @test first.(ctx.events) == [1, 1, 1, 2, 2, 2]
        end
    end

    @testset "dispatch and coordinator fallback publish no partial batch" begin
        with_fixture() do ctx
            hs = [11.0, 12.0, 13.0, 14.0]
            ranges = UnitRange{Int}[1:2, 3:4]
            args = ([myid(), myid()], ranges, hs, zeros(4), zeros(4), zeros(4), true, 200.0)
            result = PP.density_service_dispatch(args...; constructor_kwargs=recipe(seed=29))
            @test result !== nothing
            @test vcat((part[1] for part in result)...) == fill(29.0, 4)
            bad = merge(recipe(), Dict(:fail_query => true))
            @test_logs (:warn, r"Density service worker failed on a batch") begin
                @test PP.density_service_dispatch(args...; constructor_kwargs=bad) === nothing
            end
            # Use the current process as a test-owned pool entry: exercise the
            # actual coordinator and local dispatch paths without spawning or GRAM.
            pool = PP.density_process_pool()
            old_workers = copy(pool.workers)
            old_atexit = PP._DENSITY_ATEXIT_REGISTERED[]
            lock(pool.lock) do
                empty!(pool.workers)
                push!(pool.workers, myid())
            end
            PP._DENSITY_ATEXIT_REGISTERED[] = true
            try
                withenv("SPACEAGORA_GRAM_PROCESS_POOL" => "on",
                        "SPACEAGORA_GRAM_PROCESS_POOL_WORKERS" => "1",
                        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
                    model = EM.GRAMAtmosphereModel(nothing, ReentrantLock(), bad)
                    rho, temp = fill(-1.0, 4), fill(-2.0, 4)
                    winds = fill(SVector(-3.0, -3.0, -3.0), 4)
                    p = (args=(environment_model=(planet=(T_ref=200.0,),),),)
                    @test_logs (:warn, r"Density service worker failed on a batch") begin
                        @test !CB._gram_process_pool_batch_eval!(rho, temp, winds, model,
                            hs, zeros(4), zeros(4), 0.0, true, p)
                    end
                    @test rho == fill(-1.0, 4)
                    @test temp == fill(-2.0, 4)
                    @test winds == fill(SVector(-3.0, -3.0, -3.0), 4)

                    # The worker can be reassigned between ensure(A) and dispatch.
                    # Only the actual coordinator's recipe keyword restores A.
                    a, b = recipe(seed=17), recipe(seed=29)
                    before = length(ctx.builds)
                    switched = Ref(false)
                    switch_recipe = function ()
                        @test isequal(PP._WORKER_DENSITY_RECIPE[], a)
                        @test PP.install_worker_density_model!(b)
                    end
                    intervened_hs = RecipeSwitchVector(copy(hs), switch_recipe, switched)
                    good_model = EM.GRAMAtmosphereModel(nothing, ReentrantLock(), a)
                    @test CB._gram_process_pool_batch_eval!(rho, temp, winds, good_model,
                        intervened_hs, fill(2.0, 4), fill(3.0, 4), 4.0, true, p)
                    @test switched[]
                    @test rho == fill(17.0, 4)
                    @test temp == 200.0 .+ hs
                    @test winds == fill(SVector(2.0, 3.0, 4.0), 4)
                    @test [r[:seed] for r in ctx.builds[before+1:end]] == [17, 29, 17]
                    @test isequal(PP._WORKER_DENSITY_RECIPE[], a)
                end
            finally
                lock(pool.lock) do
                    empty!(pool.workers)
                    append!(pool.workers, old_workers)
                end
                PP._DENSITY_ATEXIT_REGISTERED[] = old_atexit
            end
        end
    end

    @testset "persistent setup failure is remembered; recovery is explicit" begin
        with_fixture() do ctx
            hs = [11.0, 12.0, 13.0, 14.0]
            rho, temp = fill(-1.0, 4), fill(-2.0, 4)
            winds = fill(SVector(-3.0, -3.0, -3.0), 4)
            p = (args=(environment_model=(planet=(T_ref=200.0,),),),)
            batch_eval(model) = CB._gram_process_pool_batch_eval!(rho, temp, winds, model,
                hs, zeros(4), zeros(4), 0.0, true, p)
            pool = PP.density_process_pool()
            old_workers = copy(pool.workers)
            old_atexit = PP._DENSITY_ATEXIT_REGISTERED[]
            lock(pool.lock) do
                empty!(pool.workers)
                push!(pool.workers, myid())
            end
            PP._DENSITY_ATEXIT_REGISTERED[] = true
            try
                withenv("SPACEAGORA_GRAM_PROCESS_POOL" => "on",
                        "SPACEAGORA_GRAM_PROCESS_POOL_WORKERS" => "1",
                        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
                    r = recipe(seed=31)
                    model = EM.GRAMAtmosphereModel(nothing, ReentrantLock(), r)
                    @test isempty(PP.density_service_failures())

                    # Persistent worker setup fault: the first batch pays one
                    # construction attempt and warns once; later batches decline
                    # without contacting the worker and without a log line.
                    ctx.warmup_fail[] = true
                    before = length(ctx.builds)
                    @test_logs (:warn, r"could not build its GRAM instance") (:warn, r"No density worker could build") begin
                        @test !batch_eval(model)
                    end
                    @test length(ctx.builds) == before + 1
                    @test PP.density_service_failed(r)
                    @test haskey(PP.density_service_failures(), r)
                    @test_logs min_level=Logging.Warn begin
                        @test !batch_eval(model)
                        @test !batch_eval(model)
                    end
                    @test length(ctx.builds) == before + 1
                    @test rho == fill(-1.0, 4) && temp == fill(-2.0, 4)
                    # The RHS prefill gate declines the remembered recipe up front.
                    @test !CB._rhs_density_service_candidate((args=(environment_model=(density_model=model,),),), 4)

                    # A changed configuration is a new key and is tried once more.
                    ctx.warmup_fail[] = false
                    other = EM.GRAMAtmosphereModel(nothing, ReentrantLock(), recipe(seed=32))
                    @test batch_eval(other)
                    @test rho == fill(32.0, 4)
                    @test PP.density_service_failed(r)

                    # Fixing the worker setup alone does not retry: recovery is explicit.
                    @test !batch_eval(model)
                    @test length(ctx.builds) == before + 2
                    PP.clear_density_service_failures!()
                    @test !PP.density_service_failed(r)
                    @test batch_eval(model)
                    @test rho == fill(31.0, 4)
                    @test length(ctx.builds) == before + 3

                    # A batch failure is remembered the same way, and a pool restart clears it.
                    bad = merge(recipe(seed=33), Dict(:fail_query => true))
                    bad_model = EM.GRAMAtmosphereModel(nothing, ReentrantLock(), bad)
                    @test_logs (:warn, r"Density service worker failed on a batch") begin
                        @test !batch_eval(bad_model)
                    end
                    @test PP.density_service_failed(bad)
                    events_before = length(ctx.events)
                    @test_logs min_level=Logging.Warn begin
                        @test !batch_eval(bad_model)
                    end
                    @test length(ctx.events) == events_before
                    lock(pool.lock) do
                        empty!(pool.workers) # nothing to remove; the restart only has to clear the memo
                    end
                    PP.shutdown_density_workers!()
                    @test isempty(PP.density_service_failures())
                end
            finally
                lock(pool.lock) do
                    empty!(pool.workers)
                    append!(pool.workers, old_workers)
                end
                PP._DENSITY_ATEXIT_REGISTERED[] = old_atexit
            end
        end
    end

    @testset "raw-core wrapper declines before worker acquisition" begin
        with_fixture() do ctx
            raw = EM.GRAMAtmosphereModel((planet_name="earth",))
            rho, temp, winds = [-1.0], [-2.0], [SVector(-3.0, -3.0, -3.0)]
            existing = copy(PP.density_process_pool().workers)
            processes = procs()
            withenv("SPACEAGORA_GRAM_PROCESS_POOL" => "on",
                    "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
                @test_logs min_level=Logging.Warn begin
                    @test !CB._gram_process_pool_batch_eval!(rho, temp, winds, raw,
                        [11.0], [0.0], [0.0], 0.0, true, nothing)
                end
                # No callback cache fields exist: the recipe guard precedes them.
                @test !CB._rhs_density_service_candidate((args=(environment_model=(density_model=raw,),),), 1)
            end
            @test isempty(ctx.builds)
            @test PP.density_process_pool().workers == existing
            @test procs() == processes
            @test rho == [-1.0] && temp == [-2.0]
            @test winds == [SVector(-3.0, -3.0, -3.0)]
        end
    end
end

end # module DensityServiceRecipeTests
