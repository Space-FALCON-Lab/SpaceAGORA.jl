# Actual one-worker service acceptance. Requires native GRAM; the coverage driver
# runs it only under its native gate and reports completion after it returns.
module GRAMDensityServiceProbes
using Test, Distributed, SpaceAGORA, StaticArrays
if Base.find_package("GRAMSuite") === nothing
    pushfirst!(LOAD_PATH, joinpath(dirname(dirname(pathof(SpaceAGORA))), "data", "GRAMSuite.jl"))
end
import GRAMSuite
const SM = SpaceAGORA.SimulationModel
const EM = SM.EnvironmentModels
const CB = SM.SimulationCallbacks
const PP = SpaceAGORA.ParallelProcess
isempty(PP.density_process_pool().workers) || error("Native density probe requires an empty density pool")
# Load the dynamically included native wrapper before compiling the testset,
# so its wind-state bindings are visible to the local reference evaluator.
EM.GRAMAtmosphereModel(planet_name="mars")

function worker_snapshot(worker)
    remotecall_fetch(Core.eval, worker, Main, quote
        let pp = SpaceAGORA.ParallelProcess, m = SpaceAGORA.ParallelProcess._WORKER_DENSITY_MODEL[]
            m === nothing && return (recipe=deepcopy(pp._WORKER_DENSITY_RECIPE[]), handle=nothing)
            (recipe=deepcopy(pp._WORKER_DENSITY_RECIPE[]),
             constructed_recipe=deepcopy(m.constructor_kwargs),
             handle=objectid(m.gram_atmosphere), epoch=m.initial_time,
             paths=(m.gram_root, m.gram_data_root, m.spice_root))
        end
    end)
end

@testset "Native density service preserves and replaces recipes" begin
    @test isempty(PP.density_process_pool().workers) # this probe owns only its pool
    prior_workers = Set(Distributed.workers())
    withenv("SPACEAGORA_GRAM_PROCESS_POOL"=>"on", "SPACEAGORA_GRAM_PROCESS_POOL_WORKERS"=>"1",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE"=>nothing, "SPACEAGORA_GRAM_STATIC_GRID"=>"0",
            "SPACEAGORA_GRAM_STATIC_GRID_PREBUILD_ALL_PLANETS"=>"0",
            "SPACEAGORA_GRAM_OFFLINE_SURROGATE"=>"off", "SPACEAGORA_GRAM_WIND_MODE"=>"nominal") do
        # Mean-state Mars F10.7 and epoch discriminate an incorrect default rebuild.
        a = EM.GRAMAtmosphereModel(planet_name="mars", mars_f107=80.0, seed=7319,
            gram_min_relative_step_size=0.02, gram_perturbation_scales=(0.0,0.0,0.0,0.0),
            initial_time=GRAMSuite.InitialTime(year=2001, month=11, day=6))
        SM.Mars("", a.spice_root)
        recipe_a = deepcopy(a.constructor_kwargs)
        recipe_b = merge(recipe_a, Dict{Symbol,Any}(:mars_f107=>180.0, :seed=>42,
            :gram_min_relative_step_size=>0.03,
            :initial_time=>GRAMSuite.InitialTime(year=2002, month=7, day=1)))
        b = EM.GRAMAtmosphereModel(; recipe_b...)
        p = (args=(environment_model=(planet=SM.Mars(),),),)
        hs, lats, lons, times = [150e3, 180e3], [0.1, -0.2], [0.2, 0.5], [10.0, 30.0]
        rhos, Ts = fill(-1.0, 2), fill(-2.0, 2)
        winds = fill(SVector{3,Float64}(-3.0,-3.0,-3.0), 2)
        batch(model) = CB._gram_process_pool_batch_eval!(rhos, Ts, winds,
            model, hs, lats, lons, times, true, p)
        reference(model) = [PP._DENSITY_SERVICE_EVAL_FN[](model,
            hs[i],lats[i],lons[i],times[i],true,p.args.environment_model.planet.T_ref) for i in eachindex(hs)]
        function matches(values)
            all(isapprox(rhos[i], values[i][1]; rtol=1e-12, atol=0.0) &&
                isapprox(Ts[i], values[i][2]; rtol=1e-12, atol=0.0) &&
                isapprox(winds[i], values[i][3]; rtol=1e-12, atol=1e-12) for i in eachindex(hs))
        end
        try
            raw = EM.GRAMAtmosphereModel(a.core)
            @test !batch(raw)
            @test isempty(PP.density_process_pool().workers)
            @test rhos == [-1.0,-1.0] && Ts == [-2.0,-2.0]
            @test !CB._rhs_density_service_candidate((args=(environment_model=(density_model=raw,),),), 2)
            ref_a, ref_b = reference(a), reference(b)
            @test all(x -> isfinite(x[1]) && x[1] > 0.0, vcat(ref_a, ref_b))
            @test any(!isapprox(ref_a[i][1],ref_b[i][1];rtol=1e-6) for i in eachindex(hs))
            # Process bootstrap furnishes Earth kernels only. Non-Earth and
            # custom SPICE setup remains the caller's existing responsibility.
            worker = only(PP.ensure_process_workers!(PP.density_process_pool(), 1))
            remotecall_fetch(Core.eval, worker, Main,
                :(SpaceAGORA.SimulationModel.Mars("", $(a.spice_root)); nothing))
            @test batch(a)
            first_a = worker_snapshot(worker)
            @test isequal(first_a.recipe, recipe_a)
            @test isequal(first_a.constructed_recipe, recipe_a)
            @test matches(ref_a)
            @test batch(a)
            @test worker_snapshot(worker).handle == first_a.handle
            @test batch(b)
            first_b = worker_snapshot(worker)
            @test first_b.handle != first_a.handle
            @test isequal(first_b.recipe, b.constructor_kwargs)
            @test first_b.epoch == b.initial_time
            @test first_b.paths == (b.gram_root,b.gram_data_root,b.spice_root)
            @test matches(ref_b)
            # Another setup call can replace the model before dispatch. The
            # batch must carry and restore A instead of silently evaluating B.
            values = PP.density_service_dispatch([worker],[1:2],hs,lats,lons,times,true,
                p.args.environment_model.planet.T_ref; constructor_kwargs=recipe_a)
            @test values !== nothing
            if values !== nothing
                @test all(isapprox(values[1][1][i],ref_a[i][1];rtol=1e-12) for i in eachindex(hs))
            end
            second_a = worker_snapshot(worker)
            @test second_a.handle != first_b.handle
            @test isequal(second_a.recipe,recipe_a)
            @test isequal(a.constructor_kwargs,recipe_a)
            # Retain the old native handle on this test-owned worker so failure
            # recovery cannot reuse a collected object's address. Compare live
            # objects below instead of objectid values after cache invalidation.
            remotecall_fetch(Core.eval, worker, Main, quote
                global _density_probe_prior_handle =
                    SpaceAGORA.ParallelProcess._WORKER_DENSITY_MODEL[].gram_atmosphere
                nothing
            end)
            # Unsupported planet selection throws after library selection and
            # native initialization. Recovery must reconstruct the old recipe.
            bad = merge(recipe_a, Dict{Symbol,Any}(:planet_name=>"unsupported"))
            failed = PP.density_service_dispatch([worker],[1:2],hs,lats,lons,times,true,
                p.args.environment_model.planet.T_ref; constructor_kwargs=bad)
            @test failed === nothing
            @test worker_snapshot(worker).handle === nothing
            @test worker_snapshot(worker).recipe === nothing
            remotecall_fetch(GC.gc, worker)
            @test batch(a)
            @test remotecall_fetch(Core.eval, worker, Main, quote
                SpaceAGORA.ParallelProcess._WORKER_DENSITY_MODEL[].gram_atmosphere !==
                    _density_probe_prior_handle
            end)
            @test isequal(worker_snapshot(worker).recipe,recipe_a)
            @test matches(ref_a)
        finally
            PP.shutdown_density_workers!()
        end
    end
    @test isempty(PP.density_process_pool().workers)
    @test Set(Distributed.workers()) == prior_workers
end
end
