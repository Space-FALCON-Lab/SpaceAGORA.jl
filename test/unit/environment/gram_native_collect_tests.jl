module GramNativeCollectTests
# Discarded native GRAM atmospheres and the collection that releases them.
#
# A native GRAM atmosphere is freed by its wrapper's finalizer, and the Julia
# collector cannot see the native memory it holds, so a process that discards
# one per run accumulates them (measured at about 106 MB resident per Earth
# atmosphere; docs/architecture/gram_thread_scaling.md, "Pool-worker memory
# growth"). `run_simulation` calls `collect_unreferenced_gram_atmospheres!`
# after every run, and the GRAM extension runs a full collection once enough
# native atmospheres have accumulated since its previous one.
#
# The first testset needs no GRAM: it checks that the engine calls the hook
# after a run, whether the run succeeds or throws. The rest need a native GRAM
# on the host and are skipped without one.
using Test
using SpaceAGORA
using StaticArrays

const SM = SpaceAGORA.SimulationModel
const EM = SM.EnvironmentModels

const NC_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const NC_GRAM_ROOT = joinpath(NC_REPO, "data", "GRAMSuite.jl")
const NC_SPICE_PATH = joinpath(NC_GRAM_ROOT, "GRAM Suite 2.0", "SPICE")
const NC_GRAM_LIB = joinpath(NC_GRAM_ROOT, "GRAM Suite 2.0", "Build", "lib", "libGRAM.so")
const NC_GRAM_READY = isfile(NC_GRAM_LIB) && isdir(NC_SPICE_PATH)
const EPOCH = SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)

# A density model whose epoch alignment can be made to throw, so a run can be
# failed before it solves anything.
struct HookDensity <: SpaceAGORA.AbstractDensityModel
    fail::Bool
end
function EM.with_density_model_epoch(model::HookDensity, epoch)
    model.fail && error("alignment failure requested by the test")
    return model
end
EM.getDensity(::HookDensity, h::Float64, lat::Float64, lon::Float64,
    elapsed::Float64, wind::Bool, p) = (0.0, 300.0, SVector(0.0, 0.0, 0.0))

function configuration(density)
    planet = SM.make_no_gram_planet(:earth)
    root = SM.Link(root=true, m=500.0, ref_area=12.0)
    ic = SM.InitialCondition(ra=planet.Rp_e + 550e3, rp=planet.Rp_e + 550e3,
        i=53.0, ω=0.0, Ω=10.0, ν=0.0)
    spacecraft = SM.SpacecraftModel(SM.Joint[], [root], root, true, 500.0,
        0.0, root.inertia, 0, 0, ic, 1)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime,
            mission_time=1.0, num_steps_to_save=4, data_rate=2.0),
        environment_model=SM.EnvironmentModel(planet=planet, EI=300.0,
            density_model=density, thermal_model=SM.MaxwellianHeat(
                thermal_accomodation_factor=1.0, planet=planet),
            topography=false, topo_degree=12, topo_order=8, wind=false,
            ephemerides_model=SM.SimpleEphemeridesModel()),
        dynamics_model=SM.DynamicsModel([spacecraft], (SM.InverseSquaredGravityModel(),)),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=EPOCH,
        integration_tolerances=SM.IntegrationTolerances(dt_max_orbit=0.5),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))
end

@testset "run_simulation calls the native-atmosphere collection after every run" begin
    calls = Ref(0)
    saved = EM._COLLECT_UNREFERENCED_GRAM_ATMOSPHERES_FN[]
    EM._COLLECT_UNREFERENCED_GRAM_ATMOSPHERES_FN[] = () -> (calls[] += 1; false)
    try
        sol = SpaceAGORA.run_simulation(configuration(HookDensity(false));
            isolate_state=false, return_solution=true, visualization=false)
        @test string(sol.retcode) == "Success"
        @test calls[] == 1
        SpaceAGORA.run_simulation(configuration(HookDensity(false)); visualization=false)
        @test calls[] == 2
        # A run that throws still discards whatever it built.
        @test_throws ErrorException SpaceAGORA.run_simulation(
            configuration(HookDensity(true)); visualization=false)
        @test calls[] == 3
    finally
        EM._COLLECT_UNREFERENCED_GRAM_ATMOSPHERES_FN[] = saved
    end
end

if !NC_GRAM_READY
    @info "Skipping native GRAM collection tests: no libGRAM on this host." lib = NC_GRAM_LIB
else
    if Base.find_package("GRAMSuite") === nothing
        pushfirst!(LOAD_PATH, NC_GRAM_ROOT)
    end
    @eval import GRAMSuite
    const EXT = Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt)

    # Sets the accounting to a known state; the counters are process-global.
    function set_counts!(created, finalized, baseline)
        EXT._NATIVE_ATMOSPHERES_CREATED[] = created
        EXT._NATIVE_ATMOSPHERES_FINALIZED[] = finalized
        EXT._NATIVE_ATMOSPHERES_BASELINE[] = baseline
    end

    @testset "collection threshold, baseline and override" begin
        @test EM._COLLECT_UNREFERENCED_GRAM_ATMOSPHERES_FN[] === EXT._collect_unreferenced_native_atmospheres!
        saved = (EXT._NATIVE_ATMOSPHERES_CREATED[], EXT._NATIVE_ATMOSPHERES_FINALIZED[],
                 EXT._NATIVE_ATMOSPHERES_BASELINE[])
        try
            withenv("SPACEAGORA_GRAM_NATIVE_COLLECT_LIMIT" => nothing) do
                @test EXT._native_collect_limit() == EXT._NATIVE_COLLECT_LIMIT_DEFAULT == 8
                # Seven more alive than after the last collection: no collection.
                set_counts!(1007, 1000, 0)
                @test !EM.collect_unreferenced_gram_atmospheres!()
                @test EXT._NATIVE_ATMOSPHERES_BASELINE[] == 0
                # Eight: one collection, after which the baseline is the live count,
                # so atmospheres that are still referenced do not re-trigger it.
                set_counts!(1008, 1000, 0)
                @test EM.collect_unreferenced_gram_atmospheres!()
                @test EXT._NATIVE_ATMOSPHERES_BASELINE[] == EXT._native_atmospheres_live()
                @test !EM.collect_unreferenced_gram_atmospheres!()
            end
            withenv("SPACEAGORA_GRAM_NATIVE_COLLECT_LIMIT" => "3") do
                @test EXT._native_collect_limit() == 3
                set_counts!(1003, 1000, 0)
                @test EM.collect_unreferenced_gram_atmospheres!()
            end
            for off in ("0", "off", "OFF")
                withenv("SPACEAGORA_GRAM_NATIVE_COLLECT_LIMIT" => off) do
                    @test EXT._native_collect_limit() == 0
                    set_counts!(5000, 0, 0)
                    @test !EM.collect_unreferenced_gram_atmospheres!()
                end
            end
            for bad in ("-2", "many")
                withenv("SPACEAGORA_GRAM_NATIVE_COLLECT_LIMIT" => bad) do
                    @test EXT._native_collect_limit() == EXT._NATIVE_COLLECT_LIMIT_DEFAULT
                end
            end
        finally
            set_counts!(saved...)
        end
    end

    # Builds `n` Earth models and lets them go. Each keyword construction and each
    # deepcopy is a native atmosphere of its own.
    function discard_models!(n)
        for _ in 1:n
            model = EM.GRAMAtmosphereModel(planet_name="earth", initial_time=EPOCH)
            deepcopy(model)
        end
        return nothing
    end

    @testset "discarded native atmospheres are counted and released" begin
        withenv("SPACEAGORA_GRAM_NATIVE_COLLECT_LIMIT" => "4") do
            GC.gc(true)
            EXT._NATIVE_ATMOSPHERES_BASELINE[] = EXT._native_atmospheres_live()
            created0 = EXT._NATIVE_ATMOSPHERES_CREATED[]
            finalized0 = EXT._NATIVE_ATMOSPHERES_FINALIZED[]
            # On a task of its own, so no stack slot of this one still holds
            # the last model when the collection runs.
            wait(@async discard_models!(2))
            @test EXT._NATIVE_ATMOSPHERES_CREATED[] - created0 == 4
            @test EM.collect_unreferenced_gram_atmospheres!()
            # The collection ran the wrappers' finalizers: all four are freed.
            @test EXT._NATIVE_ATMOSPHERES_FINALIZED[] - finalized0 == 4
            @test EXT._NATIVE_ATMOSPHERES_BASELINE[] == EXT._native_atmospheres_live()
        end
    end

    # The parallelization_performance harness builds its cached GRAM model at
    # the epoch its cases fly, so a run's epoch alignment keeps that model
    # instead of rebuilding a native atmosphere for every sample.
    const HARNESS = Module(:HarnessEpoch)
    Core.eval(HARNESS, :(using SpaceAGORA))
    for file in ("cli.jl", "modes.jl", "cases.jl")
        Base.include(HARNESS, joinpath(NC_REPO, "benchmarks", "studies",
            "parallelization_performance", file))
    end

    @testset "harness GRAM model is built at the harness epoch" begin
        H = HARNESS
        model = H.ppc_gram_atmosphere_model("earth")
        @test model === H.ppc_gram_atmosphere_model("earth")
        args = H.ppc_single_config("aero_4096sat_l50_gram_process_100s",
            H.PPCConfig(profile="full"); seed=1, mc_index=7)
        @test args.environment_model.density_model === model
        @test EM._density_epoch_key(args.initial_time) == EM._density_epoch_key(H.ppc_initial_time())
        created0 = EXT._NATIVE_ATMOSPHERES_CREATED[]
        aligned = SpaceAGORA.SimulationEngine._with_density_model_epoch(args)
        @test aligned.environment_model.density_model === model
        @test EXT._NATIVE_ATMOSPHERES_CREATED[] == created0
    end
end
end # module GramNativeCollectTests
