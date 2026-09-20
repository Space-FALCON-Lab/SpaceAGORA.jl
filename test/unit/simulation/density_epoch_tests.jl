module DensityEpochTests
using Test, SpaceAGORA, StaticArrays
const SM = SpaceAGORA.SimulationModel
const EM = SM.EnvironmentModels
const SE = SpaceAGORA.SimulationEngine
const OLD = SM.InitialTime(year=2020, month=1, day=2, hour=3, minute=4, second=5.5)
const NEW = SM.InitialTime(year=2021, month=2, day=3, hour=4, minute=5, second=6.25)

# Exercise the actual epoch implementation with a local constructor stand-in.
# No methods of the real GRAM constructor or extension are replaced.
module RecipeFixture
using SpaceAGORA
const InitialTime = SpaceAGORA.SimulationModel.InitialTime
abstract type AbstractDensityModel end
struct PlainDensity <: AbstractDensityModel end
struct GRAMAtmosphereModel <: AbstractDensityModel
    core
    instance_lock::ReentrantLock
    constructor_kwargs::Union{Nothing, Dict{Symbol, Any}}
end
struct GRAMAtmosphereModelSurrogate{M} <: AbstractDensityModel
    base_model::M
    surrogate_file::String
    point_fallback_below_m::Union{Nothing, Float64}
end
_gram_not_loaded_error(name) = error("Unexpected extension fallback: $name")
include(joinpath(dirname(pathof(SpaceAGORA)), "environment", "atmosphere", "density_epoch.jl"))
const CONSTRUCTIONS = Ref(0)
const FAIL_CONSTRUCTION = Ref(false)
function _rebuild_gram_epoch_model(recipe::Dict{Symbol, Any})
    CONSTRUCTIONS[] += 1
    FAIL_CONSTRUCTION[] && error("fixture constructor failure")
    return GRAMAtmosphereModel((initial_time=deepcopy(recipe[:initial_time]),),
        ReentrantLock(), deepcopy(recipe))
end
end

@testset "density epoch recipe and fixed-table contract (native-free)" begin
    F = RecipeFixture
    recipe = Dict{Symbol, Any}(:initial_time => OLD, :planet_name => "mars",
        :gram_root_directory => "/owned/root", :gram_data_directory => "/owned/data",
        :spice_directory => "/owned/spice", :gram_library_path => "/owned/libGRAM.so",
        :gram_perturbation_scales => [0.0, 0.25, 0.5, 0.75],
        :gram_min_relative_step_size => 0.02, :mars_f107 => 111.0)
    model = F.GRAMAtmosphereModel((initial_time=OLD,), ReentrantLock(), deepcopy(recipe))
    F.CONSTRUCTIONS[] = 0
    @test F.with_density_model_epoch(model, deepcopy(OLD)) === model
    @test F.CONSTRUCTIONS[] == 0
    # Field-compatible epochs are compared at the native stored precision.
    equivalent = (; year=2020, month=1, day=2, hour=3, minute=4, second=5.5 + eps(5.5))
    @test F.with_density_model_epoch(model, equivalent) === model
    @test F.CONSTRUCTIONS[] == 0
    nearby = SM.InitialTime(year=2020, month=1, day=2, hour=3, minute=4,
        second=nextfloat(OLD.second))
    rebuilt = F.with_density_model_epoch(model, nearby)
    @test F.CONSTRUCTIONS[] == 1
    @test F._density_epoch_key(rebuilt.core.initial_time) == F._density_epoch_key(nearby)
    @test rebuilt.instance_lock !== model.instance_lock
    @test rebuilt.constructor_kwargs !== model.constructor_kwargs
    for name in keys(recipe)
        name === :initial_time && continue
        @test rebuilt.constructor_kwargs[name] == recipe[name]
    end
    @test rebuilt.constructor_kwargs[:gram_perturbation_scales] !== model.constructor_kwargs[:gram_perturbation_scales]
    rebuilt.constructor_kwargs[:gram_perturbation_scales][1] = 99.0
    @test model.constructor_kwargs == recipe
    @test model.core.initial_time == OLD
    @test_throws ArgumentError F.with_density_model_epoch(model, (;year=2020))
    @test_throws ArgumentError F.with_density_model_epoch(model,
        (;year=2020, month=1, day=2, hour=3, minute=4, second=NaN))
    @test F.CONSTRUCTIONS[] == 1

    raw = F.GRAMAtmosphereModel(model.core, ReentrantLock(), nothing)
    @test F.with_density_model_epoch(raw, OLD) === raw
    @test_throws ArgumentError F.with_density_model_epoch(raw, NEW)
    @test F.CONSTRUCTIONS[] == 1
    plain = F.PlainDensity()
    @test F.with_density_model_epoch(plain, NEW) === plain
    for fallback in (nothing, 1234.0), base in (model, raw)
        surrogate = F.GRAMAtmosphereModelSurrogate(base, "fixed.arrow", fallback)
        @test F.with_density_model_epoch(surrogate, OLD) === surrogate
        @test surrogate.surrogate_file == "fixed.arrow"
        @test surrogate.point_fallback_below_m === fallback
        @test_throws ArgumentError F.with_density_model_epoch(surrogate, NEW)
    end
    @test F.CONSTRUCTIONS[] == 1  # Reject before rebuilding any native fallback.
    custom_surrogate = F.GRAMAtmosphereModelSurrogate(plain, "custom.arrow", nothing)
    @test F.with_density_model_epoch(custom_surrogate, NEW) === custom_surrogate
    F.FAIL_CONSTRUCTION[] = true
    try
        @test_throws ErrorException F.with_density_model_epoch(model, NEW)
    finally
        F.FAIL_CONSTRUCTION[] = false
    end
    @test model.constructor_kwargs == recipe
    @test model.core.initial_time == OLD

    # The real owner/root API is present without loading GRAMSuite.
    @test SpaceAGORA.with_density_model_epoch === EM.with_density_model_epoch
    @test EM.with_density_model_epoch(SM.NoAtmosphereModel(), NEW) isa SM.NoAtmosphereModel
    actual_raw = EM.GRAMAtmosphereModel((initial_time=OLD,))
    @test EM.with_density_model_epoch(actual_raw, OLD) === actual_raw
    @test_throws ArgumentError EM.with_density_model_epoch(actual_raw, NEW)
    for fallback in (nothing, 1234.0)
        actual_surrogate = EM.GRAMAtmosphereModelSurrogate(actual_raw, "fixed.arrow", fallback)
        @test EM.with_density_model_epoch(actual_surrogate, OLD) === actual_surrogate
        @test_throws ArgumentError EM.with_density_model_epoch(actual_surrogate, NEW)
    end
end

# A custom epoch-aware atmosphere verifies the real engine's public hook and
# copy ordering, without native assets or altering global GRAM methods.
struct EpochDensity <: SpaceAGORA.AbstractDensityModel
    initial_time::SM.InitialTime
    owned::Vector{Int}
end
const ALIGN_INPUTS = EpochDensity[]
const ALIGN_OUTPUTS = EpochDensity[]
function EM.with_density_model_epoch(model::EpochDensity, epoch)
    push!(ALIGN_INPUTS, model)
    aligned = EM._density_epoch_key(model.initial_time) == EM._density_epoch_key(epoch) ?
        model : EpochDensity(deepcopy(epoch), deepcopy(model.owned))
    push!(ALIGN_OUTPUTS, aligned)
    return aligned
end
EM.getDensity(::EpochDensity, h::Float64, lat::Float64, lon::Float64,
    elapsed::Float64, wind::Bool, p) = (0.0, 300.0, SVector(0.0, 0.0, 0.0))

function configuration(density; epoch=NEW)
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
        initial_time=epoch,
        integration_tolerances=SM.IntegrationTolerances(dt_max_orbit=0.5),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))
end

@testset "epoch alignment preserves configuration and run isolation" begin
    original = configuration(EpochDensity(OLD, [7]))
    aligned = SE._with_density_model_epoch(original)
    @test aligned !== original
    @test aligned.environment_model.density_model.initial_time == NEW
    @test original.environment_model.density_model.initial_time == OLD
    for name in fieldnames(typeof(original))
        name === :environment_model && continue
        @test getfield(aligned, name) === getfield(original, name)
    end
    for name in fieldnames(typeof(original.environment_model))
        name === :density_model && continue
        @test getfield(aligned.environment_model, name) === getfield(original.environment_model, name)
    end
    @test SE._with_density_model_epoch(aligned) === aligned
    plain = configuration(SM.NoAtmosphereModel())
    @test SE._with_density_model_epoch(plain) === plain

    for same_epoch in (false, true), isolate in (false, true)
        model = EpochDensity(same_epoch ? NEW : OLD, [7])
        args = configuration(model)
        empty!(ALIGN_INPUTS); empty!(ALIGN_OUTPUTS)
        sol = SpaceAGORA.run_simulation(args; isolate_state=isolate,
            return_solution=true, visualization=false)
        @test string(sol.retcode) == "Success"
        @test sol.t[end] == 1.0
        @test length(ALIGN_INPUTS) == length(ALIGN_OUTPUTS) == 1
        @test only(ALIGN_INPUTS) === model  # Hook precedes isolation.
        rebuilt = only(ALIGN_OUTPUTS)
        used = sol.prob.p.args.environment_model.density_model
        @test used.initial_time == NEW
        @test model.initial_time == (same_epoch ? NEW : OLD)
        @test (rebuilt === model) == same_epoch
        @test (used === rebuilt) == !isolate
        @test (used.owned === model.owned) == (same_epoch && !isolate)
        push!(used.owned, 9)
        @test model.owned == (same_epoch && !isolate ? [7, 9] : [7])
        @test sol.prob.p.args.solver_config == args.solver_config
    end
end
end # module DensityEpochTests
