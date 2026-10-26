module ConfigFieldCopyTests
# Copies of configuration structs name only the fields they change
# (`SimConfig._with_fields`): every other field reaches the copy as it was, the
# same object or the same bits, including fields added to a type later.
using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SM = SpaceAGORA.SimulationModel
const SC = SpaceAGORA.SimulationCampaigns
const SE = SpaceAGORA.SimulationEngine
const SV = SM.SceneVisualization
const TV = SpaceAGORA.TelemetryVerification
const MW = SpaceAGORA.MissionWorkflow
const CF = SpaceAGORA.Conformance
const with_fields = SM.SimConfig._with_fields

const PLANET = make_no_gram_planet(:earth)
const ROOT = Link(root=true, m=500.0, ref_area=4.0)
const IC = InitialCondition(a=PLANET.Rp_e + 4.0e5, e=0.0, i=53.0, Ω=0.0, ω=0.0, ν=0.0)

# A value away from its default for every field of every struct that a copy site
# copies. A field added to one of these types fails the first test set until it is
# listed here; from then on every copy below is checked for it.
const ODD = (
    SimulationSettings=(results=false, verbose=true, results_directory="odd_results", generate_plots=false,
        generate_filenames=true, normalize=true, save_csv=false, save_visualization_scene=true,
        checkpoint_enabled=true, checkpoint_interval_s=123.0, checkpoint_directory="odd_checkpoints",
        resume_from_checkpoint=true, articulated_live_pose_loads=true, results_layout=:long),
    MissionConfiguration=(mission_type=MissionOrbits, keplerian=false, number_of_orbits=3, mission_time=4321.0,
        orientation_sim=true, num_steps_to_save=500, data_rate=2.5),
    EnvironmentModel=(planet=PLANET, EI=250.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(), topography=true, topo_degree=12, topo_order=11, wind=false,
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=0.5, planet=PLANET)),
    SpacecraftModel=(joints=Joint[], links=[ROOT], root=ROOT, instant_actuation=false, dry_mass=500.0,
        prop_mass=7.25, inertia_tensor=ROOT.inertia, initial_condition=IC, id=7,
        # Non-empty GNC, so a copy without GNC that fails to strip it is caught.
        guidance=GuidanceModel((:g,), [2.0]), navigation=NavigationModel((:n,), [3.0]),
        control=ControlModel((:c,), [4.0]), attachments=SM.CompliantAttachment[], n_terminal=2,
        battery_energy_index=3.5, tempurature_index=4.5, constellation_slot=ConstellationSlot(2, 3)),
    MonteCarloSampleResult=(index=3, seed=17, success=false, elapsed_s=1.25, value=:value, error=:error,
        backtrace=:backtrace, finished_ns=99.0),
)
const TYPES = (SimulationSettings=SimulationSettings, MissionConfiguration=MissionConfiguration,
    EnvironmentModel=EnvironmentModel, SpacecraftModel=SpacecraftModel,
    MonteCarloSampleResult=SC.MonteCarloSampleResult)
# Each type built with its defaults, and the fields it has no default for.
const DEFAULTS = (SimulationSettings=SimulationSettings(), MissionConfiguration=MissionConfiguration(),
    EnvironmentModel=EnvironmentModel(planet=PLANET, EI=250.0, density_model=ODD.EnvironmentModel.density_model,
        thermal_model=ODD.EnvironmentModel.thermal_model),
    SpacecraftModel=SpacecraftModel(),
    MonteCarloSampleResult=SC.MonteCarloSampleResult(index=3, seed=17, success=false, elapsed_s=1.25))
const REQUIRED = (EnvironmentModel=(:planet, :EI, :density_model, :thermal_model),
    MonteCarloSampleResult=(:index, :seed, :success, :elapsed_s))

odd(name::Symbol) = name === :SpacecraftModel ? SpacecraftModel(values(ODD[name])...) :
    TYPES[name](; ODD[name]...)

plain(v) = isbits(v) || v isa AbstractString || v isa Symbol

# The fields of `copy` that do not hold what the copy should: each `changed` field
# its given value (or passing the given predicate), every other field the same as
# in `source` (`===`: the same object, or the same bits).
function wrong_fields(copy, source; changed=(;))
    wrong = Symbol[]
    for f in fieldnames(typeof(copy))
        v = getfield(copy, f)
        ok = if haskey(changed, f)
            expected = changed[f]
            expected isa Function ? expected(v) : isequal(v, expected)
        else
            v === getfield(source, f)
        end
        ok || push!(wrong, f)
    end
    return wrong
end

const SETTINGS = odd(:SimulationSettings)
const MISSION = odd(:MissionConfiguration)
const ENVIRONMENT = odd(:EnvironmentModel)

function walker4()
    factory = (plane, slot, index) -> begin
        root = Link(root=true, m=500.0, ref_area=4.0)
        SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, IC, index)
    end
    return walker_constellation(factory; total=4, planes=2, phasing=1, altitude=4.0e5, planet=PLANET, i=53)
end

function odd_configuration(settings=SETTINGS)
    spacecraft = walker4()
    return SimulationConfiguration(;
        simulation_settings=settings,
        mission_configuration=MISSION,
        environment_model=ENVIRONMENT,
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(),)),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-10, dt_max_orbit=1.0))
end
const ARGS = odd_configuration()

@testset "the table names every field of every copied struct, away from its default" begin
    for name in keys(TYPES)
        T = TYPES[name]
        @test Set(fieldnames(T)) == Set(keys(ODD[name]))
        instance = odd(name)
        @test instance isa T
        required = get(REQUIRED, name, ())
        at_default = [f for f in fieldnames(T) if !(f in required) && plain(getfield(instance, f)) &&
                      isequal(getfield(instance, f), getfield(DEFAULTS[name], f))]
        @test isempty(at_default)
    end
end

@testset "_with_fields changes the named fields and passes on the rest" begin
    for name in keys(TYPES)
        source = odd(name)
        @test isempty(wrong_fields(with_fields(source), source))
    end
    copy = with_fields(SETTINGS; verbose=false, results_directory="elsewhere")
    @test isempty(wrong_fields(copy, SETTINGS; changed=(verbose=false, results_directory="elsewhere")))
    # A spacecraft is rebuilt in field order: its keyword constructor computes fields.
    spacecraft = odd(:SpacecraftModel)
    @test isempty(wrong_fields(with_fields(spacecraft; id=8), spacecraft; changed=(id=8,)))
    # The constructor runs: its checks apply, and type parameters follow the new values.
    @test_throws ArgumentError with_fields(MISSION; mission_time=-1.0)
    density = ExponentialAtmosphereModel(1.0e-11, 4.0e5, 6.0e4; temperature_k=1000.0,
        valid_min_altitude_m=0.0, valid_max_altitude_m=2.0e6)
    @test with_fields(ENVIRONMENT; density_model=density).density_model === density
    @test_throws "SimulationSettings has no field `result`" with_fields(SETTINGS; result=true)
    # Negative control: a copy that names every field but one loses that one.
    dropping(s) = SimulationSettings(; (f => getfield(s, f) for f in fieldnames(SimulationSettings)
                                         if f !== :articulated_live_pose_loads)...)
    @test wrong_fields(dropping(SETTINGS), SETTINGS) == [:articulated_live_pose_loads]
    @test isempty(wrong_fields(SM.SimConfig._with_configuration(ARGS; solver_config=SolverConfig()), ARGS;
        changed=(solver_config=c -> c isa SolverConfig,)))
end

@testset "ensemble members keep every setting but their folders" begin
    tag = "sat_3_id_3"
    # Checkpoints in an explicit folder: only that folder is split.
    @test isempty(wrong_fields(SC._ensemble_member_settings(SETTINGS, tag), SETTINGS;
        changed=(checkpoint_directory=joinpath("odd_checkpoints", tag),)))
    # Results written, checkpoints in the default folder: only the results folder is split.
    results = with_fields(SETTINGS; results=true, checkpoint_directory="")
    @test isempty(wrong_fields(SC._ensemble_member_settings(results, tag), results;
        changed=(results_directory=joinpath("odd_results", tag),)))
    nothing_to_split = with_fields(SETTINGS; checkpoint_enabled=false, resume_from_checkpoint=false)
    @test SC._ensemble_member_settings(nothing_to_split, tag) === nothing_to_split
    for (index, spacecraft) in enumerate(ARGS.dynamics_model.spacecraft)
        member = SC._ensemble_member_configuration(ARGS, spacecraft, index, 4)
        member_tag = SC._ensemble_member_tag(index, 4, spacecraft.id)
        @test isempty(wrong_fields(member.simulation_settings, SETTINGS;
            changed=(checkpoint_directory=joinpath("odd_checkpoints", member_tag),)))
        @test isempty(wrong_fields(member, ARGS; changed=(simulation_settings=s -> true,
            dynamics_model=d -> d.spacecraft == [spacecraft])))
    end
end

@testset "scene and viewer-ensemble copies" begin
    @test SV.with_visualization_scene(ARGS, true) === ARGS
    @test isempty(wrong_fields(SV.with_visualization_scene(ARGS, false).simulation_settings, SETTINGS;
        changed=(save_visualization_scene=false,)))
    for visualization in (true, false)
        moved = SV.with_results_directory(ARGS, "elsewhere"; visualization=visualization)
        @test isempty(wrong_fields(moved.simulation_settings, SETTINGS;
            changed=(results=true, results_directory="elsewhere", save_visualization_scene=visualization)))
    end
end

@testset "the example smoke copy keeps every setting it does not mean to change" begin
    for (keep, output) in ((true, ""), (false, "smoke_output"))
        withenv("SPACEAGORA_EXAMPLE_SMOKE" => "1", "SPACEAGORA_EXAMPLE_SMOKE_RESULTS" => keep ? "1" : "0",
                "SPACEAGORA_CLI_OUTPUT_DIR" => output, "SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME" => "120.0") do
            smoke = TV._example_smoke_args(ARGS)
            @test isempty(wrong_fields(smoke.simulation_settings, SETTINGS; changed=(results=keep, verbose=false,
                results_directory=isempty(output) ? joinpath(pwd(), "output") : abspath(output),
                generate_plots=false, normalize=false, save_csv=keep, save_visualization_scene=false,
                checkpoint_enabled=false, resume_from_checkpoint=false)))
            # articulated_live_pose_loads, a physics setting, is among the fields kept.
            @test smoke.simulation_settings.articulated_live_pose_loads
            @test isempty(wrong_fields(smoke.mission_configuration, MISSION;
                changed=(number_of_orbits=1, mission_time=120.0, num_steps_to_save=200)))
        end
    end
    withenv("SPACEAGORA_EXAMPLE_SMOKE" => "0") do
        @test TV._example_smoke_args(ARGS) === ARGS
    end
end

@testset "telemetry verification copies" begin
    @test isempty(wrong_fields(TV._with_environment_wind(ARGS, true).environment_model, ENVIRONMENT;
        changed=(wind=true,)))
    @test isempty(wrong_fields(TV._with_orbit_mission(ARGS, 5, 999.0).mission_configuration, MISSION;
        changed=(mission_type=MissionOrbits, number_of_orbits=5, mission_time=999.0)))
    @test isempty(wrong_fields(TV._with_study_settings(ARGS; quick=true).mission_configuration, MISSION;
        changed=(num_steps_to_save=2000,)))
end

@testset "mission workflow and conformance copies" begin
    settings = with_fields(SETTINGS; checkpoint_enabled=false, resume_from_checkpoint=false)
    args = odd_configuration(settings)
    plan = MW.MissionPlan("field-copy", (; kwargs...) -> args)
    destination = joinpath(mktempdir(), "mission")
    @test isempty(wrong_fields(MW._configuration(plan, (;), destination, true).simulation_settings, settings;
        changed=(results=true, results_directory=destination, save_visualization_scene=true)))
    for resume in (true, false)
        @test isempty(wrong_fields(CF._checkpointed(ARGS, "ck"; resume=resume).simulation_settings, SETTINGS;
            changed=(checkpoint_enabled=true, checkpoint_interval_s=MISSION.mission_time / 2,
                checkpoint_directory="ck", resume_from_checkpoint=resume)))
    end
end

@testset "a spacecraft without its GNC keeps everything else" begin
    spacecraft = odd(:SpacecraftModel)
    @test isempty(wrong_fields(SE._without_gnc(spacecraft), spacecraft; changed=(
        guidance=g -> isempty(g.guidance_effectors), navigation=n -> isempty(n.navigation_effectors),
        control=c -> isempty(c.control_effectors))))
end

@testset "Monte Carlo sample copies" begin
    sample = odd(:MonteCarloSampleResult)
    @test isempty(wrong_fields(SC._stamp_finished(sample), sample;
        changed=(finished_ns=t -> isfinite(t) && t != 99.0,)))
    @test isempty(wrong_fields(SC._reindexed_sample(sample, 9), sample; changed=(index=9,)))
end

end # module ConfigFieldCopyTests
