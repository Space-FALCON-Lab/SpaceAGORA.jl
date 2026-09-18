module TerrainTouchdownTests
using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA
const SM = SpaceAGORA.SimulationModel
const CB = SM.SimulationCallbacks
const CH = SM.ControlHooks
const SE = SpaceAGORA.SimulationEngine
const ZERO3 = SVector(0.0, 0.0, 0.0)
const DIRECTION = SVector(cosd(30.0), 0.0, sind(30.0))

mutable struct ContactRecorder <: SpaceAGORA.AbstractControlEffectorModel
    selected::Int
    terrain::SM.AbstractTerrainModel
    radius::Float64
    clearance::Float64
    events::Vector{Any}
end
function CH.touchdown_spec(model::ContactRecorder, i::Int)
    i == model.selected || return nothing
    on_touchdown = (t, r_p, v_p, idx) -> push!(model.events, (t=t, r_p=r_p, v_p=v_p, i=idx))
    return (terrain=model.terrain, reference_radius_m=model.radius,
            height_m=model.clearance, on_touchdown=on_touchdown)
end
CH.calcControlEffect!(::ContactRecorder, u, p, t, i) = nothing
CH.calcControlForceTorque(::ContactRecorder, u, p, i, t) = (ZERO3, ZERO3)

struct FixedSpec <: SpaceAGORA.AbstractControlEffectorModel
    value::Any
end
CH.touchdown_spec(model::FixedSpec, i::Int) = i == 2 ? model.value : nothing

function configuration(; contact=true, second_altitude=52_000.0, mission_s=60.0)
    planet = SM.Earth()
    radius = planet.Rp_e - 2_000.0
    grid = SM.DEMGrid(Float32[3500 3500; 2500 2500], 20, 40, -5, 5;
                      reference_radius_m=radius)
    terrain = SM.DEMTerrainModel([grid]; reference_radius_m=radius)
    controller = ContactRecorder(2, terrain, radius, 5.0, Any[])
    spacecraft = SM.SpacecraftModel[]
    for (id, altitude) in ((41, 53_000.0), (97, second_altitude))
        root = SM.Link(root=true, m=500.0, ref_area=1.0)
        ic = SM.CartesianInitialCondition((planet.Rp_e + altitude) * DIRECTION, -1000.0 * DIRECTION)
        push!(spacecraft, SM.SpacecraftModel(SM.Joint[], [root], root, true,
            500.0, 0.0, root.inertia, 0, 0, ic, id))
    end
    args = SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime,
            mission_time=mission_s, number_of_orbits=1, keplerian=true,
            orientation_sim=false, num_steps_to_save=20),
        environment_model=SM.EnvironmentModel(planet=planet, EI=120.0,
            density_model=SM.NoAtmosphereModel(), topography=false, wind=false,
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            ephemerides_model=SM.SimpleEphemeridesModel(prime_meridian_at_reference_rad=0.0)),
        dynamics_model=SM.DynamicsModel(spacecraft, ()),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=contact ? (controller,) : (),
            control_rates=contact ? [1.0] : Float64[]),
        initial_time=SM.InitialTime(year=2000, month=1, day=1, hour=12, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-10,
            abstol_orbit=1e-10, dt_max_orbit=0.5),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))
    return args, controller
end

@testset "Terrain touchdown selection, event geometry and isolation" begin
    args, controller = configuration()
    initial_positions = [copy(sc.initial_condition.pos) for sc in args.dynamics_model.spacecraft]
    sol = SpaceAGORA.run_simulation(args; return_solution=true)
    planet = args.environment_model.planet
    # No force acts, so radial descent is linear. The north-south DEM is
    # h(latitude)=100*latitude: at planetocentric 30 degrees it is 3000 m.
    expected_time = (52_000.0 - (-2_000.0 + 3_000.0 + 5.0)) / 1000.0
    @test string(sol.retcode) == "Terminated"
    @test sol.t[end] ≈ expected_time atol=2e-7 rtol=0
    @test all(!, sol.prob.p.is_active)
    @test norm(SE._state_position_ii(sol.u[end], 1)) - planet.Rp_e ≈ 50_000.0 atol=2e-6
    @test norm(SE._state_position_ii(sol.u[end], 2)) - planet.Rp_e ≈ 1_005.0 atol=2e-6
    used = only(sol.prob.p.args.control_model.control_effectors)
    @test used !== controller
    @test length(used.events) == 1
    event = only(used.events)
    @test event.i == 2
    @test event.t ≈ expected_time atol=2e-7 rtol=0
    # Independent analytic rotating-frame oracle, including the spin transport
    # term. This detects returning inertial or merely rotated velocity.
    theta = planet.ω[3] * event.t
    rotation = @SMatrix [cos(theta) sin(theta) 0.0; -sin(theta) cos(theta) 0.0; 0.0 0.0 1.0]
    expected_r = rotation * ((planet.Rp_e + 1_005.0) * DIRECTION)
    expected_v = rotation * (-1000.0 * DIRECTION) - cross(planet.ω, expected_r)
    @test event.r_p ≈ expected_r atol=2e-6 rtol=0
    @test event.v_p ≈ expected_v atol=2e-8 rtol=0
    @test abs(event.v_p[2]) > 100.0
    @test isempty(controller.events)
    @test [sc.initial_condition.pos for sc in args.dynamics_model.spacecraft] == initial_positions

    # Once inactive, a craft contributes no new root and the contact hook is
    # not called again, even if the solver restarts its callbacks.
    callback = CB.get_touchdown_callback(CB._touchdown_specs(sol.prob.p.args, 2))
    integrator = (p=sol.prob.p, u=sol.u[end], t=sol.t[end])
    out = zeros(2)
    callback.condition(out, sol.u[end], sol.t[end], integrator)
    @test out == [1.0, 1.0]
    callback.affect_neg!(integrator, 2)
    @test length(used.events) == 1
end

@testset "No touchdown specification preserves the existing impact stop" begin
    args, _ = configuration(contact=false)
    @test all(isnothing, CB._touchdown_specs(args, 2))
    sol = SpaceAGORA.run_simulation(args; return_solution=true)
    @test string(sol.retcode) == "Terminated"
    @test sol.t[end] ≈ 3.0 atol=2e-7 rtol=0
    for i in 1:2
        @test norm(SE._state_position_ii(sol.u[end], i)) - args.environment_model.planet.Rp_e ≈ 50_000.0 atol=2e-6
    end
end

@testset "Touchdown specification validation and initial clearance" begin
    args, controller = configuration()
    original = CH.touchdown_spec(controller, 2)
    @test CH.touchdown_spec(controller, 1) === nothing
    @test CH.touchdown_spec(nothing, 2) === nothing
    duplicate = SM.SimConfig._with_configuration(args; control_model=SM.ControlModel(
        control_effectors=(controller, deepcopy(controller)), control_rates=[1.0, 1.0]))
    @test_throws ArgumentError CB._touchdown_specs(duplicate, 2)
    for bad in (42, (terrain=controller.terrain,),
                merge(original, (terrain=nothing,)),
                merge(original, (reference_radius_m=NaN,)),
                merge(original, (reference_radius_m=0.0,)),
                merge(original, (reference_radius_m=controller.radius + 1.0,)),
                merge(original, (height_m=-1.0,)),
                merge(original, (height_m=Inf,)),
                merge(original, (on_touchdown=nothing,)))
        invalid = SM.SimConfig._with_configuration(args; control_model=SM.ControlModel(
            control_effectors=(FixedSpec(bad),), control_rates=[1.0]))
        @test_throws ArgumentError CB._touchdown_specs(invalid, 2)
    end
    invalid_args, _ = configuration(second_altitude=500.0)
    @test_throws ArgumentError SpaceAGORA.run_simulation(invalid_args)
end
end # module TerrainTouchdownTests
