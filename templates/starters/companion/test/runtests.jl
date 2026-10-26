using Test, SpaceAGORA, StarterName
using SpaceAGORA.Conformance: test_service_facade

const EP = SpaceAGORA.ExternalPropagationServices
const MU = 3.986004418e14

# One ordinary spacecraft and one on a circular orbit that the companion owns.
function owned_configuration(; mission_time_s=20.0)
    planet = make_no_gram_planet(:earth)
    spacecraft = map(1:2) do id
        root = Link(root=true, m=50.0, ref_area=1.0)
        radius = planet.Rp_e + 600e3 + 100e3 * id
        speed = sqrt(MU / radius)
        initial = CartesianInitialCondition([radius, 0.0, 0.0], [0.0, speed, 0.0])
        SpacecraftModel(joints=Joint[], links=[root], root=root, prop_mass=0.0,
            inertia_tensor=root.inertia, initial_condition=initial, id=id)
    end
    configuration = SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true,
            number_of_orbits=1, mission_time=mission_time_s, orientation_sim=false, num_steps_to_save=20),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
    return EP.with_external_propagators(configuration, (CircularOrbitPropagator([2]; step_s=0.5, mu=MU),))
end

@testset "StarterName" begin
    # The facade this companion serves has the shape the adapter standard
    # requires, at the contract version the companion was written for.
    test_service_facade(EP; version=StarterName.CORE_CONTRACT_VERSION)

    @testset "the owned spacecraft follows its circle in a run" begin
        solution = run_simulation(owned_configuration(); return_solution=true)
        @test last(solution.t) ≈ 20.0
        @test all(isfinite, solution.u[end])
        final = solution.u[end].sc[2]
        @test isapprox(sqrt(sum(abs2, final.pos)), make_no_gram_planet(:earth).Rp_e + 800e3; rtol=1e-9)
    end

    @testset "a core-only session loads no companion code" begin
        script = joinpath(@__DIR__, "core_only.jl")
        @test success(pipeline(`$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $script`;
            stdout=stdout, stderr=stderr))
    end
end
