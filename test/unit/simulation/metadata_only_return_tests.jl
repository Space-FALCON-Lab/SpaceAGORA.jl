using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SM = SpaceAGORA.SimulationModel

# `return_solver_metadata` used to be reachable only under `return_solution`,
# and used to force per-step solution storage on. A caller that wanted the
# retcode and nothing else therefore paid for a whole trajectory and got
# `nothing` back for it. These tests pin the contract that replaced that:
# metadata is available on its own, it names the retcode directly, and the
# bundled form is unchanged for the callers that already use it.
function _metadata_test_config(; n_sats::Int=2, mission_s::Float64=120.0)
    planet = SM.Earth()
    spacecraft = SpacecraftModel[]
    for i in 1:n_sats
        root = Link(root=true, m=500.0, ref_area=12.0)
        ic = InitialCondition(
            ra=planet.Rp_e + 550_000.0 + 100.0 * i,
            rp=planet.Rp_e + 550_000.0 + 100.0 * i,
            i=53.0, ω=0.0, Ω=10.0, ν=360.0 * (i - 1) / n_sats
        )
        push!(spacecraft, SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0,
                                          root.inertia, 0, 0, ic, i))
    end
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_s, orientation_sim=false, num_steps_to_save=20
        ),
        environment_model=EnvironmentModel(
            planet=planet, EI=300.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false,
            ephemerides_model=SimpleEphemeridesModel()
        ),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=2.0)
    )
end

@testset "metadata-only return" begin
    args = _metadata_test_config()

    metadata = run_simulation(args; isolate_state=false, return_solver_metadata=true)

    # The whole point: a tuple, not `nothing`, and no trajectory inside it.
    @test metadata !== nothing
    @test metadata.solution === nothing
    @test metadata.retcode == "Success"
    @test metadata.solver_mode isa String
    @test metadata.solver_trace isa Vector
    @test !isempty(metadata.solver_trace)
    @test haskey(metadata, :parallel_policy)
    @test haskey(metadata, :spice_counters)
end

@testset "bundled form still carries the solution, and agrees" begin
    args = _metadata_test_config()

    bundled = run_simulation(args; isolate_state=false,
                             return_solution=true, return_solver_metadata=true)
    metadata = run_simulation(args; isolate_state=false, return_solver_metadata=true)

    @test bundled.solution !== nothing
    @test string(bundled.solution.retcode) == "Success"
    # `retcode` is new; it must say the same thing the solution did, so callers
    # can move off `result.solution.retcode` without changing meaning.
    @test bundled.retcode == string(bundled.solution.retcode)
    @test metadata.retcode == bundled.retcode
    @test metadata.solver_mode == bundled.solver_mode
end

@testset "the other two return shapes are untouched" begin
    args = _metadata_test_config()

    @test run_simulation(args; isolate_state=false) === nothing
    solution = run_simulation(args; isolate_state=false, return_solution=true)
    @test solution !== nothing
    @test string(solution.retcode) == "Success"
end
