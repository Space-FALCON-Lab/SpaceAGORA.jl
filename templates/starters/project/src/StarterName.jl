"""
    StarterName

A research study on SpaceAGORA's public API. `LowThrustPush` is the study's own
force model, an example to replace; `run_case` runs the study's case. Every
SpaceAGORA name here is a root export.
"""
module StarterName

using SpaceAGORA

export LowThrustPush, case_configuration, run_case

"""
    LowThrustPush(thrust_n)

A constant force of `thrust_n` newtons along the inertial velocity, with no
torque: an example force model to replace with the study's own.
"""
struct LowThrustPush <: AbstractForceTorqueModel
    thrust_n::Float64
end

# The force (inertial frame, N) and torque (body frame, N m) on the spacecraft.
function SpaceAGORA.wrench(model::LowThrustPush, x::StateSample, env::EnvironmentSample, t::Float64)
    speed = sqrt(sum(abs2, x.vel_ii))
    return (model.thrust_n / speed) * x.vel_ii, zero(x.vel_ii)
end

"""
    case_configuration(; mission_time_s=60.0, thrust_n=1.0e-3)

The study's case: one 100 kg spacecraft in a 650 x 700 km Earth orbit under J2
gravity and `LowThrustPush`, with no atmosphere. It needs no native GRAM, SPICE
kernels or private data.
"""
function case_configuration(; mission_time_s=60.0, thrust_n=1.0e-3)
    planet = make_no_gram_planet(:earth)
    root = Link(root=true, m=100.0, ref_area=2.0)
    initial = InitialCondition(ra=planet.Rp_e + 700e3, rp=planet.Rp_e + 650e3, i=45.0, ω=0.0, Ω=0.0, ν=10.0)
    spacecraft = SpacecraftModel(joints=Joint[], links=[root], root=root, prop_mass=0.0,
        inertia_tensor=root.inertia, initial_condition=initial, id=1)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true,
            number_of_orbits=1, mission_time=mission_time_s, orientation_sim=false, num_steps_to_save=20),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(SpacecraftModel[spacecraft],
            (InverseSquaredJ2GravityModel(), LowThrustPush(thrust_n))),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

"""
    run_case(; kwargs...)

Run [`case_configuration`](@ref) and return the solution.
"""
run_case(; kwargs...) = run_simulation(case_configuration(; kwargs...); return_solution=true)

end # module StarterName
