# The complete script of "Adding an Actuator of Your Own"
# (docs/src/user/custom_actuator.md), run by CI.
# The region between `# snippet: script` and its `# end snippet` is the page's
# script, and each nested region one of its shorter excerpts; the page shows each
# without the marker lines, and the HPC and extensibility docs gate compares them.
# snippet: script
using SpaceAGORA
using StaticArrays
using LinearAlgebra

# --- 1. The actuator kind ---------------------------------------------------------
struct ColdGasNozzle <: SpaceAGORA.AbstractActuator
    axis::SVector{3, Float64}        # thrust direction in the link frame (unit vector)
    location::SVector{3, Float64}    # nozzle position in the link frame, m
    max_thrust_n::Float64
    isp_s::Float64
end

SpaceAGORA.actuator_limits(n::ColdGasNozzle) = (0.0, n.max_thrust_n)

function SpaceAGORA.actuator_wrench(n::ColdGasNozzle, level, ctx::SpaceAGORA.ActuatorContext)
    force_body = ctx.link_to_body * (level * n.axis)            # newtons, body frame
    arm_body = ctx.link_to_body * n.location + ctx.link_r_body  # metres from the body reference point
    return force_body, cross(arm_body, force_body)
end

SpaceAGORA.actuator_consumption(n::ColdGasNozzle, level, ctx) =
    (mass_rate_kg_s = propellant_mass_rate(level, n.isp_s), power_w = 0.0)

# --- 2. A controller that commands it ----------------------------------------------
mutable struct BurnFor <: SpaceAGORA.AbstractControlEffectorModel
    thrust_n::Float64           # commanded thrust while burning, N
    burn_s::Float64             # burn duration, s
    levels::Vector{Float64}     # one level per mounted actuator, held between ticks
    initial_levels::Vector{Float64}  # the levels each run starts from
end

# Built with one vector of levels: each run starts from a copy of it.
BurnFor(thrust_n, burn_s, levels) = BurnFor(thrust_n, burn_s, collect(Float64, levels), collect(Float64, levels))

# snippet: levels
SpaceAGORA.actuation_mode(::BurnFor) = :mounted
SpaceAGORA.actuator_levels(c::BurnFor, spacecraft_index::Int) = c.levels
# end snippet

function SpaceAGORA.update_control!(c::BurnFor, u, p, t::Float64, spacecraft_index::Int)
    c.levels[1] = t < c.burn_s ? c.thrust_n : 0.0
    return nothing
end

# At each run's start, the levels the controller was built with.
function SpaceAGORA.initialize_control!(c::BurnFor, u, p, t)
    copy!(c.levels, c.initial_levels)
    return nothing
end

# The held levels through a checkpoint: saved with it, put back on a resumed run.
SpaceAGORA.checkpoint_state(c::BurnFor) = (levels=copy(c.levels),)
function SpaceAGORA.restore_checkpoint_state!(c::BurnFor, data)
    copy!(c.levels, data.levels)
    return nothing
end

# --- 3. A spacecraft carrying the nozzle, far from any planet ---------------------
root = Link(root=true, m=500.0, dims=MVector(1.0, 1.0, 1.0))
ic = CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), SVector(0.0, 0.0, 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
spacecraft = SpacecraftModel(links=[root], root=root, prop_mass=20.0, inertia_tensor=root.inertia,
                             initial_condition=ic, id=1)
nozzle = ColdGasNozzle(SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 0.5, 65.0)
add_actuator!(spacecraft, root, nozzle)

thrust_n, burn_s = 0.4, 300.0
planet = make_no_gram_planet(:earth)
config = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=true, verbose=false, results_directory=joinpath(@__DIR__, "output", "actuator_walkthrough"),
        generate_plots=false, generate_filenames=false, normalize=false, save_csv=true),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime, mission_time=600.0, number_of_orbits=1,
        keplerian=false, orientation_sim=true, num_steps_to_save=600, data_rate=1.0),
    environment_model=EnvironmentModel(
        planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false),
    dynamics_model=DynamicsModel([spacecraft], ()),              # no gravity: only the nozzle acts
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(BurnFor(thrust_n, burn_s, [thrust_n]),), control_rates=[1.0]),
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances())

results = run_simulation(config; return_results=true)

# --- 4. Check it by hand -----------------------------------------------------------------
last = results.table[end, :]
m0 = spacecraft.dry_mass + spacecraft.prop_mass
used_kg = propellant_mass_rate(thrust_n, nozzle.isp_s) * burn_s
dv_rocket = nozzle.isp_s * 9.80665 * log(m0 / (m0 - used_kg))
println("propellant used: $(round(m0 - last.sc1_mass, digits=6)) kg, expected $(round(used_kg, digits=6)) kg")
println("speed gained:    $(round(last.sc1_vel_1, digits=6)) m/s, rocket equation $(round(dv_rocket, digits=6)) m/s")
# end snippet

using Test
# CI's checks of what the page prints, outside the region the page shows: the
# propellant the engine took out of the vehicle and the speed it gained agree with
# the mass rate and the rocket equation to the six decimals printed.
@test isapprox(m0 - last.sc1_mass, used_kg; atol=1e-6)
@test isapprox(last.sc1_vel_1, dv_rocket; atol=1e-6)
