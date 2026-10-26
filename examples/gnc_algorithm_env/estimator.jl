# Worked example B of "Writing a GNC Algorithm" (docs/src/user/gnc_feedback.md):
# an estimator that holds a rate estimate from a gyro, and a controller that reads
# it through a declared link. Each region between a `# snippet:` line and
# `# end snippet` appears on that page as it is here; the HPC and extensibility
# docs gate compares them.
using SpaceAGORA
using SpaceAGORA.Conformance: conformance_fixture, test_control_model, test_navigation_model
using StaticArrays
using Test

# snippet: link
# The estimate the navigation model writes and the controller reads: one shared object.
mutable struct RateEstimate
    ω_body::SVector{3, Float64}     # estimated body rate, rad/s
    time_s::Float64                 # time of the measurement behind it, s into the run
end
RateEstimate() = RateEstimate(zero(SVector{3, Float64}), NaN)

# One estimator writes it; any number of controllers read it. Each run starts empty.
SpaceAGORA.link_cardinality(::RateEstimate) = (producers = 1:1, consumers = 0:typemax(Int))
function SpaceAGORA.link_initialize!(e::RateEstimate, producers, consumers, u, p, t)
    e.ω_body, e.time_s = zero(SVector{3, Float64}), NaN
    return nothing
end
# end snippet

# snippet: estimator
struct GyroRateEstimator <: SpaceAGORA.AbstractNavigationModel
    spacecraft_index::Int
    gyro::Int                       # the gyro's position among the spacecraft's sensors
    estimate::RateEstimate
end
SpaceAGORA.spacecraft_scope(n::GyroRateEstimator) = BoundSpacecraft((n.spacecraft_index,))
SpaceAGORA.gnc_links(n::GyroRateEstimator) = (Produces(n.estimate),)

# The tick: read the latest measurement and hold it as the estimate. Sensors
# sample before navigation ticks at the same time.
function SpaceAGORA.update_navigation!(n::GyroRateEstimator, u, p, t, spacecraft_index)
    spacecraft = p.args.dynamics_model.spacecraft[spacecraft_index]
    n.estimate.ω_body = SVector{3, Float64}(sensor_measurement(spacecraft, n.gyro))
    n.estimate.time_s = sensor_measurement_time(spacecraft, n.gyro)
    return nothing
end

# The producer saves the estimate. Restore its values into the existing link so
# every consumer still reads the same shared object.
SpaceAGORA.checkpoint_state(n::GyroRateEstimator) =
    (ω_body=n.estimate.ω_body, time_s=n.estimate.time_s)
function SpaceAGORA.restore_checkpoint_state!(n::GyroRateEstimator, data)
    n.estimate.ω_body, n.estimate.time_s = data.ω_body, data.time_s
    return nothing
end
# end snippet

# snippet: consumer
# The rate damper of example A, reading the estimate instead of the true state.
Base.@kwdef mutable struct EstimatedRateDamper <: SpaceAGORA.AbstractControlEffectorModel
    spacecraft_index::Int = 1
    gain_n_m_s::Float64 = 2.0
    max_torque_n_m::Float64 = 1.0
    period_s::Float64 = 1.0
    estimate::RateEstimate
    torque_body::SVector{3, Float64} = zero(SVector{3, Float64})
    tick_s::Float64 = NaN
end
SpaceAGORA.spacecraft_scope(c::EstimatedRateDamper) = BoundSpacecraft((c.spacecraft_index,))
SpaceAGORA.gnc_period_s(c::EstimatedRateDamper) = c.period_s
SpaceAGORA.gnc_links(c::EstimatedRateDamper) = (Consumes(c.estimate),)
SpaceAGORA.held_command(c::EstimatedRateDamper, spacecraft_index) = c.torque_body
SpaceAGORA.held_command_time(c::EstimatedRateDamper, spacecraft_index) = c.tick_s
function SpaceAGORA.release_command!(c::EstimatedRateDamper, t, spacecraft_index)
    c.torque_body, c.tick_s = zero(SVector{3, Float64}), t
    return nothing
end
function SpaceAGORA.initialize_control!(c::EstimatedRateDamper, u, p, t)
    c.torque_body, c.tick_s = zero(SVector{3, Float64}), NaN
    return nothing
end

function SpaceAGORA.update_control!(c::EstimatedRateDamper, u, p, t, spacecraft_index)
    c.torque_body = clamp.(-c.gain_n_m_s .* c.estimate.ω_body, -c.max_torque_n_m, c.max_torque_n_m)
    c.tick_s = t
    return nothing
end

function SpaceAGORA.control_wrench(c::EstimatedRateDamper, x::StateSample, p, t)
    zero3 = zero(SVector{3, Float64})
    x.spacecraft_index == c.spacecraft_index || return zero3, zero3
    return zero3, c.torque_body
end

# The controller saves its held command; the producer owns the shared estimate.
SpaceAGORA.checkpoint_state(c::EstimatedRateDamper) =
    (torque_body=c.torque_body, tick_s=c.tick_s)
function SpaceAGORA.restore_checkpoint_state!(c::EstimatedRateDamper, data)
    c.torque_body, c.tick_s = data.torque_body, data.tick_s
    return nothing
end
# end snippet

# snippet: run
root = Link(root=true, m=500.0, dims=MVector(1.0, 1.0, 1.0))
ω0 = 0.05
ic = CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), SVector(0.0, 0.0, 0.0);
                               q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.0, 0.0, ω0))
spacecraft = SpacecraftModel(links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia,
                             initial_condition=ic, id=1)
add_sensor!(spacecraft, root, RateGyro(period_s=1.0))      # no noise and no bias: exact rates
estimate = RateEstimate()
estimator = GyroRateEstimator(1, 1, estimate)
damper = EstimatedRateDamper(estimate=estimate)
mission_s = 120.0
planet = make_no_gram_planet(:earth)
estimator_configuration(mission_time_s) = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime, mission_time=mission_time_s, number_of_orbits=1,
        keplerian=false, orientation_sim=true, num_steps_to_save=120, data_rate=1.0),
    environment_model=EnvironmentModel(
        planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false),
    dynamics_model=DynamicsModel([spacecraft], ()),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(estimator,), navigation_rates=[1.0]),
    control_model=ControlModel(control_effectors=(damper,), control_rates=[damper.period_s]),
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances())

config = estimator_configuration(mission_s)
solution = run_simulation(config; return_solution=true)

# With an exact gyro sampled at the control ticks, the controller sees the true
# rate at each tick, and the run matches example A's closed form.
inertia = root.inertia[3, 3]
ticks = floor(Int, mission_s / damper.period_s) - 1
expected = ω0 * (1 - damper.gain_n_m_s * damper.period_s / inertia)^ticks
final = solution.u[end].sc[1].ω[3]
println("body rate after $(mission_s) s with the estimate: $(final) rad/s, by hand $(expected) rad/s")
@test final ≈ expected rtol=1e-6
# end snippet

# snippet: conformance
@testset "Estimator and controller conformance" begin
    # Keep the student demonstration at 120 s. This separate short configuration
    # checkpoints at 6.5 s, between ticks, while retaining the same nonzero-spin pair.
    conformance_config = estimator_configuration(13.0)
    fixture = conformance_fixture(configuration=conformance_config)
    reports = (test_navigation_model(estimator; fixture=fixture),
               test_control_model(damper; fixture=fixture))
    for report in reports
        @test only(filter(check -> check.id === :checkpoint_round_trip,
                          report.checks)).status === :pass
    end
end
# end snippet
