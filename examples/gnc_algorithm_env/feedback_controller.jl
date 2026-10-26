# Worked example A of "Writing a GNC Algorithm" (docs/src/user/gnc_feedback.md):
# a feedback controller that damps a spacecraft's body rate. Each region between a
# `# snippet:` line and `# end snippet` appears on that page as it is here; the HPC
# and extensibility docs gate compares them.
using SpaceAGORA
using SpaceAGORA.Conformance: conformance_fixture, test_control_model
using StaticArrays
using Test

# snippet: controller
Base.@kwdef mutable struct RateDamper <: SpaceAGORA.AbstractControlEffectorModel
    spacecraft_index::Int = 1           # the spacecraft it acts on, its position in the run
    gain_n_m_s::Float64 = 2.0           # torque per unit body rate, N·m per rad/s
    max_torque_n_m::Float64 = 1.0       # torque limit on each body axis, N·m
    period_s::Float64 = 1.0             # the tick period the law assumes, s
    # Run state: the held command and the time of the tick that set it.
    torque_body::SVector{3, Float64} = zero(SVector{3, Float64})
    tick_s::Float64 = NaN
end

# Declarations: what the controller acts on and what it assumes.
SpaceAGORA.spacecraft_scope(c::RateDamper) = BoundSpacecraft((c.spacecraft_index,))
SpaceAGORA.gnc_period_s(c::RateDamper) = c.period_s
SpaceAGORA.held_command(c::RateDamper, spacecraft_index) = c.torque_body
SpaceAGORA.held_command_time(c::RateDamper, spacecraft_index) = c.tick_s
# While its spacecraft is inactive (after an impact or a touchdown), the engine calls
# this at each tick in place of the tick (`skip_inactive`): the off command.
function SpaceAGORA.release_command!(c::RateDamper, t, spacecraft_index)
    c.torque_body, c.tick_s = zero(SVector{3, Float64}), t
    return nothing
end

# At each run's start: fresh run state, so two runs of one configuration start equal.
function SpaceAGORA.initialize_control!(c::RateDamper, u, p, t)
    c.torque_body, c.tick_s = zero(SVector{3, Float64}), NaN
    return nothing
end

# The tick (shape T), once per control period: run the law, limit it, hold the command.
function SpaceAGORA.update_control!(c::RateDamper, u, p, t, spacecraft_index)
    ω_body = SVector{3, Float64}(u.sc[spacecraft_index].ω)
    c.torque_body = clamp.(-c.gain_n_m_s .* ω_body, -c.max_torque_n_m, c.max_torque_n_m)
    c.tick_s = t
    return nothing
end

# The output (shape R), at every evaluation of the equations of motion: the held
# command, never the law. The engine calls it for the spacecraft in the scope
# only; direct callers may ask about others, so it checks the index.
function SpaceAGORA.control_wrench(c::RateDamper, x::StateSample, p, t)
    zero3 = zero(SVector{3, Float64})
    x.spacecraft_index == c.spacecraft_index || return zero3, zero3
    return zero3, c.torque_body             # force_ii (N), torque_body (N·m)
end

# A checkpoint must carry the held output between ticks, including its timestamp.
# These values are immutable; restoring them cannot share mutable buffers.
SpaceAGORA.checkpoint_state(c::RateDamper) =
    (torque_body=c.torque_body, tick_s=c.tick_s)
function SpaceAGORA.restore_checkpoint_state!(c::RateDamper, data)
    c.torque_body, c.tick_s = data.torque_body, data.tick_s
    return nothing
end

# Declared on a spacecraft: a copy for that spacecraft, with fresh run state.
SpaceAGORA.bind_spacecraft(c::RateDamper, spacecraft_index::Int) = RateDamper(;
    spacecraft_index, gain_n_m_s=c.gain_n_m_s, max_torque_n_m=c.max_torque_n_m, period_s=c.period_s)
# end snippet

# snippet: run
# A 500 kg cube far from any planet, spinning about body z; only the controller acts.
root = Link(root=true, m=500.0, dims=MVector(1.0, 1.0, 1.0))
ω0 = 0.05                                    # initial rate, rad/s
ic = CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), SVector(0.0, 0.0, 0.0);
                               q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.0, 0.0, ω0))
spacecraft = SpacecraftModel(links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia,
                             initial_condition=ic, id=1)
damper = RateDamper()
mission_s = 120.0
planet = make_no_gram_planet(:earth)
config = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime, mission_time=mission_s, number_of_orbits=1,
        keplerian=false, orientation_sim=true, num_steps_to_save=120, data_rate=1.0),
    environment_model=EnvironmentModel(
        planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false),
    dynamics_model=DynamicsModel([spacecraft], ()),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(damper,), control_rates=[damper.period_s]),
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances())

solution = run_simulation(config; return_solution=true)
# end snippet

# snippet: check
# By hand: the first tick is one period after the start, and each held torque
# -k ω(t_n) acts for one period on a cube (no gyroscopic torque), so after the
# ticks at 1, 2, ..., N s the rate is ω0 (1 - k Δt / I)^N.
inertia = root.inertia[3, 3]
ticks = floor(Int, mission_s / damper.period_s) - 1        # ticks before the end
expected = ω0 * (1 - damper.gain_n_m_s * damper.period_s / inertia)^ticks
final = solution.u[end].sc[1].ω[3]
println("body rate after $(mission_s) s: $(final) rad/s, by hand $(expected) rad/s")
@test final ≈ expected rtol=1e-6
# end snippet

# snippet: conformance
# The engine's own checks on the controller: hooks, reproducibility, copies, the
# reset at each run's start, output-hook purity and the copies `bind_spacecraft` makes.
@testset "RateDamper conformance" begin
    # Use the spinning case above: a zero-rate fixture would not expose a lost torque.
    report = test_control_model(damper;
        fixture=conformance_fixture(configuration=config))
    # Require the checkpoint check itself to pass on this nonzero-spin fixture.
    @test only(filter(check -> check.id === :checkpoint_round_trip,
                      report.checks)).status === :pass
end
# end snippet
