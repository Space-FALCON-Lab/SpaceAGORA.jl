module IntegratorDeprecationTests

# Engine callbacks must not call deprecated integrator functions. Under
# --depwarn=error, as CI's depwarn smoke runs and as users who forbid deprecations
# run, the warning becomes an error and ends the simulation. A child process runs,
# with --depwarn=error, one external-propagation case (the owner's sync callback)
# and one two-segment mission (the segment entry and transition callbacks).

using Test

const CHILD = raw"""
using SpaceAGORA, StaticArrays
using SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel
const EP = SpaceAGORA.ExternalPropagationServices
const MU = 3.986004418e14
const RADIUS = 7.0e6
const RATE = sqrt(MU / RADIUS^3)
orbit_pos(t) = RADIUS * SVector(cos(RATE * t), sin(RATE * t), 0.0)
orbit_vel(t) = RADIUS * RATE * SVector(-sin(RATE * t), cos(RATE * t), 0.0)

function configuration()
    body = SM.Link(root=true, m=10.0, dims=MVector(1.0, 1.0, 1.0), r=MVector(0.0, 0.0, 0.0),
        q=MVector(0.0, 0.0, 0.0, 1.0))
    sc = SM.SpacecraftModel(; joints=SM.Joint[], links=[body], root=body, inertia_tensor=body.inertia,
        initial_condition=SM.CartesianInitialCondition(orbit_pos(0.0), orbit_vel(0.0)))
    return make_example_config(planet=SM.make_no_gram_planet(:earth), spacecraft=sc, mission_time=4.0,
        initial_time=SM.InitialTime(year=2014, month=5, day=27), dynamic_effectors=(SM.InverseSquaredGravityModel(),),
        density_model=SM.NoAtmosphereModel(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=false,
        keplerian=true, verbose=false, results=false, results_directory=mktempdir(),
        solver_config=SM.SolverConfig(solver_mode=:dp8))
end

# An owner on the analytic circle, synced every 0.5 s: its sync callback writes the shadow entry.
struct CircleOwner <: EP.AbstractExternalPropagator end
mutable struct CircleRuntime
    n::Int
end
EP.external_spacecraft(::CircleOwner) = [1]
EP.external_step(::CircleOwner) = 0.5
EP.external_prepare(::CircleOwner, args, u0) = CircleRuntime(0)
function EP.external_sync!(rt::CircleRuntime, t::Float64)
    target = floor(Int, t / 0.5 + 1e-6)
    steps = target - rt.n
    rt.n = max(rt.n, target)
    return max(steps, 0)
end
EP.external_time(rt::CircleRuntime) = rt.n * 0.5
EP.external_state(rt::CircleRuntime, t::Float64, slot::Int) =
    (pos=orbit_pos(t), vel=orbit_vel(t), q=SVector(0.0, 0.0, 0.0, 1.0), ω=SVector(0.0, 0.0, 0.0))
EP.external_acceleration(rt::CircleRuntime, t::Float64, slot::Int) = -RATE^2 * orbit_pos(t)

owned = EP.with_external_propagators(configuration(), (CircleOwner(),))
result = SpaceAGORA.run_simulation(owned; return_results=true)
println("EXTERNAL_PROPAGATION_OK rows=", size(result.table, 1))

# Two segments: the entry callback and the transition both reach the solver's discontinuity flag.
plan = MissionPlan("two segments", configuration;
    segments=(MissionSegment(:first; until=MissionDuration(2.0)), MissionSegment(:second)))
mission = run_mission(plan; output_dir=mktempdir(), save_scene=false, save_payload=false)
println("MISSION_SEGMENTS_OK transitions=", length(mission.history.transitions))
"""

@testset "External propagation and mission segments run under --depwarn=error" begin
    command = `$(Base.julia_cmd()) --startup-file=no --depwarn=error --project=$(Base.active_project()) -e $(CHILD)`
    output = IOBuffer()
    ok = success(pipeline(ignorestatus(command); stdout=output, stderr=output))
    text = String(take!(output))
    ok || println(text)
    @test ok
    @test occursin("EXTERNAL_PROPAGATION_OK", text)
    @test occursin("MISSION_SEGMENTS_OK transitions=2", text)
end

end # module IntegratorDeprecationTests
