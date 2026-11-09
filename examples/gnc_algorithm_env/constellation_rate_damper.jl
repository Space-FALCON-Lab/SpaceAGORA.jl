# Rate damping of one target spacecraft in a constellation of N spacecraft, run
# together in one simulation. A RateDamper (worked example A, its controller code
# taken from `feedback_controller.jl`) is bound to the target and damps its body rate; the
# other N - 1 spacecraft must receive no torque from it. Every spacecraft is a
# 500 kg, 1 m cube spinning about body z, so each is checked by hand: the target's
# rate decays as ω0 (1 - k Δt / I)^n, the others keep ω0. The script also checks
# that binding the damper to an index outside 1:N is refused, and runs the
# conformance checks on the N-spacecraft configuration.
#
# From the SpaceAGORA checkout, prepare the environment once, then run the script
# with two or more threads so that the conformance checks also compare concurrent
# runs with sequential ones (with one thread that check is skipped):
#
#   julia --startup-file=no --compiled-modules=existing --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/setup.jl
#   julia --startup-file=no --compiled-modules=existing --depwarn=error --threads=2 --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/constellation_rate_damper.jl
#
# Every setting below can be edited here or overridden from the shell, e.g.
#
#   N_SPACECRAFT=10 TARGET=7 OUT_OF_RANGE=12 julia --startup-file=no --compiled-modules=existing --depwarn=error --threads=2 --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/constellation_rate_damper.jl
using SpaceAGORA
using SpaceAGORA.Conformance: conformance_fixture, test_control_model
using Printf
using StaticArrays
using Test

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
setting(name, default) = parse(Int, get(ENV, name, string(default)))

const N_SPACECRAFT = setting("N_SPACECRAFT", 3)                  # spacecraft in the run
const TARGET = setting("TARGET", 1)                              # run position the damper is bound to, 1:N
const OUT_OF_RANGE = setting("OUT_OF_RANGE", N_SPACECRAFT + 1)   # an index outside 1:N, for the refusal
const WALKER_PLANES = setting("WALKER_PLANES", 1)                # Walker planes, must divide N
const WALKER_PHASING = setting("WALKER_PHASING", 0)              # Walker phasing, 0:(planes - 1)
const RUN_CONFORMANCE = setting("RUN_CONFORMANCE", 1) == 1       # conformance checks (~20 s)
const MISSION_S = 120.0
const ω0 = 0.05                                                  # initial rate about body z, rad/s

N_SPACECRAFT >= 1 || error("N_SPACECRAFT = $N_SPACECRAFT: the run needs at least one spacecraft.")
1 <= TARGET <= N_SPACECRAFT ||
    error("TARGET = $TARGET is not a run position of $N_SPACECRAFT spacecraft; choose 1 to $N_SPACECRAFT. " *
          "(To test an out-of-range binding, set OUT_OF_RANGE instead.)")
1 <= OUT_OF_RANGE <= N_SPACECRAFT &&
    error("OUT_OF_RANGE = $OUT_OF_RANGE is a valid run position of $N_SPACECRAFT spacecraft, so the refusal " *
          "would not be tested; choose an index below 1 or above $N_SPACECRAFT (default $(N_SPACECRAFT + 1)).")
N_SPACECRAFT % WALKER_PLANES == 0 ||
    error("WALKER_PLANES = $WALKER_PLANES does not divide N_SPACECRAFT = $N_SPACECRAFT.")

# =============================================================================
# The controller from feedback_controller.jl
# =============================================================================

# Only the controller is taken from `feedback_controller.jl`: the region between its
# `# snippet: controller` and `# end snippet` lines (the RateDamper type and its
# hooks), evaluated here. The rest of that file (its own run and checks) is not run.
function include_snippet(path, name)
    lines = readlines(path)
    first_line = findfirst(==("# snippet: $name"), lines)
    first_line === nothing && error("$path has no `# snippet: $name` region.")
    last_line = findnext(==("# end snippet"), lines, first_line)
    return include_string(@__MODULE__, join(lines[first_line+1:last_line-1], "\n"), path)
end
include_snippet(joinpath(@__DIR__, "feedback_controller.jl"), "controller")

# A test instrument, not part of the controller: it wraps any control model, passes
# every hook through unchanged, and logs each tick as the engine delivered it (the
# time, the position the engine ticked the model for, the rate the law read and the
# command the model then holds, read through the `held_command` declaration).
struct TickRecord
    t::Float64
    spacecraft_index::Int
    ω_body::SVector{3, Float64}
    torque_body::SVector{3, Float64}
end

mutable struct TickRecorder{M} <: SpaceAGORA.AbstractControlEffectorModel
    inner::M
    ticks::Vector{TickRecord}
end
TickRecorder(inner) = TickRecorder(inner, TickRecord[])

SpaceAGORA.spacecraft_scope(r::TickRecorder) = spacecraft_scope(r.inner)
SpaceAGORA.gnc_period_s(r::TickRecorder) = gnc_period_s(r.inner)
SpaceAGORA.skip_inactive(r::TickRecorder) = skip_inactive(r.inner)
SpaceAGORA.held_command(r::TickRecorder, i) = held_command(r.inner, i)
SpaceAGORA.held_command_time(r::TickRecorder, i) = held_command_time(r.inner, i)
SpaceAGORA.release_command!(r::TickRecorder, t, i) = release_command!(r.inner, t, i)
SpaceAGORA.control_wrench(r::TickRecorder, x::StateSample, p, t) = control_wrench(r.inner, x, p, t)
SpaceAGORA.bind_spacecraft(r::TickRecorder, i::Int) = TickRecorder(SpaceAGORA.bind_spacecraft(r.inner, i))

function SpaceAGORA.initialize_control!(r::TickRecorder, u, p, t)
    SpaceAGORA.initialize_control!(r.inner, u, p, t)
    empty!(r.ticks)
    return nothing
end

function SpaceAGORA.update_control!(r::TickRecorder, u, p, t, i)
    SpaceAGORA.update_control!(r.inner, u, p, t, i)
    push!(r.ticks, TickRecord(t, i, SVector{3, Float64}(u.sc[i].ω), held_command(r.inner, i)))
    return nothing
end

# The inner model's checkpoint plus the log, as plain arrays.
function SpaceAGORA.checkpoint_state(r::TickRecorder)
    return (inner=SpaceAGORA.checkpoint_state(r.inner),
            tick_t=Float64[x.t for x in r.ticks], tick_index=Int[x.spacecraft_index for x in r.ticks],
            tick_ω=SVector{3, Float64}[x.ω_body for x in r.ticks],
            tick_torque=SVector{3, Float64}[x.torque_body for x in r.ticks])
end
function SpaceAGORA.restore_checkpoint_state!(r::TickRecorder, data)
    SpaceAGORA.restore_checkpoint_state!(r.inner, data.inner)
    r.ticks = TickRecord.(data.tick_t, data.tick_index, data.tick_ω, data.tick_torque)
    return nothing
end

# =============================================================================
# The configuration
# =============================================================================

planet = make_no_gram_planet(:earth)

# Example A's spacecraft for every position: a 500 kg, 1 m cube spinning at ω0
# about body z. A cube has no gyroscopic torque, so without a torque its rate stays
# ω0, and under the damper it follows example A's hand calculation.
function factory(plane, slot, index)
    root = Link(root=true, m=500.0, dims=MVector(1.0, 1.0, 1.0))
    # The builder replaces the orbit and keeps this attitude and body rate.
    ic = CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(0.0, 7.5e3, 0.0);
                                   q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.0, 0.0, ω0))
    return SpacecraftModel(links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia,
                           initial_condition=ic)
end

# Walker delta i: total/planes/phasing (docs/src/user/constellations.md) at 510 km.
# The builder numbers the spacecraft plane by plane and slot by slot (the run
# positions) and gives them ids 1 to N.
const WALKER = (total=N_SPACECRAFT, planes=WALKER_PLANES, phasing=WALKER_PHASING, altitude=510e3, i=35)

# A configuration of fresh spacecraft with the given control effectors. Gravity
# (J2) acts on every spacecraft and applies no torque (gravity gradient is a
# separate effector, not included).
function constellation_configuration(control_effectors::Tuple)
    spacecraft = walker_constellation(factory; WALKER..., planet=planet)
    rates = Float64[gnc_period_s(m) for m in control_effectors]
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, mission_time=MISSION_S, number_of_orbits=1,
            keplerian=false, orientation_sim=true, num_steps_to_save=120, data_rate=1.0),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=rates),
        initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances())
end

damper = RateDamper(spacecraft_index=TARGET)
recorder = TickRecorder(damper)
config_on = constellation_configuration((recorder,))
config_off = constellation_configuration(())

println("== Configuration")
println("  Walker $(WALKER.i)°: $(WALKER.total)/$(WALKER.planes)/$(WALKER.phasing) at $(WALKER.altitude / 1e3) km; ",
        "every spacecraft a 500 kg 1 m cube, ω0 = $ω0 rad/s about z")
println("  RateDamper bound to position $TARGET: scope = ", spacecraft_scope(damper), ", gain = ",
        damper.gain_n_m_s, " N·m·s, limit = ", damper.max_torque_n_m, " N·m, period = ", damper.period_s, " s")

# =============================================================================
# Runs: with the damper, and with it removed
# =============================================================================

# `isolate_state=false` runs `recorder` and `damper` themselves, so the tick log and
# the final held command are readable afterwards.
sol_on = run_simulation(config_on; return_solution=true, isolate_state=false)
sol_off = run_simulation(config_off; return_solution=true)

# =============================================================================
# Recorded output
# =============================================================================

println("\n== Ticks and commands")
ticks = recorder.ticks
tick_positions = sort(unique(r.spacecraft_index for r in ticks))
println("  ", length(ticks), " ticks, t = ", first(ticks).t, " to ", last(ticks).t,
        " s, all for run position(s) ", tick_positions)
println("     t (s)   ω_z (rad/s)   torque_z (N·m)")
for r in ticks[[1:3; length(ticks)-1:length(ticks)]]
    @printf("  %7.1f   % .4e   % .4e\n", r.t, r.ω_body[3], r.torque_body[3])
end

# By hand (example A): the first tick is one period after the start, and each held
# torque -k ω(tₙ) acts for one period, so after the ticks at 1, 2, ..., N s the
# target's rate is ω0 (1 - k Δt / I)^N. Every other spacecraft keeps ω0.
inertia = config_on.dynamics_model.spacecraft[TARGET].inertia_tensor[3, 3]
n_ticks = floor(Int, MISSION_S / damper.period_s) - 1
damped = ω0 * (1 - damper.gain_n_m_s * damper.period_s / inertia)^n_ticks
by_hand(i) = i == TARGET ? damped : ω0
ω_on(i) = SVector{3}(sol_on.u[end].sc[i].ω)
ω_off(i) = SVector{3}(sol_off.u[end].sc[i].ω)

println("\n== Body rate after $(MISSION_S) s, against by hand calculation")
println("  pos   ω_z with damper    by hand            ω_z without damper  |ω_x|+|ω_y| with damper")
for i in 1:N_SPACECRAFT
    @printf("  %3d   %.10e   %.10e   %.10e    %.1e%s\n", i, ω_on(i)[3], by_hand(i), ω_off(i)[3],
            abs(ω_on(i)[1]) + abs(ω_on(i)[2]), i == TARGET ? "   <- target" : "")
end

println("\n== The output hook asked directly about each spacecraft (final state)")
wrenches = map(1:N_SPACECRAFT) do i
    s = sol_on.u[end].sc[i]
    x = StateSample(SVector{3}(s.pos), SVector{3}(s.vel), Float64(s.mass), SVector{4}(s.q), SVector{3}(s.ω),
                    config_on.dynamics_model.spacecraft[i], i)
    control_wrench(damper, x, nothing, MISSION_S)
end
for (i, (f, τ)) in enumerate(wrenches)
    @printf("  pos %3d: force_ii = %s N, torque_body = (% .3e, % .3e, % .3e) N·m\n", i, f, τ...)
end

@testset "RateDamper bound to position $TARGET of $N_SPACECRAFT" begin
    # Ticks: only the target, once per period from 1 s to the last period before the end.
    @test tick_positions == [TARGET]
    @test [r.t for r in ticks] ≈ collect(1.0:damper.period_s:(MISSION_S - damper.period_s))
    @test all(r.torque_body ≈ clamp.(-damper.gain_n_m_s .* r.ω_body, -1, 1) for r in ticks)
    # Every spacecraft against the hand calculation, with the damper; without it, all keep ω0.
    for i in 1:N_SPACECRAFT
        @test ω_on(i)[3] ≈ by_hand(i) rtol=1e-6
        @test abs(ω_on(i)[1]) + abs(ω_on(i)[2]) < 1e-12
        @test ω_off(i)[3] ≈ ω0 rtol=1e-9
    end
    # Direct queries: a torque for the target only.
    @test all(iszero(wrenches[i][2]) for i in setdiff(1:N_SPACECRAFT, TARGET))
    @test !iszero(wrenches[TARGET][2])
end

# =============================================================================
# Refusal: a scope outside the run
# =============================================================================

println("\n== RateDamper(spacecraft_index=$OUT_OF_RANGE) in a run of $N_SPACECRAFT spacecraft")
stray = RateDamper(spacecraft_index=OUT_OF_RANGE, tick_s=-1.0)    # -1: a sentinel that initialize_control! would clear
refusal = try
    run_simulation(constellation_configuration((stray,)); return_solution=true, isolate_state=false)
    nothing
catch err
    err
end
if refusal === nothing
    println("  NOT refused: the run completed")
else
    println("  refused with ", typeof(refusal), ":")
    println("    ", sprint(showerror, refusal))
    println("  refused before the run started (initialize_control! never called): ", stray.tick_s == -1.0)
end
@testset "Out-of-range binding is refused" begin
    @test refusal isa ArgumentError
    message = refusal === nothing ? "" : sprint(showerror, refusal)
    @test occursin("RateDamper", message)                 # which model
    @test occursin("spacecraft $OUT_OF_RANGE", message)   # which index is wrong
    @test occursin("1 to $N_SPACECRAFT", message)         # the valid range
    @test stray.tick_s == -1.0                            # nothing ran
end

# =============================================================================
# Conformance on the N-spacecraft configuration
# =============================================================================

if RUN_CONFORMANCE
    println("\n== Conformance")
    Threads.nthreads() < 2 &&
        @warn "Julia runs with $(Threads.nthreads()) thread: concurrent_equals_sequential is skipped. Start Julia with --threads=2."
    @testset "RateDamper conformance, $N_SPACECRAFT spacecraft" begin
        # Example A's RateDamper alone, without the recorder, in the N-spacecraft run.
        # Includes the check that the output hooks return zero for spacecraft outside the scope.
        # positive_consumption is skipped: the damper uses no propellant.
        bare = RateDamper(spacecraft_index=TARGET)
        report = test_control_model(bare; fixture=conformance_fixture(configuration=constellation_configuration((bare,))))
        for check in report.checks
            println("  ", rpad(String(check.id), 40), check.status)
        end
        status(id) = only(filter(check -> check.id === id, report.checks)).status
        # With two or more threads, the concurrent check must run and pass.
        Threads.nthreads() >= 2 && @test status(:concurrent_equals_sequential) === :pass
        # The round trip covers the RateDamper's whole run state: its held torque and the
        # time of the tick that set it. A pass is not evidence that a controller with more
        # run state (a plan, an estimate, a solver's warm start) is serialized; each such
        # state needs its own checkpoint_state entry and a test whose result depends on it.
        @test status(:checkpoint_round_trip) === :pass
    end
else
    println("\n== Conformance skipped (RUN_CONFORMANCE=0)")
end
