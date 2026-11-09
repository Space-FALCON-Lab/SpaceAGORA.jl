# Estimated rate damping of one target spacecraft in a constellation of N
# spacecraft, run together in one simulation, and the pair's run state across a
# checkpoint. Worked example B (its link, estimator and consumer code taken from
# `estimator.jl`) is bound to the target: a gyro estimator writes the rate estimate
# into a shared link and a rate damper reads it. Every spacecraft is a 500 kg, 1 m
# cube spinning about body z and only the target carries a gyro, so each is checked
# by hand: with an exact gyro sampled at the control ticks the target's rate decays
# as ω0 (1 - k Δt / I)^n, the others keep ω0.
#
# The checkpoint checks write a real checkpoint file between two ticks and resume
# from it: once for a short leg with no tick, so the estimate and the held command
# the resumed run holds can only have come from the file, and once to the end of
# the run, compared with a run checkpointed at the same times and never interrupted.
# The script ends with the conformance checks of both models on the N-spacecraft
# configuration.
#
# From the SpaceAGORA checkout, prepare the environment once, then run the script
# with two or more threads so that the conformance checks also compare concurrent
# runs with sequential ones (with one thread that check is skipped):
#
#   julia --startup-file=no --compiled-modules=existing --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/setup.jl
#   julia --startup-file=no --compiled-modules=existing --depwarn=error --threads=2 --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/constellation_estimator.jl
#
# Every setting below can be edited here or overridden from the shell, e.g.
#
#   N_SPACECRAFT=10 TARGET=7 CHECKPOINT_S=40.5 julia --startup-file=no --compiled-modules=existing --depwarn=error --threads=2 --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/constellation_estimator.jl
using SpaceAGORA
using SpaceAGORA.Conformance: conformance_fixture, test_control_model, test_navigation_model
using Printf
using StaticArrays
using Test

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
setting(name, default) = parse(Int, get(ENV, name, string(default)))
setting_s(name, default) = parse(Float64, get(ENV, name, string(default)))

const N_SPACECRAFT = setting("N_SPACECRAFT", 3)                  # spacecraft in the run
const TARGET = setting("TARGET", 1)                              # run position the pair is bound to, 1:N
const WALKER_PLANES = setting("WALKER_PLANES", 1)                # Walker planes, must divide N
const WALKER_PHASING = setting("WALKER_PHASING", 0)              # Walker phasing, 0:(planes - 1)
const RUN_CONFORMANCE = setting("RUN_CONFORMANCE", 1) == 1       # conformance checks (~30 s)
const MISSION_S = 120.0
const PERIOD_S = 1.0                                             # gyro, estimator and damper period, s
const CHECKPOINT_S = setting_s("CHECKPOINT_S", 60.5)             # checkpoint time, between two ticks
const ω0 = 0.05                                                  # initial rate about body z, rad/s

N_SPACECRAFT >= 1 || error("N_SPACECRAFT = $N_SPACECRAFT: the run needs at least one spacecraft.")
1 <= TARGET <= N_SPACECRAFT ||
    error("TARGET = $TARGET is not a run position of $N_SPACECRAFT spacecraft; choose 1 to $N_SPACECRAFT.")
N_SPACECRAFT % WALKER_PLANES == 0 ||
    error("WALKER_PLANES = $WALKER_PLANES does not divide N_SPACECRAFT = $N_SPACECRAFT.")
# Between ticks, so the resumed run must carry the estimate and the held command to
# the next tick; and at least one period before the end, so a tick follows it.
(PERIOD_S < CHECKPOINT_S < MISSION_S - PERIOD_S && !isinteger(CHECKPOINT_S / PERIOD_S)) ||
    error("CHECKPOINT_S = $CHECKPOINT_S must fall between two ticks (not a multiple of $PERIOD_S s), " *
          "after the first and at least one period before the end ($MISSION_S s), e.g. 60.5.")

# =============================================================================
# The link, estimator and consumer from estimator.jl
# =============================================================================

# Only the models are taken from `estimator.jl`: the regions between its
# `# snippet: link`, `# snippet: estimator` and `# snippet: consumer` lines and
# `# end snippet`, evaluated here. The rest of that file (its own run and checks)
# is not run.
function include_snippet(path, name)
    lines = readlines(path)
    first_line = findfirst(==("# snippet: $name"), lines)
    first_line === nothing && error("$path has no `# snippet: $name` region.")
    last_line = findnext(==("# end snippet"), lines, first_line)
    return include_string(@__MODULE__, join(lines[first_line+1:last_line-1], "\n"), path)
end
for name in ("link", "estimator", "consumer")
    include_snippet(joinpath(@__DIR__, "estimator.jl"), name)
end

# =============================================================================
# The configuration
# =============================================================================

planet = make_no_gram_planet(:earth)

# Example B's spacecraft for every position: a 500 kg, 1 m cube spinning at ω0
# about body z. Only the target carries the gyro (no noise and no bias: exact
# rates), its first sensor, which the estimator reads.
function factory(plane, slot, index)
    root = Link(root=true, m=500.0, dims=MVector(1.0, 1.0, 1.0))
    # The builder replaces the orbit and keeps this attitude and body rate.
    ic = CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(0.0, 7.5e3, 0.0);
                                   q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.0, 0.0, ω0))
    spacecraft = SpacecraftModel(links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia,
                                 initial_condition=ic)
    index == TARGET && add_sensor!(spacecraft, root, RateGyro(period_s=PERIOD_S))
    return spacecraft
end

# Walker delta i: total/planes/phasing (docs/src/user/constellations.md) at 510 km.
# The builder numbers the spacecraft plane by plane and slot by slot (the run
# positions) and gives them ids 1 to N.
const WALKER = (total=N_SPACECRAFT, planes=WALKER_PLANES, phasing=WALKER_PHASING, altitude=510e3, i=35)

# A configuration of fresh spacecraft and, with `pair=true`, a fresh estimator and
# damper bound to the target at configuration level, sharing one new estimate.
# `settings` go to SimulationSettings (the checkpoint fields). Gravity (J2) acts on
# every spacecraft and applies no torque (gravity gradient is a separate effector,
# not included).
function constellation_configuration(; mission_s=MISSION_S, pair=true, settings...)
    spacecraft = walker_constellation(factory; WALKER..., planet=planet)
    estimate = RateEstimate()
    navigation = pair ? (GyroRateEstimator(TARGET, 1, estimate),) : ()
    control = pair ? (EstimatedRateDamper(spacecraft_index=TARGET, period_s=PERIOD_S, estimate=estimate),) : ()
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(;
            results=false, verbose=false, generate_plots=false, normalize=false, settings...),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, mission_time=mission_s, number_of_orbits=1,
            keplerian=false, orientation_sim=true, num_steps_to_save=120, data_rate=1.0),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=navigation, navigation_rates=fill(PERIOD_S, length(navigation))),
        control_model=ControlModel(control_effectors=control, control_rates=fill(PERIOD_S, length(control))),
        initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances())
end

# The estimator and damper of the configuration that ran (`return_results=true`
# returns it, with the models' state at the end of the run).
ran_pair(res) = (estimator=only(res.configuration.navigation_model.navigation_effectors),
                 damper=only(res.configuration.control_model.control_effectors))

println("== Configuration")
println("  Walker $(WALKER.i)°: $(WALKER.total)/$(WALKER.planes)/$(WALKER.phasing) at $(WALKER.altitude / 1e3) km; ",
        "every spacecraft a 500 kg 1 m cube, ω0 = $ω0 rad/s about z; gyro on position $TARGET only")
println("  GyroRateEstimator and EstimatedRateDamper bound to position $TARGET, linked by one RateEstimate; ",
        "period $PERIOD_S s, gain 2.0 N·m·s, limit 1.0 N·m")

# =============================================================================
# Runs: with the pair, and with it removed
# =============================================================================

res_on = run_simulation(constellation_configuration(); return_results=true, return_solution=true)
res_off = run_simulation(constellation_configuration(pair=false); return_solution=true)
on = ran_pair(res_on)

# =============================================================================
# Body rates and the link
# =============================================================================

# By hand (example A): the first tick is one period after the start, and each held
# torque -k ω(tₙ) acts for one period, so after the ticks at 1, 2, ..., n s the
# target's rate is ω0 (1 - k Δt / I)^n. Every other spacecraft keeps ω0. The gyro
# samples before the navigation tick, and navigation ticks before control at the
# same time, so the damper reads the exact rate at each tick.
inertia = res_on.configuration.dynamics_model.spacecraft[TARGET].inertia_tensor[3, 3]
decay = 1 - on.damper.gain_n_m_s * PERIOD_S / inertia
rate_after(n) = ω0 * decay^n                    # the target's rate after n ticks
n_ticks = floor(Int, MISSION_S / PERIOD_S) - 1
by_hand(i) = i == TARGET ? rate_after(n_ticks) : ω0
ω_on(i) = SVector{3}(res_on.solution.u[end].sc[i].ω)
ω_off(i) = SVector{3}(res_off.u[end].sc[i].ω)

println("\n== Body rate after $(MISSION_S) s, against the hand calculation")
println("  pos   ω_z with pair      by hand            ω_z without pair    |ω_x|+|ω_y| with pair")
for i in 1:N_SPACECRAFT
    @printf("  %3d   %.10e   %.10e   %.10e    %.1e%s\n", i, ω_on(i)[3], by_hand(i), ω_off(i)[3],
            abs(ω_on(i)[1]) + abs(ω_on(i)[2]), i == TARGET ? "   <- target" : "")
end

# At the end of the run the estimate holds the gyro's measurement at the last tick,
# the rate after one tick fewer, and the damper holds -k times it.
last_tick = n_ticks * PERIOD_S
println("\n== The link at the end of the run")
@printf("  estimate: ω_z = %.10e rad/s at t = %.1f s (by hand %.10e at %.1f s)\n",
        on.estimator.estimate.ω_body[3], on.estimator.estimate.time_s, rate_after(n_ticks - 1), last_tick)
@printf("  damper:   torque_z = %.10e N·m, held since t = %.1f s; reads the estimator's object: %s\n",
        on.damper.torque_body[3], on.damper.tick_s, on.damper.estimate === on.estimator.estimate)

@testset "Estimator and damper bound to position $TARGET of $N_SPACECRAFT" begin
    # Every spacecraft against the hand calculation, with the pair; without it, all keep ω0.
    for i in 1:N_SPACECRAFT
        @test ω_on(i)[3] ≈ by_hand(i) rtol=1e-6
        @test abs(ω_on(i)[1]) + abs(ω_on(i)[2]) < 1e-12
        @test ω_off(i)[3] ≈ ω0 rtol=1e-9
    end
    # The link: one object, written at the last tick, read by the damper.
    @test on.damper.estimate === on.estimator.estimate
    @test on.estimator.estimate.time_s == last_tick
    @test on.estimator.estimate.ω_body[3] ≈ rate_after(n_ticks - 1) rtol=1e-6
    @test on.damper.tick_s == last_tick
    @test on.damper.torque_body ≈ -on.damper.gain_n_m_s .* on.estimator.estimate.ω_body
end

# =============================================================================
# Run state across a checkpoint
# =============================================================================

# Each run state the pair keeps, and who saves it (estimator.jl): the estimate (its
# rate and measurement time) is the shared link's, saved by its producer, the
# estimator, and restored into the existing link so the damper still reads the same
# object; the held torque and its tick time are the damper's own. A run with
# checkpoints is solved in segments of `checkpoint_interval_s`; one that ends at a
# checkpoint can be resumed to a later mission time (docs/src/user/outputs.md,
# Checkpointing).
checkpointed(dir; mission_s=MISSION_S, resume=false) = constellation_configuration(; mission_s,
    checkpoint_enabled=true, checkpoint_interval_s=CHECKPOINT_S, checkpoint_directory=dir,
    resume_from_checkpoint=resume)

# The first leg: from the start to CHECKPOINT_S, which writes the checkpoint. Its
# rows are kept beside the checkpoint (`return_results=true`), so a resumed run's
# results hold the whole run.
leg_dir = mktempdir()
leg = run_simulation(checkpointed(leg_dir; mission_s=CHECKPOINT_S); return_results=true)
saved = ran_pair(leg)
last_tick_before = floor(CHECKPOINT_S / PERIOD_S) * PERIOD_S

# A short leg from the checkpoint to a time before the next tick: nothing ticks,
# so the estimate and the held command at its end came from the file. A pair
# without these restores would end it with the link's reset (zero rate, NaN time)
# and the damper's reset (zero torque, NaN time).
short_dir = mktempdir()
cp(leg_dir, short_dir; force=true)
short = ran_pair(run_simulation(checkpointed(short_dir; mission_s=(CHECKPOINT_S + last_tick_before + PERIOD_S) / 2,
                                             resume=true); return_results=true))

# To the end: the resumed run against a run checkpointed at the same times and never
# interrupted. Both keep every saved row, so whole histories are compared.
resumed = run_simulation(checkpointed(leg_dir; resume=true); return_results=true)
uninterrupted = run_simulation(checkpointed(mktempdir()); return_results=true)
numeric(table) = Matrix{Float64}(table[:, [c for c in names(table) if eltype(table[!, c]) <: Real]])
rows_resumed, rows_uninterrupted = numeric(resumed.table), numeric(uninterrupted.table)
same_shape = size(rows_resumed) == size(rows_uninterrupted)
largest_difference = same_shape ? maximum(abs.(rows_resumed .- rows_uninterrupted) ./ max.(abs.(rows_uninterrupted), 1.0)) : Inf
r, u = ran_pair(resumed), ran_pair(uninterrupted)

println("\n== Run state across a checkpoint at $(CHECKPOINT_S) s")
@printf("  saved:    estimate ω_z = %.10e rad/s at %.1f s, torque_z = %.10e N·m at %.1f s (by hand %.10e rad/s)\n",
        saved.estimator.estimate.ω_body[3], saved.estimator.estimate.time_s, saved.damper.torque_body[3],
        saved.damper.tick_s, rate_after(Int(last_tick_before / PERIOD_S) - 1))
@printf("  restored: estimate ω_z = %.10e rad/s at %.1f s, torque_z = %.10e N·m at %.1f s (short leg, no tick)\n",
        short.estimator.estimate.ω_body[3], short.estimator.estimate.time_s, short.damper.torque_body[3], short.damper.tick_s)
println("  resumed to $(MISSION_S) s against uninterrupted: ", size(rows_resumed, 1), " and ", size(rows_uninterrupted, 1),
        " rows, every saved value equal: ", rows_resumed == rows_uninterrupted,
        @sprintf(" (largest relative difference %.1e)", largest_difference))
println("  final pair state equal: ", r.estimator.estimate.ω_body == u.estimator.estimate.ω_body &&
        r.estimator.estimate.time_s == u.estimator.estimate.time_s && r.damper.torque_body == u.damper.torque_body &&
        r.damper.tick_s == u.damper.tick_s)

@testset "Pair run state across a checkpoint" begin
    # The checkpoint falls between ticks: the state saved is the last tick's.
    @test saved.estimator.estimate.time_s == last_tick_before
    @test saved.damper.tick_s == last_tick_before
    # Restored from the file, without a tick: the values saved, in one shared link.
    @test short.estimator.estimate.ω_body == saved.estimator.estimate.ω_body
    @test short.estimator.estimate.time_s == saved.estimator.estimate.time_s
    @test short.damper.torque_body == saved.damper.torque_body
    @test short.damper.tick_s == saved.damper.tick_s
    @test short.damper.estimate === short.estimator.estimate
    # Resumed to the end: the same rows and the same final pair state as never interrupted.
    @test same_shape && largest_difference <= 1e-12
    @test r.estimator.estimate.ω_body == u.estimator.estimate.ω_body && r.estimator.estimate.time_s == u.estimator.estimate.time_s
    @test r.damper.torque_body == u.damper.torque_body && r.damper.tick_s == u.damper.tick_s
end

# =============================================================================
# Conformance on the N-spacecraft configuration
# =============================================================================

if RUN_CONFORMANCE
    println("\n== Conformance")
    Threads.nthreads() < 2 &&
        @warn "Julia runs with $(Threads.nthreads()) thread: concurrent_equals_sequential is skipped. Start Julia with --threads=2."
    @testset "Estimator and damper conformance, $N_SPACECRAFT spacecraft" begin
        # As in estimator.jl: a separate short configuration whose checkpoint round trip
        # writes its checkpoint halfway, at 6.5 s, between two ticks. Skipped:
        # positive_consumption (the damper uses no propellant) and bind_fresh_state
        # (example B declares no bind_spacecraft; constellation_ensemble.jl adds it and
        # checks its copies).
        config = constellation_configuration(mission_s=13.0)
        fixture = conformance_fixture(configuration=config)
        estimator = only(config.navigation_model.navigation_effectors)
        damper = only(config.control_model.control_effectors)
        for report in (test_navigation_model(estimator; fixture=fixture), test_control_model(damper; fixture=fixture))
            println("  ", report.target, " (", report.kind, ")")
            for check in report.checks
                println("    ", rpad(String(check.id), 38), check.status)
            end
            status(id) = only(filter(check -> check.id === id, report.checks)).status
            Threads.nthreads() >= 2 && @test status(:concurrent_equals_sequential) === :pass
            @test status(:checkpoint_round_trip) === :pass
        end
        # The trajectory alone cannot show that the estimate was restored: after a resume
        # the navigation tick refreshes it before the control tick at the same time reads
        # it. The short leg above checks the restored values themselves.
    end
else
    println("\n== Conformance skipped (RUN_CONFORMANCE=0)")
end
