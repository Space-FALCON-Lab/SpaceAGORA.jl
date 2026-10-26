module SensorTests
using Test, StaticArrays, LinearAlgebra, Statistics, Random, SpaceAGORA
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel

# A navigation effector that keeps the readings of the first sensor of the
# spacecraft of the run at its ticks.
mutable struct ReadingLog <: AbstractNavigationModel
    readings::Vector{Vector{Float64}}
    times::Vector{Float64}
end
function SpaceAGORA.update_navigation!(n::ReadingLog, u, p, t::Float64, i::Int)
    sc = p.args.dynamics_model.spacecraft[i]
    push!(n.readings, collect(sensor_measurement(sc, 1))); push!(n.times, sensor_measurement_time(sc, 1))
    return nothing
end

# A user-defined sensor kind: it records the time of every sample and
# measures one random number from its stream.
mutable struct SampleClock <: AbstractSensor
    period_s::Float64
    seed::Int
    rng::Any
    latest::SVector{1, Float64}
    latest_time_s::Float64
    times::Vector{Float64}
end
SampleClock(period_s; seed=1) = SampleClock(period_s, seed, nothing, SVector(NaN), NaN, Float64[])
SpaceAGORA.sensor_dimension(::SampleClock) = 1
function SpaceAGORA.sensor_measure(s::SampleClock, truth::SensorTruth, ctx, rng)
    push!(s.times, truth.t_s)
    return SVector(rand(rng))
end

# A one-component sensor whose measurement is a number: the body's mass.
mutable struct MassGauge <: AbstractSensor
    period_s::Float64
    seed::Int
    rng::Any
    latest::Float64
    latest_time_s::Float64
end
MassGauge() = MassGauge(1.0, 1, nothing, NaN, NaN)
SpaceAGORA.sensor_dimension(::MassGauge) = 1
SpaceAGORA.sensor_measure(::MassGauge, truth::SensorTruth, ctx, rng) = truth.mass_kg

const ω0 = SVector(0.01, -0.02, 0.005)

# A controller that holds levels for the mounted actuators of spacecraft 1.
mutable struct HoldLevels <: AbstractControlEffectorModel
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::HoldLevels) = :mounted
SpaceAGORA.actuator_levels(c::HoldLevels, i::Int) = i == 1 ? c.levels : nothing
SpaceAGORA.update_control!(::HoldLevels, u, p, t::Float64, i::Int64) = nothing

# A gravitational model of the user's own: a uniform field along -z.
struct UniformGravity <: SpaceAGORA.AbstractForceTorqueModel end
SpaceAGORA.is_gravitational(::UniformGravity) = true
SpaceAGORA.wrench(::UniformGravity, x::SpaceAGORA.StateSample, env::SpaceAGORA.EnvironmentSample, t::Float64) =
    (x.mass_kg * SVector(0.0, 0.0, -0.5), SVector(0.0, 0.0, 0.0))

const EPOCH = SM.InitialTime(year=2024, month=3, day=20, hour=0, minute=0, second=0.0)
sun_hat() = (e = SM.SimpleEphemeridesModel(); SVector{3, Float64}(SM.ephemerides_sun_ii(e, SM.Earth(), SM.ephemerides_epoch_seconds(e, EPOCH))))

function sensor_spacecraft(; sensors=(RateGyro(period_s=0.5), StarTracker(period_s=1.0), Magnetometer(period_s=0.5)), panel_q=nothing,
                            direction=SVector(1.0, 0.0, 0.0), spin=ω0)
    root = SM.Link(root=true, m=50.0, dims=MVector(1.0, 1.0, 1.0))
    r = SM.Earth().Rp_e + 500e3
    d = normalize(SVector{3, Float64}(direction))
    along = normalize(cross(SVector(0.0, 0.0, 1.0), d))
    ic = SM.CartesianInitialCondition(r * d, sqrt(SM.Earth().μ / r) * along; q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=spin)
    links = [root]
    mount = root
    if panel_q !== nothing
        mount = SM.Link(m=1.0, r=MVector(0.0, 1.0, 0.0), q=MVector(panel_q...))
        push!(links, mount)
    end
    sc = SM.SpacecraftModel(links=links, root=root, inertia_tensor=root.inertia, initial_condition=ic, id=3)
    for s in sensors
        add_sensor!(sc, mount, s)
    end
    return sc
end

function sensor_args(sc; mission_time=10.0, navigation=(), settings=(;), solver=nothing, orientation_sim=true, effectors=(), control=())
    base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=mission_time,
        initial_time=EPOCH,
        dynamic_effectors=effectors, ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=orientation_sim,
        keplerian=false, EI_km=100.0, verbose=false, results=false,
        solver_config=solver === nothing ? nothing : SM.SolverConfig(solver_mode=solver))
    extra = isempty(navigation) ? (;) : (; navigation_model=SM.NavigationModel(navigation_effectors=navigation, navigation_rates=fill(1.0, length(navigation))))
    isempty(control) || (extra = merge(extra, (; control_model=SM.ControlModel(control_effectors=control, control_rates=fill(1.0, length(control))))))
    if !isempty(settings)
        st = base.simulation_settings
        fields = NamedTuple{fieldnames(typeof(st))}(map(n -> getfield(st, n), fieldnames(typeof(st))))
        extra = merge(extra, (; simulation_settings=SM.SimulationSettings(; merge(fields, settings)...)))
    end
    return SM.SimConfig._with_configuration(base; extra...)
end

sensor_run(sc; kw...) = run_simulation(sensor_args(sc; kw...); isolate_state=false, return_results=true, return_solution=true)

# The row of `table` at time `t` (a sample time that is also a saved time).
row_at(table, t) = findfirst(==(t), table.time)

@testset "Sensors" begin
    @testset "records and mounting" begin
        @test_throws ArgumentError RateGyro(period_s=0.0)
        @test_throws ArgumentError StarTracker(cross_boresight_noise_rad=-1.0)
        @test_throws ArgumentError Magnetometer(noise_t=NaN)
        root = SM.Link(root=true, m=10.0)
        sc = SM.SpacecraftModel(links=[root], root=root)
        @test_throws ArgumentError add_sensor!(sc, SM.Link(m=1.0), RateGyro())          # not one of its links
        g = add_sensor!(sc, root, RateGyro())
        @test_throws ArgumentError add_sensor!(sc, root, g)                               # mounted twice
        @test [e.sensor for e in mounted_sensors(sc)] == [g] && all(isnan, sensor_measurement(g))
        @test sensor_measurement(sc, 1) === sensor_measurement(g) && isnan(sensor_measurement_time(sc, 1))
        @test sensor_dimension(g) == 3 && sensor_dimension(StarTracker()) == 4 && sensor_requirements(Magnetometer()).magnetic_field
        @test sensor_period_s(g) == 0.1
        for bad in (0.0, -1.0, NaN, Inf)                                                # a user kind's period is checked
            @test_throws ArgumentError add_sensor!(sc, root, SampleClock(bad))
        end
    end

    @testset "the noise models follow their parameters" begin
        truth = SensorTruth(SVector(7e6, 0.0, 0.0), SVector(0.0, 7.5e3, 0.0), SVector(0.0, 0.0, 0.0, 1.0), ω0, 50.0, 0.0,
                            SVector(2e-5, -1e-5, 3e-5))
        ctx = actuator_context(SM.Link(root=true, m=1.0))
        rng = Xoshiro(11)
        gyro = RateGyro(period_s=0.1, angle_random_walk=1e-3, initial_bias=SVector(1e-3, 0.0, 0.0))
        errors = [sensor_measure(gyro, truth, ctx, rng) - ω0 for _ in 1:4000]
        σ = 1e-3 / sqrt(0.1)
        @test std(getindex.(errors, 2)) ≈ σ rtol=0.1
        @test mean(getindex.(errors, 1)) ≈ 1e-3 atol=4σ / sqrt(4000)                      # the bias
        walker = RateGyro(period_s=1.0, rate_random_walk=1e-4)
        steps = Float64[]
        for _ in 1:2000
            before = walker.bias[3]
            sensor_measure(walker, truth, ctx, rng)
            push!(steps, walker.bias[3] - before)
        end
        @test std(steps) ≈ 1e-4 * sqrt(1.0) rtol=0.1                                      # a random walk of the bias
        mag = Magnetometer(noise_t=1e-7, bias_t=SVector(0.0, 0.0, 5e-7))
        readings = [sensor_measure(mag, truth, ctx, rng) for _ in 1:4000]
        @test std(getindex.(readings, 1)) ≈ 1e-7 rtol=0.1
        @test mean(getindex.(readings, 3)) ≈ 3e-5 + 5e-7 atol=1e-8
        tracker = StarTracker(cross_boresight_noise_rad=1e-4, about_boresight_noise_rad=1e-3)
        angles = [2 * acos(clamp(abs(dot(sensor_measure(tracker, truth, ctx, rng), truth.q_ib)), 0.0, 1.0)) for _ in 1:2000]
        # A rotation vector with two components of 1e-4 and one of 1e-3: its RMS angle.
        @test sqrt(mean(angles .^ 2)) ≈ sqrt(2 * 1e-4^2 + 1e-3^2) rtol=0.1
    end

    @testset "noise-free sensors measure the truth in their frame" begin
        turn = SVector(0.0, 0.0, sin(π / 4), cos(π / 4))                                 # the panel turned 90 degrees about z
        sc = sensor_spacecraft(panel_q=turn)
        res = sensor_run(sc; mission_time=4.0)
        gyro, tracker, mag = (e.sensor for e in mounted_sensors(sc))
        @test sensor_measurement(gyro) ≈ SM.rot(turn) * ω0 atol=1e-12                      # torque-free spin of a cube
        k = row_at(res.table, sensor_measurement_time(tracker))
        @test k !== nothing
        q_row = SVector(res.table.sc1_q_1[k], res.table.sc1_q_2[k], res.table.sc1_q_3[k], res.table.sc1_q_4[k])
        @test sensor_measurement(tracker) ≈ q_row atol=1e-12
        @test 1e-5 < norm(sensor_measurement(mag)) < 1e-4                                   # tens of microtesla at 500 km
        @test sensor_measurement_time(gyro) == 4.0 && sensor_measurement_time(tracker) == 4.0
        @test "sc1_sensor_1" in names(res.table) && "sc1_sensor_10" in names(res.table)     # 3 + 4 + 3 components
        @test !("sc1_sensor_11" in names(res.table))
    end

    @testset "samples fall at the run's start and every period, also across checkpoint segments" begin
        clock = SampleClock(0.25; seed=4)
        sensor_run(sensor_spacecraft(sensors=(clock,)); mission_time=5.0)
        @test clock.times ≈ collect(0.0:0.25:5.0) atol=1e-12
        straight_last = sensor_measurement(clock)
        split = SampleClock(0.25; seed=4)
        sensor_run(sensor_spacecraft(sensors=(split,)); mission_time=5.0,
            settings=(checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=mktempdir()))
        @test split.times ≈ clock.times atol=1e-12                                          # each sample once
        @test sensor_measurement(split) == straight_last                                     # the same stream draws
        other = SampleClock(0.25; seed=5)
        sensor_run(sensor_spacecraft(sensors=(other,)); mission_time=5.0)
        @test sensor_measurement(other) != straight_last
        # A run resumed from the file samples on the same grid after the last
        # sample taken, continuing the same random stream.
        resume_dir = mktempdir()
        sensor_run(sensor_spacecraft(sensors=(SampleClock(0.25; seed=4),)); mission_time=2.0,
            settings=(checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=resume_dir))
        resumed = SampleClock(0.25; seed=4)
        sensor_run(sensor_spacecraft(sensors=(resumed,)); mission_time=5.0,
            settings=(checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=resume_dir, resume_from_checkpoint=true))
        @test resumed.times ≈ collect(2.25:0.25:5.0) atol=1e-12
        @test sensor_measurement(resumed) == straight_last
        # Another sensor configuration, or none, is refused.
        @test_throws ArgumentError sensor_run(sensor_spacecraft(sensors=(SampleClock(0.25; seed=5),)); mission_time=5.0,
            settings=(resume_from_checkpoint=true, checkpoint_directory=resume_dir))
        @test_throws ArgumentError sensor_run(sensor_spacecraft(sensors=()); mission_time=5.0,
            settings=(resume_from_checkpoint=true, checkpoint_directory=resume_dir))
    end

    @testset "a seed repeats a noisy gyro and another seed changes it" begin
        noisy(seed) = RateGyro(period_s=0.1, angle_random_walk=1e-3, rate_random_walk=1e-4, seed=seed)
        a, b, c = noisy(7), noisy(7), noisy(8)
        for g in (a, b, c)
            sensor_run(sensor_spacecraft(sensors=(g,)); mission_time=3.0)
        end
        @test sensor_measurement(a) == sensor_measurement(b) && a.bias == b.bias
        @test sensor_measurement(a) != sensor_measurement(c)
        @test norm(sensor_measurement(a) - ω0) < 0.1                                         # noise of a few 1e-3
    end

    @testset "navigation effectors read the latest measurement" begin
        log = ReadingLog(Vector{Float64}[], Float64[])
        sensor_run(sensor_spacecraft(sensors=(RateGyro(period_s=0.5),)); mission_time=5.0, navigation=(log,))
        @test !isempty(log.readings) && all(r -> norm(r - ω0) < 1e-12, log.readings)
        # Navigation ticks every second; each reads the sample taken at its own time,
        # since sensors sample before navigation effectors at the same time.
        @test length(log.times) >= 4 && log.times ≈ collect(1.0:1.0:length(log.times))
    end

    @testset "a one-component sensor may measure a number" begin
        gauge = MassGauge()
        res = sensor_run(sensor_spacecraft(sensors=(gauge,)); mission_time=3.0)
        @test sensor_measurement(gauge) isa Float64 && sensor_measurement(gauge) ≈ 50.0
        @test res.table.sc1_sensor_1[end] ≈ 50.0 && !("sc1_sensor_2" in names(res.table))
    end

    @testset "solver modes that cannot sample sensors refuse them" begin
        err = try
            sensor_run(sensor_spacecraft(sensors=(RateGyro(),)); mission_time=1.0, solver=:multirate)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("sensors", err.msg)
        reason = SpaceAGORA.SimulationEngine._gravity_backbone_reject_reason(
            sensor_args(sensor_spacecraft(sensors=(RateGyro(),)); orientation_sim=false))
        @test reason !== nothing && occursin("sensors", reason)
        @test !occursin("sensors", something(SpaceAGORA.SimulationEngine._gravity_backbone_reject_reason(
            sensor_args(sensor_spacecraft(sensors=()); orientation_sim=false)), ""))
    end

    @testset "GNSS, sun sensor and accelerometer records and noise" begin
        @test_throws ArgumentError GNSSReceiver(position_noise_m=-1.0)
        @test_throws ArgumentError SunSensor(field_of_view_rad=0.0)
        @test_throws ArgumentError SunSensor(min_lit_fraction=1.5)
        @test_throws ArgumentError Accelerometer(noise_mps2=Inf)
        @test sensor_dimension(GNSSReceiver()) == 6 && sensor_requirements(SunSensor()).solar && sensor_needs_specific_force(Accelerometer())
        @test !sensor_needs_specific_force(RateGyro()) && is_gravitational(SM.InverseSquaredJ2GravityModel()) && !is_gravitational(HoldLevels([0.0]))
        ctx = actuator_context(SM.Link(root=true, m=1.0))
        rng = Xoshiro(21)
        truth = SensorTruth(SVector(7e6, 0.0, 0.0), SVector(0.0, 7.5e3, 0.0), SVector(0.0, 0.0, 0.0, 1.0), ω0, 50.0, 0.0,
                            SVector(NaN, NaN, NaN), SVector(0.0, 0.0, 1.0), 1.0, SVector(0.01, 0.0, 0.0))
        fixes = [sensor_measure(GNSSReceiver(position_noise_m=5.0, velocity_noise_mps=0.1), truth, ctx, rng) for _ in 1:4000]
        @test std(getindex.(fixes, 2)) ≈ 5.0 rtol=0.1
        @test std(getindex.(fixes, 6)) ≈ 0.1 rtol=0.1
        @test sensor_measure(GNSSReceiver(), truth, ctx, rng) == SVector(7e6, 0.0, 0.0, 0.0, 7.5e3, 0.0)
        sun = SunSensor(noise_rad=1e-3)
        angles = [acos(clamp(dot(sensor_measure(sun, truth, ctx, rng), SVector(0.0, 0.0, 1.0)), -1.0, 1.0)) for _ in 1:4000]
        @test sqrt(mean(angles .^ 2)) ≈ sqrt(2) * 1e-3 rtol=0.1                          # two axes turn the direction
        side = SensorTruth(truth.pos_ii, truth.vel_ii, truth.q_ib, ω0, 50.0, 0.0, truth.magnetic_field_body, SVector(1.0, 0.0, 0.0), 1.0, truth.specific_force_body)
        @test all(isnan, sensor_measure(SunSensor(field_of_view_rad=π / 4), side, ctx, rng))    # outside the field of view
        shadow = SensorTruth(truth.pos_ii, truth.vel_ii, truth.q_ib, ω0, 50.0, 0.0, truth.magnetic_field_body, SVector(0.0, 0.0, 1.0), 0.2, truth.specific_force_body)
        @test all(isnan, sensor_measure(SunSensor(), shadow, ctx, rng))                  # eclipsed
        acc = Accelerometer(noise_mps2=1e-4, bias_mps2=SVector(0.0, 2e-4, 0.0))
        reads = [sensor_measure(acc, truth, ctx, rng) for _ in 1:4000]
        @test mean(getindex.(reads, 1)) ≈ 0.01 atol=2e-5
        @test mean(getindex.(reads, 2)) ≈ 2e-4 atol=2e-5
        @test std(getindex.(reads, 3)) ≈ 1e-4 rtol=0.1
        @test SensorTruth(truth.pos_ii, truth.vel_ii, truth.q_ib, ω0, 50.0, 0.0, truth.magnetic_field_body).sun_lit_fraction === NaN
    end

    @testset "noise-free GNSS and sun sensor measure the truth; the sun sensor sees the eclipse" begin
        gnss, sun = GNSSReceiver(period_s=1.0), SunSensor(period_s=1.0, field_of_view_rad=π)
        sc = sensor_spacecraft(sensors=(gnss, sun), direction=sun_hat())
        res = sensor_run(sc; mission_time=4.0)
        x = res.solution.u[end].sc[1]
        @test sensor_measurement(gnss) == SVector(x.pos..., x.vel...)
        e = SM.SimpleEphemeridesModel()
        s_ii = SVector{3, Float64}(SM.ephemerides_sun_ii(e, SM.Earth(), SM.ephemerides_epoch_seconds(e, EPOCH) + 4.0)) * SM.mean_solar_distance_m(SM.Earth())
        expected = SM.rot(SVector{4}(x.q)) * normalize(s_ii - SVector{3}(x.pos))
        @test sensor_measurement(sun) ≈ expected atol=1e-12
        night = SunSensor(period_s=1.0, field_of_view_rad=π)
        sensor_run(sensor_spacecraft(sensors=(night,), direction=-sun_hat()); mission_time=2.0)
        @test all(isnan, sensor_measurement(night)) && sensor_measurement_time(night) == 2.0   # in Earth's shadow
    end

    @testset "an accelerometer reads zero in free fall and the thrust under thrust" begin
        free = Accelerometer(period_s=1.0)
        sensor_run(sensor_spacecraft(sensors=(free,)); mission_time=5.0, effectors=(SM.InverseSquaredJ2GravityModel(),))
        @test norm(sensor_measurement(free)) < 1e-7                                      # gravity is not felt
        pushed = Accelerometer(period_s=1.0)
        sc = sensor_spacecraft(sensors=(pushed,))
        add_actuator!(sc, sc.root, Thruster(max_thrust=2.0, direction=MVector(1.0, 0.0, 0.0)))
        sensor_run(sc; mission_time=5.0, effectors=(SM.InverseSquaredGravityModel(),), control=(HoldLevels([1.0]),))
        @test sensor_measurement(pushed) ≈ SVector(1.0 / 50.0, 0.0, 0.0) atol=1e-7           # 1 N on 50 kg, along body x
        own = Accelerometer(period_s=1.0)                                                # gravity by a model of the user's
        sensor_run(sensor_spacecraft(sensors=(own,)); mission_time=5.0, effectors=(SM.InverseSquaredGravityModel(), UniformGravity()))
        @test norm(sensor_measurement(own)) < 1e-7
    end

    @testset "an accelerometer away from the centre of mass reads the rotation" begin
        # Torque-free spin about z: ω × (ω × r) at 1 m along x.
        spun = Accelerometer(period_s=1.0, position_m=SVector(1.0, 0.0, 0.0))
        sensor_run(sensor_spacecraft(sensors=(spun,), spin=SVector(0.0, 0.0, 0.1)); mission_time=5.0,
            effectors=(SM.InverseSquaredGravityModel(),))
        @test sensor_measurement(spun) ≈ SVector(-0.01, 0.0, 0.0) atol=1e-7
        # A wheel spins the body up about z at a constant rate: ω̇ × r joins the centripetal term.
        pushed = Accelerometer(period_s=1.0, position_m=SVector(1.0, 0.0, 0.0))
        sc = sensor_spacecraft(sensors=(pushed,), spin=SVector(0.0, 0.0, 0.0))
        add_actuator!(sc, sc.root, ReactionWheel(spin_axis=SVector(0.0, 0.0, 1.0)))
        x = sensor_run(sc; mission_time=5.0, effectors=(SM.InverseSquaredGravityModel(),),
            control=(HoldLevels([0.01]),)).solution.u[end].sc[1]
        α = x.ω[3] / 5.0
        @test sensor_measurement(pushed) ≈ SVector(-x.ω[3]^2, α, 0.0) rtol=1e-6 atol=1e-12
    end

    @testset "GNSS in the planet-fixed frame and above its coverage" begin
        fixed = GNSSReceiver(period_s=1.0, frame=:planet_fixed)
        above = GNSSReceiver(period_s=1.0, max_altitude_m=400e3)                       # the orbit is at 500 km
        below = GNSSReceiver(period_s=1.0, max_altitude_m=600e3)
        @test_throws ArgumentError GNSSReceiver(frame=:ecef)
        @test sensor_requirements(fixed).planet_frame && !sensor_requirements(GNSSReceiver()).planet_frame
        x = sensor_run(sensor_spacecraft(sensors=(fixed, above, below)); mission_time=3.0).solution.u[end].sc[1]
        m = sensor_measurement(fixed)
        @test norm(m[1:3]) ≈ norm(x.pos) rtol=1e-12
        @test norm(m[1:3] - collect(x.pos)) > 1.0e3                                     # the planet's rotation
        @test all(isnan, sensor_measurement(above))
        @test sensor_measurement(below) == SVector(x.pos..., x.vel...)
    end

    @testset "a sun sensor with albedo leans toward the lit planet" begin
        ctx = actuator_context(SM.Link(root=true, m=1.0))
        rng = Xoshiro(3)
        R = SM.Earth().Rp_e
        r = R + 500e3
        nan3 = SVector(NaN, NaN, NaN)
        truth(sun) = SensorTruth(SVector(r, 0.0, 0.0), SVector(0.0, 7.6e3, 0.0), SVector(0.0, 0.0, 0.0, 1.0), ω0, 50.0, 0.0,
                                 nan3, sun, 1.0, nan3, nan3, nan3, nan3, NaN, R)
        bright = SunSensor(field_of_view_rad=π, albedo=0.3)
        @test_throws ArgumentError SunSensor(albedo=1.5)
        η = 0.3 * (2 / 3) * (R / r)^2 / π                                              # α = 90°: Φ = 1/π
        @test sensor_measure(bright, truth(SVector(0.0, 1.0, 0.0)), ctx, rng) ≈ normalize(SVector(-η, 1.0, 0.0))
        @test sensor_measure(SunSensor(field_of_view_rad=π), truth(SVector(0.0, 1.0, 0.0)), ctx, rng) == SVector(0.0, 1.0, 0.0)
        @test sensor_measure(bright, truth(SVector(-1.0, 0.0, 0.0)), ctx, rng) == SVector(-1.0, 0.0, 0.0)   # over the night side
        narrow = SunSensor(field_of_view_rad=0.2, albedo=0.3)                          # the planet's disk out of view
        @test sensor_measure(narrow, truth(SVector(0.0, 0.0, 1.0)), ctx, rng) == SVector(0.0, 0.0, 1.0)
        @test truth(SVector(0.0, 0.0, 1.0)).planet_radius_m == R
        @test isnan(SensorTruth(SVector(r, 0.0, 0.0), SVector(0.0, 7.6e3, 0.0), SVector(0.0, 0.0, 0.0, 1.0), ω0, 50.0, 0.0,
            nan3, SVector(0.0, 0.0, 1.0), 1.0, nan3).planet_radius_m)                   # the ten-argument form
    end

    @testset "sensors are opt-in and do not change the dynamics" begin
        # Sampling adds solver stops, so the comparison is to tolerance, not bit for bit.
        with = sensor_run(sensor_spacecraft(); mission_time=5.0).solution.u[end].sc[1]
        without = sensor_run(sensor_spacecraft(sensors=()); mission_time=5.0)
        @test !any(n -> occursin("sensor", n), names(without.table))
        @test collect(with.ω) ≈ collect(without.solution.u[end].sc[1].ω) rtol=1e-8
        @test collect(with.pos) ≈ collect(without.solution.u[end].sc[1].pos) rtol=1e-12
    end
end

end # module SensorTests
