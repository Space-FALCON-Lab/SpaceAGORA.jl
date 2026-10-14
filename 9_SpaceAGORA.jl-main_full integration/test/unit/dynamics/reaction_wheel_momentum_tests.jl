module ReactionWheelMomentumTests
using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA

const RWM = SpaceAGORA.SimulationModel

@testset "ReactionWheelMomentum" begin

    # =======================================================================
    # The spline the wheel speeds are read through
    # =======================================================================
    @testset "wheel-speed spline" begin
        @testset "a straight line is reproduced exactly, value and derivative" begin
            # A natural cubic spline sets the second derivative to zero at both
            # ends, so a straight line is the one polynomial it reproduces
            # everywhere with no end effect at all.
            t = collect(0.0:0.25:10.0)
            s = RWM.wheel_speed_spline(t, 3.0 .- 1.5 .* t)
            for x in (0.0, 0.1, 3.7, 6.25, 9.9, 10.0)
                @test RWM.wheel_spline_value(s, x) ≈ 3.0 - 1.5 * x atol = 1e-12
                @test RWM.wheel_spline_derivative(s, x) ≈ -1.5 atol = 1e-12
            end
        end

        @testset "it interpolates its knots and is smooth between them" begin
            t = collect(0.0:0.5:20.0)
            y = sin.(0.3 .* t)
            s = RWM.wheel_speed_spline(t, y)
            for k in eachindex(t)
                @test RWM.wheel_spline_value(s, t[k]) ≈ y[k] atol = 1e-12
            end
            # Away from the ends a half-second sampling of a 0.3 rad/s sine is
            # resolved to well under a percent, and the analytic derivative
            # agrees with a central difference of the spline itself.
            for x in (4.3, 9.1, 15.7)
                @test RWM.wheel_spline_value(s, x) ≈ sin(0.3 * x) atol = 1e-4
                h = 1e-5
                fd = (RWM.wheel_spline_value(s, x + h) - RWM.wheel_spline_value(s, x - h)) / (2h)
                @test RWM.wheel_spline_derivative(s, x) ≈ fd atol = 1e-7
            end
        end

        @testset "queries outside the table clamp instead of extrapolating" begin
            t = collect(0.0:0.5:5.0)
            s = RWM.wheel_speed_spline(t, 2.0 .* t)
            @test RWM.wheel_spline_value(s, -3.0) ≈ 0.0 atol = 1e-12
            @test RWM.wheel_spline_value(s, 12.0) ≈ 10.0 atol = 1e-12
            @test RWM.wheel_spline_derivative(s, -3.0) == 0.0
            @test RWM.wheel_spline_derivative(s, 12.0) == 0.0
        end

        @testset "non-uniform knots are supported and bad input is rejected" begin
            t = [0.0, 0.3, 0.9, 1.0, 2.5, 4.0]
            s = RWM.wheel_speed_spline(t, 1.0 .+ 0.5 .* t)
            @test RWM.wheel_spline_value(s, 1.7) ≈ 1.0 + 0.5 * 1.7 atol = 1e-12
            @test_throws ArgumentError RWM.wheel_speed_spline([0.0, 1.0], [0.0, 1.0])
            @test_throws ArgumentError RWM.wheel_speed_spline([0.0, 1.0, 1.0, 2.0], zeros(4))
            @test_throws ArgumentError RWM.wheel_speed_spline([0.0, 1.0, 2.0], zeros(2))
        end
    end

    # =======================================================================
    # The torque law, against cases whose answer is known in closed form
    # =======================================================================
    @testset "wheel reaction torque" begin
        axes = Matrix{Float64}(I, 3, 3)
        inertia_wheel = 2.0e-5
        t = collect(0.0:0.25:100.0)

        @testset "constant wheel speed on a non-rotating body gives exactly zero" begin
            speeds = repeat([400.0 -250.0 120.0], length(t), 1)
            m = RWM.ReactionWheelMomentumModel(t, speeds, axes, inertia_wheel)
            for x in (0.0, 13.0, 57.5, 100.0)
                h = RWM.wheel_momentum_body(m, x)
                @test h ≈ SVector{3, Float64}(400.0, -250.0, 120.0) .* inertia_wheel atol = 1e-15
                @test RWM.wheel_momentum_rate_body(m, x) ≈ SVector{3, Float64}(0, 0, 0) atol = 1e-15
                τ = RWM.wheel_reaction_torque(h, RWM.wheel_momentum_rate_body(m, x), SVector{3, Float64}(0, 0, 0))
                @test τ ≈ SVector{3, Float64}(0, 0, 0) atol = 1e-15
            end
        end

        @testset "constant wheel speed on a rotating body gives the gyroscopic term alone" begin
            speeds = repeat([400.0 0.0 0.0], length(t), 1)
            m = RWM.ReactionWheelMomentumModel(t, speeds, axes, inertia_wheel)
            ω = SVector{3, Float64}(0.0, 1.0e-3, -4.0e-4)
            h = RWM.wheel_momentum_body(m, 40.0)
            τ = RWM.wheel_reaction_torque(h, RWM.wheel_momentum_rate_body(m, 40.0), ω)
            @test τ ≈ -cross(ω, h) atol = 1e-18
        end

        @testset "a linear speed ramp gives a constant torque along the wheel axis" begin
            rate = 3.0                       # rad/s per second
            speeds = hcat(rate .* t, zeros(length(t)), zeros(length(t)))
            m = RWM.ReactionWheelMomentumModel(t, speeds, axes, inertia_wheel)
            expected = SVector{3, Float64}(-rate * inertia_wheel, 0.0, 0.0)
            for x in (5.0, 30.0, 72.0, 95.0)
                τ = RWM.wheel_reaction_torque(RWM.wheel_momentum_body(m, x),
                    RWM.wheel_momentum_rate_body(m, x), SVector{3, Float64}(0, 0, 0))
                @test τ ≈ expected atol = 1e-15
            end
        end

        @testset "a skewed, non-unit axis set scales the momentum by its own columns" begin
            skew = [0.8 0.0 0.3; 0.0 1.2 -0.4; 0.6 0.5 0.9]
            speeds = repeat([100.0 -50.0 25.0], length(t), 1)
            m = RWM.ReactionWheelMomentumModel(t, speeds, skew, inertia_wheel)
            @test RWM.wheel_momentum_body(m, 20.0) ≈
                SVector{3, Float64}(skew * (inertia_wheel .* [100.0, -50.0, 25.0])) atol = 1e-15
        end

        @testset "the time offset shifts the table's clock" begin
            speeds = hcat(2.0 .* t, zeros(length(t)), zeros(length(t)))
            plain = RWM.ReactionWheelMomentumModel(t, speeds, axes, inertia_wheel)
            shifted = RWM.ReactionWheelMomentumModel(t, speeds, axes, inertia_wheel; time_offset_s=30.0)
            @test RWM.wheel_momentum_body(shifted, 10.0) ≈ RWM.wheel_momentum_body(plain, 40.0) atol = 1e-15
        end

        @testset "malformed construction is rejected" begin
            @test_throws ArgumentError RWM.ReactionWheelMomentumModel(t, zeros(length(t), 3), zeros(3, 2), inertia_wheel)
            @test_throws ArgumentError RWM.ReactionWheelMomentumModel(t, zeros(length(t) - 1, 3), axes, inertia_wheel)
            @test_throws ArgumentError RWM.ReactionWheelMomentumModel(t, zeros(length(t), 3), axes, -1.0)
        end
    end

    # =======================================================================
    # The invariant the whole scenario rests on
    # =======================================================================
    @testset "a torque-free body conserves total angular momentum" begin
        # Integrate the toolkit's own rigid-body right-hand side with the wheel
        # reaction torque as its only forcing, and check that the TOTAL
        # momentum, body momentum plus wheel momentum carried into the inertial
        # frame, does not move. This is the statement the CYGNSS slew scenario
        # is built on: the wheels exchange momentum with the body and neither
        # creates it.
        inertia = SMatrix{3, 3, Float64}(1.4, -0.0171, 0.00808, -0.0171, 0.819, -0.00535, 0.00808, -0.00535, 1.95)
        inertia_wheel = 2.8648e-5
        axes = [0.9 0.1 -0.2; -0.1 1.0 0.3; 0.2 -0.3 0.95]
        t = collect(0.0:0.25:200.0)
        # Wheel speeds that actually move: a ramp, a sine and a step-like arctan.
        speeds = hcat(300.0 .- 4.0 .* t, 250.0 .* sin.(0.05 .* t), 100.0 .* atan.(0.2 .* (t .- 100.0)))
        model = RWM.ReactionWheelMomentumModel(t, speeds, axes, inertia_wheel)

        q = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
        ω = SVector{3, Float64}(1.0e-3, -5.0e-4, 2.0e-4)
        total0 = inertia * ω + RWM.wheel_momentum_body(model, 0.0)   # inertial, since q is identity

        function derivative(q, ω, time)
            h = RWM.wheel_momentum_body(model, time)
            hdot = RWM.wheel_momentum_rate_body(model, time)
            τ = RWM.wheel_reaction_torque(h, hdot, ω)
            return (RWM.DynamicsRotational.quaternion_derivative(ω, q),
                RWM.DynamicsRotational.angular_acceleration(ω, inertia, τ))
        end

        dt = 0.01
        drift = 0.0
        time = 0.0
        while time < 200.0 - dt / 2
            k1q, k1w = derivative(q, ω, time)
            k2q, k2w = derivative(q + 0.5dt * k1q, ω + 0.5dt * k1w, time + 0.5dt)
            k3q, k3w = derivative(q + 0.5dt * k2q, ω + 0.5dt * k2w, time + 0.5dt)
            k4q, k4w = derivative(q + dt * k3q, ω + dt * k3w, time + dt)
            q = q + (dt / 6) * (k1q + 2k2q + 2k3q + k4q)
            ω = ω + (dt / 6) * (k1w + 2k2w + 2k3w + k4w)
            q = q / norm(q)
            time += dt
            total = RWM.rot(q)' * (inertia * ω + RWM.wheel_momentum_body(model, time))
            drift = max(drift, norm(total - total0))
        end
        # The wheels swing the body momentum by several times its own initial
        # size over these 200 s, so a relative drift of 1e-6 is a real test of
        # the torque law rather than of the integrator's step size.
        @test drift / norm(total0) < 1.0e-6
    end
end

# Additional tests for the bounded wheel extraction. API/fixture baseline:
# main 62199228; original wheel model cc891c76. Run after applying the candidate.
# No native GRAM, SPICE assets, output files or performance measurements required.
module ReactionWheelMomentumIntegrationTests
using Test
using LinearAlgebra
using StaticArrays
using ComponentArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const ZERO3 = SVector{3,Float64}(0.0, 0.0, 0.0)
const BODY_I = SMatrix{3,3,Float64}(2.0,0.1,-0.05, 0.1,3.0,0.07, -0.05,0.07,4.0)
const WHEEL_AXES = [1.0 0.2 -0.1 0.35; -0.1 0.9 0.3 -0.4; 0.15 -0.2 1.1 0.9]
const WHEEL_I = 0.01
const SPEED0 = SVector{4,Float64}(20.0, 10.0, -15.0, 12.0)
const SPEED_RATE = SVector{4,Float64}(5.0, -2.0, 1.0, -3.0)
const TABLE_OFFSET = 0.5
const Q0 = normalize(SVector{4,Float64}(0.1, -0.2, 0.15, 1.0))
const OMEGA0 = SVector{3,Float64}(0.03, -0.02, 0.01)
const MISSION_S = 12.0

expected_h(t) = SVector{3,Float64}(WHEEL_AXES * (WHEEL_I .* (SPEED0 + SPEED_RATE .* (t + TABLE_OFFSET))))
expected_hdot() = SVector{3,Float64}(WHEEL_AXES * (WHEEL_I .* SPEED_RATE))

function wheel_fixture(; include_wheels=true)
    planet = Earth()
    spacecraft = SpacecraftModel[]
    for (index, id) in enumerate((41, 97))
        root = Link(root=true, m=500.0, ref_area=12.0, inertia=BODY_I)
        ic = InitialCondition(
            ra=planet.Rp_e + 550_000.0 + 100.0 * index,
            rp=planet.Rp_e + 550_000.0 + 100.0 * index,
            i=53.0, ω=0.0, Ω=10.0, ν=180.0 * (index - 1), q=Q0, ang_vel=OMEGA0)
        push!(spacecraft, SpacecraftModel(Joint[], [root], root, true,
            500.0, 0.0, BODY_I, 0, 0, ic, id))
    end
    knots = collect(-2.0:2.0:16.0)
    speeds = [SPEED0[k] + SPEED_RATE[k] * t for t in knots, k in 1:4]
    wheel = SM.ReactionWheelMomentumModel(knots, speeds, WHEEL_AXES, WHEEL_I;
        spacecraft_index=2, num_spacecraft=2, time_offset_s=TABLE_OFFSET)
    effectors = include_wheels ? (InverseSquaredGravityModel(), wheel) :
                                (InverseSquaredGravityModel(),)
    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            keplerian=true, number_of_orbits=1, mission_time=MISSION_S,
            orientation_sim=true, num_steps_to_save=32, data_rate=0.5),
        environment_model=EnvironmentModel(planet=planet, EI=300.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false, ephemerides_model=SimpleEphemeridesModel()),
        dynamics_model=DynamicsModel(spacecraft, effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-10, abstol_orbit=1e-10,
            reltol_quaternion=1e-11, abstol_quaternion=1e-12,
            reltol_angular_rate=1e-11, abstol_angular_rate=1e-12,
            dt_max_orbit=0.05),
        solver_config=SolverConfig(solver_mode=:tsit5))
    return args, wheel
end

function saved_wheel_momentum(u, t, integrator)
    effectors = integrator.p.args.dynamics_model.dynamic_effectors
    wheel_index = findfirst(e -> e isa SM.ReactionWheelMomentumModel, effectors)
    wheel_index === nothing && return fill(ZERO3, length(u.sc))
    model = effectors[wheel_index]
    # Query at the save time; diagnostics contain the last RHS stage instead.
    return [i == model.spacecraft_index ? SM.wheel_momentum_body(model, Float64(t)) : ZERO3
            for i in eachindex(u.sc)]
end

function run_recorded(args)
    fields = SaveField[
        SaveField(:q, (u,t,integrator) -> [SVector{4,Float64}(sc.q) for sc in u.sc]; per_satellite=true),
        SaveField(:omega, (u,t,integrator) -> [SVector{3,Float64}(sc.ω) for sc in u.sc]; per_satellite=true),
        SaveField(:wheel_h, saved_wheel_momentum; per_satellite=true),
    ]
    recorder = TrajectoryRecorder(args; save_fields=fields)
    metadata = withenv(
        "SPACEAGORA_RHS_EXECUTION_MODE"=>"serial",
        "SPACEAGORA_RHS_CALIBRATE"=>"off", "SPACEAGORA_RHS_IDENTIFY"=>"0",
        "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS"=>"0",
        "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST"=>"0",
    ) do
        # Intentionally omit isolate_state: this tests the public default.
        run_simulation(args; return_solver_metadata=true, visualization=false,
            extra_callbacks=(get_trajectory_recorder_callback(recorder),))
    end
    @test metadata.retcode == "Success"
    @test recorder.count >= 20
    @test first(trajectory_times(recorder)) ≈ 0.0 atol=1e-12
    @test last(trajectory_times(recorder)) ≈ MISSION_S atol=1e-12
    return collect(trajectory_times(recorder)), trajectory_save_data(recorder)
end

source_snapshot(args, wheel) = deepcopy((
    momentum=wheel.state.momentum_body, torque=wheel.state.torque_body,
    initial=collect(ComponentArrays.getdata(SE.build_initial_conditions(args))),
    roots=[(copy(sc.root.r), copy(sc.root.q), copy(sc.root.ω),
            copy(sc.root.net_force), copy(sc.root.net_torque)) for sc in args.dynamics_model.spacecraft]))

@testset "wheel wrench purity and explicit spacecraft-index selection" begin
    args, model = wheel_fixture()
    state = model.state
    fill!(state.momentum_body, SVector(1.0, 2.0, 3.0))
    fill!(state.torque_body, SVector(4.0, 5.0, 6.0))
    before = deepcopy((state.momentum_body, state.torque_body))
    sc = args.dynamics_model.spacecraft[2]
    sample = StateSample(SVector(7.0e6,0.0,0.0), SVector(0.0,7.5e3,0.0), 500.0;
        q_ib=Q0, ω_body=OMEGA0, spacecraft=sc)
    env = EnvironmentSample(args.environment_model.planet)
    time = 1.25
    expected = -expected_hdot() - cross(OMEGA0, expected_h(time))
    force, torque = SM.wrench(model, sample, env, time)
    @test force == ZERO3
    @test torque ≈ expected rtol=1e-13 atol=1e-14
    @test (state.momentum_body, state.torque_body) == before
    # The spacecraft object's ID is 97, while its integration index is 2.
    @test sc.id == 97
    @test SM.wrench_caching!(model, sample, env, time, nothing, 1) == (ZERO3, ZERO3)
    @test (state.momentum_body, state.torque_body) == before
    @test SM.wrench_caching!(model, sample, env, time, nothing, 2)[2] ≈ expected
    @test state.momentum_body[1] == before[1][1]
    @test state.torque_body[1] == before[2][1]
    @test state.momentum_body[2] ≈ expected_h(time)
    @test state.torque_body[2] ≈ expected
end

@testset "real propagation selects only its craft and conserves inertial momentum" begin
    args, model = wheel_fixture()
    @test [sc.id for sc in args.dynamics_model.spacecraft] == [41, 97]
    @test length(model.speeds) == 4
    @test isposdef(BODY_I) && !isdiag(BODY_I)
    @test all(sc.root.rw_assembly.n_wheels == 0 for sc in args.dynamics_model.spacecraft)
    before = source_snapshot(args, model)
    times, snapshots = run_recorded(args)
    @test source_snapshot(args, model) == before
    baseline_args, _ = wheel_fixture(include_wheels=false)
    baseline_times, baseline = run_recorded(baseline_args)
    @test times == baseline_times
    @test length(snapshots) == length(baseline)
    total0 = [SM.rot(Q0)' * (BODY_I * OMEGA0 + (i == 2 ? expected_h(0.0) : ZERO3)) for i in 1:2]
    relative_drift = zeros(2)
    for (t, snapshot, without) in zip(times, snapshots, baseline)
        @test snapshot[:wheel_h][1] == ZERO3
        @test snapshot[:wheel_h][2] ≈ expected_h(t) rtol=1e-13 atol=1e-14
        @test snapshot[:q][1] ≈ without[:q][1] rtol=0.0 atol=1e-8
        @test snapshot[:omega][1] ≈ without[:omega][1] rtol=0.0 atol=1e-8
        for i in 1:2
            q, omega, wheel_h = snapshot[:q][i], snapshot[:omega][i], snapshot[:wheel_h][i]
            @test abs(norm(q) - 1.0) < 1e-8
            total = SM.rot(q)' * (BODY_I * omega + wheel_h)
            relative_drift[i] = max(relative_drift[i], norm(total-total0[i]) / norm(total0[i]))
        end
    end
    @test maximum(relative_drift) < 2e-7
    @test norm(snapshots[end][:omega][2] - baseline[end][:omega][2]) > 1e-3
    @test norm(snapshots[end][:wheel_h][2] - snapshots[1][:wheel_h][2]) > 0.1
end

@testset "serial dispatch and flat assembly agree on the wheel wrench" begin
    args, _ = wheel_fixture()
    u = SE.build_initial_conditions(args)
    serial_p = ODEParams(n_sats=2, args=deepcopy(args))
    serial_du = zero(u)
    serial_p.shared_buffers.rhs_env_config[] = withenv("SPACEAGORA_RHS_EXECUTION_MODE"=>"serial") do
        SE._snapshot_rhs_plan_env_config()
    end
    SE.spacecraft_dynamics!(serial_du, u, serial_p, 1.25)
    # Test the existing flat assembly directly. Production forced-flat routing
    # rejects unlisted effectors, including the original wheel model. This test
    # does not imply a change to that allowlist or promise a threaded speedup.
    for width in unique([1, min(2, Threads.nthreads())])
        p = ODEParams(n_sats=2, args=deepcopy(args))
        du = zero(u)
        plan = (mode=:flat_constellation_effector_queue, allotment=width, scheduler=:static,
            dominant_axis=:effector, policy_applied=false,
            effector_decision=(use_threads=false, allotment=1, mode=:off, policy_applied=false))
        SE._spacecraft_dynamics_flat_constellation_effector_queue!(du, u, p, 1.25, plan; rhs_kind=:full)
        @test collect(ComponentArrays.getdata(du)) ≈ collect(ComponentArrays.getdata(serial_du)) rtol=1e-12 atol=1e-12
        wheel = p.args.dynamics_model.dynamic_effectors[2]
        @test wheel.state.momentum_body[1] == ZERO3
        @test wheel.state.torque_body[1] == ZERO3
        @test wheel.state.momentum_body[2] ≈ expected_h(1.25)
        @test wheel.state.torque_body[2] ≈ -expected_hdot() - cross(OMEGA0, expected_h(1.25))
    end
end

@testset "wheel construction rejects nonfinite inputs and both short state vectors" begin
    t, speeds, axes = [0.0,1.0,2.0], zeros(3,3), Matrix{Float64}(I,3,3)
    for bad in (NaN, Inf, -Inf)
        bad_axes = copy(axes)
        bad_axes[2,2] = bad
        @test_throws ArgumentError SM.ReactionWheelMomentumModel(t, speeds, bad_axes, 0.25)
        @test_throws ArgumentError SM.ReactionWheelMomentumModel(t, speeds, axes, 0.25; time_offset_s=bad)
    end
    for (nh, nt) in ((1,2), (2,1))
        state = SM.ReactionWheelMomentumState(fill(ZERO3,nh), fill(ZERO3,nt))
        @test_throws ArgumentError SM.ReactionWheelMomentumModel(t, speeds, axes, 0.25;
            spacecraft_index=2, state=state)
    end
end

@testset "raw and fitted spline constructors preserve their public contract" begin
    t, y = [0.0,1.0,2.0], [3.0,1.5,0.0]
    for spline in (SM.WheelSpeedSpline(t,y,zeros(3)), SM.WheelSpeedSpline(t,y))
        @test SM.wheel_spline_value(spline,0.5) ≈ 2.25
        @test SM.wheel_spline_derivative(spline,0.5) ≈ -1.5
    end
    for (tx, values, d2) in (
        (t,y[1:2],zeros(3)), (t,y,zeros(2)), (t[1:2],y[1:2],zeros(2)),
        ([0.0,1.0,1.0],y,zeros(3)), ([0.0,2.0,1.0],y,zeros(3)),
    )
        @test_throws ArgumentError SM.WheelSpeedSpline(tx,values,d2)
    end
    for bad in (NaN,Inf,-Inf)
        @test_throws ArgumentError SM.WheelSpeedSpline([0.0,bad,2.0],y,zeros(3))
        @test_throws ArgumentError SM.WheelSpeedSpline(t,[3.0,bad,0.0],zeros(3))
        @test_throws ArgumentError SM.WheelSpeedSpline(t,y,[0.0,bad,0.0])
    end
end

@testset "wheel-speed helper returns all K wheels without padding" begin
    times = [0.0,1.0,2.0]
    for k in (1,2,3,4)
        speeds = [(10.0+t)*j for t in times, j in 1:k]
        axes = zeros(3,k)
        axes[1,:] .= 1.0
        model = SM.ReactionWheelMomentumModel(times,speeds,axes,0.25; time_offset_s=0.25)
        values = SM.wheel_speeds_rad_s(model,0.5)
        @test values isa Vector{Float64}
        @test length(values) == k
        @test values ≈ [10.75*j for j in 1:k]
    end
end

end # module ReactionWheelMomentumIntegrationTests

end # module ReactionWheelMomentumTests
