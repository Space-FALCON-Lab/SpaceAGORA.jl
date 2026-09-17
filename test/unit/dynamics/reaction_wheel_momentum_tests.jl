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
