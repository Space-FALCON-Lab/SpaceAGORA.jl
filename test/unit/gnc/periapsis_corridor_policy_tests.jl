using Test

include(joinpath(@__DIR__, "../../../examples/periapsis_corridor_controller.jl"))
using .PeriapsisCorridorController

# Altitudes here are GEODETIC periapsis heights; the policy only compares numbers, the
# conversion lives in the predictor (see periapsis_corridor_runtime_tests.jl).
@testset "Physical periapsis corridor policy" begin
    @testset "near-passive coasts safely and corrects either drift direction" begin
        safe = select_corridor_maneuver(dv -> 97_500.0 + 100.0 * dv;
            max_dv_mps=20.0)
        @test safe.feasible
        @test safe.delta_v_mps == 0.0
        @test safe.predicted_periapsis_altitude_m == 97_500.0
        @test safe.reason == :inside_safety_band

        for (drift_m, expected_sign) in ((-3_000.0, 1.0), (3_000.0, -1.0))
            result = select_corridor_maneuver(dv -> 97_500.0 + drift_m + 100.0 * dv;
                target_m=97_500.0, max_dv_mps=40.0)
            @test result.feasible
            @test result.delta_v_mps * expected_sign > 0.0
            @test isapprox(result.predicted_periapsis_altitude_m, 97_500.0; atol=0.02)
            @test result.reason == :target_reached
        end
    end

    @testset "frequent policy honors its target deadband" begin
        coast = select_corridor_maneuver(dv -> 97_500.0 + 100.0 * dv;
            policy=:frequent, target_m=97_500.0, deadband_m=100.0,
            max_dv_mps=20.0)
        @test coast.delta_v_mps == 0.0
        @test coast.reason == :inside_deadband

        unsafe = select_corridor_maneuver(dv -> 94_000.0 + 100.0 * dv;
            policy=:frequent, target_m=97_500.0, deadband_m=5_000.0,
            max_dv_mps=40.0)
        @test unsafe.feasible
        @test unsafe.delta_v_mps > 0.0
        @test unsafe.reason == :target_reached

        correction = select_corridor_maneuver(dv -> 97_500.0 + 500.0 + 100.0 * dv;
            policy=:frequent, target_m=97_500.0, deadband_m=100.0,
            max_dv_mps=20.0)
        @test correction.feasible
        @test correction.reason == :target_reached
        @test isapprox(correction.predicted_periapsis_altitude_m, 97_500.0; atol=0.02)
    end

    @testset "actuator limits are reported explicitly" begin
        limited = select_corridor_maneuver(dv -> 94_000.0 + 100.0 * dv;
            target_m=97_500.0, max_dv_mps=10.0)
        @test !limited.feasible
        @test limited.reason == :actuator_limit_unsafe
        @test abs(limited.delta_v_mps) <= 10.0
        @test limited.predicted_periapsis_altitude_m <= 95_000.0
    end

    @testset "monotonic nonlinear predictor is solved without assuming linearity" begin
        predictor = dv -> 97_500.0 + 2_000.0 * tanh(dv)
        result = select_corridor_maneuver(predictor; policy=:frequent,
            target_m=98_500.0, deadband_m=0.0, max_dv_mps=1.0)
        @test result.feasible
        @test result.reason == :target_reached
        @test isapprox(result.delta_v_mps, atanh(0.5); atol=1e-5)
        @test isapprox(result.predicted_periapsis_altitude_m, 98_500.0; atol=0.02)
    end

    @testset "only the burn side toward the target is evaluated" begin
        positive = select_corridor_maneuver(dv -> begin
            dv < -0.5 && error("retrograde reset would enter the ground")
            97_500.0 + 1_000.0 * dv
        end; policy=:frequent, target_m=97_600.0, deadband_m=0.0,
            max_dv_mps=20.0)
        @test positive.feasible
        @test positive.delta_v_mps > 0.0
        @test isapprox(positive.predicted_periapsis_altitude_m, 97_600.0; atol=0.02)

        negative = select_corridor_maneuver(dv -> begin
            dv > 0.5 && error("prograde reset would enter the ground")
            97_600.0 + 1_000.0 * dv
        end; policy=:frequent, target_m=97_500.0, deadband_m=0.0,
            max_dv_mps=20.0)
        @test negative.feasible
        @test negative.delta_v_mps < 0.0
        @test isapprox(negative.predicted_periapsis_altitude_m, 97_500.0; atol=0.02)
    end

    @testset "invalid predictor outputs fail closed" begin
        nan_result = select_corridor_maneuver(_ -> NaN)
        @test !nan_result.feasible
        @test nan_result.reason == :predictor_error

        nonmonotonic = select_corridor_maneuver(dv -> 94_000.0 - 100.0 * abs(dv);
            max_dv_mps=10.0)
        @test !nonmonotonic.feasible
        @test nonmonotonic.reason == :nonmonotonic_predictor
    end

    @testset "configuration is validated" begin
        @test_throws ArgumentError select_corridor_maneuver(_ -> 97_500.0;
            lower_m=100_000.0, upper_m=95_000.0)
        @test_throws ArgumentError select_corridor_maneuver(_ -> 97_500.0;
            margin_m=3_000.0)
        @test_throws ArgumentError select_corridor_maneuver(_ -> 97_500.0;
            policy=:unknown)
    end
end
