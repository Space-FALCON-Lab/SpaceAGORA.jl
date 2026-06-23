# Implement and run ADS test suite (test/gnc/ads_tests.jl)

Create test/gnc/ads_tests.jl and run it to verify DRM 1 implementation.

Requirements: ADS-VNV-001, ADS-VNV-002, ADS-VNV-003.

IMPORTANT: After writing the test file, run it using:
  julia --project=C:/Users/Robbie/falcon/SpaceAGORA.jl test/gnc/ads_tests.jl
and report the test results.

using Test, LinearAlgebra, Statistics, Random
# Include ADS files
include("../../src/gnc/navigation/ads_sensor_models.jl")
include("../../src/gnc/navigation/ads_mekf.jl")
include("../../src/gnc/navigation/ads_mode_manager.jl")
include("../../src/gnc/navigation/ads_navigation_effector.jl")

@testset "ADS — DRM 1 Verification" begin

  @testset "ADS-VNV-001: Sensor Noise Statistics" begin
    # Star tracker: measure N times with fixed truth, verify noise std matches 5 arcsec
    # Gyro: verify ARW power spectral density
    # Magnetometer: verify 50 nT std across 1000 samples
    # Sun sensor: verify 1 deg std
    # Use chi-squared goodness-of-fit at p=0.05
    N = 1000
    # ... test implementation
  end

  @testset "ADS-VNV-002: MEKF Convergence" begin
    # Initialize with 5 deg attitude error
    q_true = [1.0, 0.0, 0.0, 0.0]  # identity quaternion
    # Create initial estimate with 5 deg rotation error around z-axis
    # Run 60 seconds of predict+update cycles at 1 Hz
    # Final attitude error must be < 0.1 deg per ADS-SYS-005
    @test final_error_deg < 0.1
  end

  @testset "ADS-VNV-003: Filter Consistency (NEES)" begin
    # 100 Monte Carlo trials
    # Each trial: random 5 deg initial error, run 100 update cycles
    # Compute NEES at each step
    # Verify: mean NEES within 95% CI bounds of chi2(3) distribution
    # CI bounds: [0.352, 7.815] for chi2(3) at 95%
    N_trials = 100
    # ... MC implementation
    fraction_in_bounds = count(nees_within_95pct) / length(all_nees)
    @test fraction_in_bounds >= 0.90  # 90% of steps within 95% CI
  end

  @testset "ADS mode transitions" begin
    mgr = ads_mode_init()
    @test mgr.current_mode == INIT
    # Test INIT -> COARSE_ACQ
    all_healthy = Dict(:star_tracker=>true,:gyro=>true,:magnetometer=>true,:sun_sensor=>true)
    ads_mode_update!(mgr, all_healthy, false, 0.0, 10.0)
    @test mgr.current_mode == COARSE_ACQ
    # Test COARSE_ACQ -> FINE_TRACK (error < 2 deg)
    ads_mode_update!(mgr, all_healthy, false, 0.0, 1.0)
    @test mgr.current_mode == FINE_TRACK
    # Test FINE_TRACK -> ECLIPSE
    ads_mode_update!(mgr, all_healthy, true, 0.0, 0.05)
    @test mgr.current_mode == ECLIPSE
    # Test -> SAFE on gyro failure
    no_gyro = Dict(:star_tracker=>true,:gyro=>false,:magnetometer=>true,:sun_sensor=>true)
    ads_mode_update!(mgr, no_gyro, false, 0.0, 0.05)
    @test mgr.current_mode == SAFE
  end

end

After running tests, report results via complete_run with test_results summary.
If tests fail, diagnose and fix the implementation files, re-run until passing.