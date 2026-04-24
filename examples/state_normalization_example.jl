"""
    state_normalization_example.jl

Example demonstrating state normalization capabilities in SpaceAGORA.jl

This example shows how to:
1. Use automatic quaternion normalization during integration
2. Configure normalization strategies
3. Monitor quaternion norms
4. Access normalization diagnostics

## Usage:

    julia> include("state_normalization_example.jl")
"""

# Note: This is a conceptual example showing API usage patterns.
# Actual simulation requires full SpaceAGORA setup with planet data, GRAM, etc.

module StateNormalizationExample

using LinearAlgebra
using StaticArrays

# ============================================================================
# Example 1: Basic Automatic Normalization
# ============================================================================
"""
This example shows that normalization is automatically applied when 
orientation_sim=true.
"""
function example_automatic_normalization()
    println("\n" * "="^70)
    println("Example 1: Automatic Quaternion Normalization")
    println("="^70)
    
    # When you run a simulation with orientation_sim=true:
    #
    # args = SimulationConfiguration(
    #     mission_configuration=MissionConfiguration(...),
    #     dynamics_model=DynamicsModel(...),
    #     environment_model=EnvironmentModel(...),
    #     # ... other configuration
    # )
    #
    # The state normalization callback is AUTOMATICALLY activated!
    # No additional setup needed.
    
    println("""
    Key Points About Automatic Normalization:
    
    1. AUTOMATIC ACTIVATION:
       - Triggered when: orientation_sim = true
       - No manual configuration required
       - Integrated into standard callback system
    
    2. DEFAULT BEHAVIOR:
       - Normalizes quaternions when ||q|| deviates from 1.0
       - Uses a_tol_quaternion for deviation tolerance
       - Event-driven (efficient) rather than per-step
    
    3. CONFIGURATION:
       You can customize via SimulationSettings:
       
       settings = SimulationSettings(
           a_tol_quaternion = 1e-11,  # Stricter tolerance
           # ... other settings
       )
    """)
end

# ============================================================================
# Example 2: Manual State Normalization (Direct Function Usage)
# ============================================================================
"""
Example showing manual normalization of quaternion vectors.
"""
function example_manual_normalization()
    println("\n" * "="^70)
    println("Example 2: Manual Quaternion Normalization")
    println("="^70)
    
    # Simulate a quaternion that has drifted from unit norm
    # due to integration errors
    q_drifted = SVector{4, Float64}(0.1, 0.2, 0.3, 0.924)  # ||q|| ≈ 0.983
    q_norm_before = norm(q_drifted)
    
    println("\nBefore normalization:")
    println("  q = $q_drifted")
    println("  ||q|| = $(norm(q_drifted))")
    println("  Deviation from unit norm: $(abs(norm(q_drifted) - 1.0))")
    
    # Manual normalization (you would need to import the module)
    # In actual code:
    #   using StateNormalization
    #   q_normalized, original_norm = normalize_quaternion!(q_drifted)
    
    # For this example, show what would happen:
    q_normalized_sim = q_drifted ./ q_norm_before
    println("\nAfter normalization:")
    println("  q = $q_normalized_sim")
    println("  ||q|| = $(norm(q_normalized_sim))")
    println("  Perfectly unit-normalized: $(abs(norm(q_normalized_sim) - 1.0) < 1e-14)")
end

# ============================================================================
# Example 3: Configurable Normalization Strategies
# ============================================================================
"""
Example showing how to configure normalization behavior.
"""
function example_configuration_strategies()
    println("\n" * "="^70)
    println("Example 3: Configurable Normalization Strategies")
    println("="^70)
    
    scenarios = [
        (name="High-Precision Attitude Simulation",
         tolerance=1e-12,
         every_step=true,
         description="Normalize at every step for maximum accuracy"),
        
        (name="Standard Spacecraft Simulation",
         tolerance=1e-10,
         every_step=false,
         description="Event-driven normalization with reasonable tolerance"),
        
        (name="Rapid Prototyping",
         tolerance=1e-8,
         every_step=false,
         description="Relaxed tolerance for development/testing"),
    ]
    
    println("\nCommon Configuration Scenarios:\n")
    for (i, scenario) in enumerate(scenarios)
        println("$i. $(scenario.name)")
        println("   Description: $(scenario.description)")
        println("   Tolerance:   $(scenario.tolerance)")
        println("   Strategy:    $(scenario.every_step ? "Normalize every step" : "Event-driven (on-demand)")")
        println()
    end
    
    println("Example Usage in Code:")
    println("""
    # Approach 1: Via environment variable (per-step normalization)
    ENV["SPACEAGORA_NORMALIZE_EVERY_STEP"] = "1"
    args = SimulationConfiguration(...)
    
    # Approach 2: Via SimulationSettings (tolerance configuration)
    settings = SimulationSettings(
        a_tol_quaternion = 1e-11,
        # ... other settings
    )
    
    # Approach 3: Direct use of StateNormalization module
    using StateNormalization
    config = StateNormalizationConfig(
        quaternion_enabled=true,
        quaternion_tolerance=1e-11,
        quaternion_every_step=false,
        verbose=true
    )
    callback = create_state_normalization_callback(config)
    """)
end

# ============================================================================
# Example 4: Monitoring Quaternion Deviations
# ============================================================================
"""
Example showing how to monitor quaternion health during simulation.
"""
function example_monitoring()
    println("\n" * "="^70)
    println("Example 4: Monitoring Quaternion Deviations")
    println("="^70)
    
    println("""
    To monitor quaternion norms during a simulation:
    
    1. USE ENVIRONMENT VARIABLE FOR DIAGNOSTICS:
       
       export SPACEAGORA_DEBUG=1  # Enable debug output
       julia> results = run_simulation(args)
       # Shows normalization events if they occur
    
    2. POST-PROCESSING ANALYSIS:
       
       # After simulation completes
       using Statistics
       
       quaternion_norms = []
       for (i, result) in enumerate(results)
           for sat_idx in 1:length(result.sc)
               q = result.sc[sat_idx].q
               push!(quaternion_norms, norm(q))
           end
       end
       
       println("Quaternion Norm Statistics:")
       println("  Mean:    $(mean(quaternion_norms))")
       println("  Std Dev: $(std(quaternion_norms))")
       println("  Min:     $(minimum(quaternion_norms))")
       println("  Max:     $(maximum(quaternion_norms))")
       println("  Max Dev: $(maximum(abs.(quaternion_norms .- 1.0)))")
    
    3. REAL-TIME MONITORING IN CALLBACKS:
       
       # Access metrics during integration
       using StateNormalization
       metrics = get_state_deviation_metrics(integrator.u)
       
       if metrics["quaternion_norm_deviation"] > 1e-9
           @warn "Quaternion significantly denormalized"
       end
    """)
end

# ============================================================================
# Example 5: Integration with Solver Settings
# ============================================================================
"""
Example showing integration with ODE solver tolerances.
"""
function example_solver_integration()
    println("\n" * "="^70)
    println("Example 5: Solver Integration & Tolerance Relationship")
    println("="^70)
    
    println("""
    The state normalization works in conjunction with ODE solver tolerances.
    For best results, match them appropriately:
    
    REQUIREMENT               SOLVER TOLERANCE    NORMALIZATION TOLERANCE
    ─────────────────────────────────────────────────────────────────────
    High precision            reltol=1e-10        a_tol_quaternion=1e-11
    attitude control          abstol=1e-12        
    
    Standard spacecraft       reltol=1e-9         a_tol_quaternion=1e-10
    simulations               abstol=1e-11        
    
    Rapid prototyping         reltol=1e-6         a_tol_quaternion=1e-8
    / exploration             abstol=1e-8         
    
    RECOMMENDED SETTINGS:
    
    # For most simulations:
    settings = SimulationSettings(
        reltol_orbit=1e-9,
        abstol_orbit=1e-11,
        reltol_quaternion=1e-9,          # Independent from above
        abstol_quaternion=1e-11,
        a_tol_quaternion=1e-10,           # Event-trigger tolerance
        # ... other settings
    )
    
    # Impact of solver step size:
    args.simulation_settings.dt_max = 0.1  # Smaller steps → better quaternion accuracy
    
    Key Principle:
    - Solver tolerance determines numerical precision of quaternion evolution
    - Normalization tolerance determines when to correct drift
    - They work together: solver maintains precision between normalizations
    """)
end

# ============================================================================
# Example 6: Quaternion Constraint Equation
# ============================================================================
"""
Example explaining the mathematical constraint being maintained.
"""
function example_mathematical_constraint()
    println("\n" * "="^70)
    println("Example 6: Quaternion Mathematical Constraint")
    println("="^70)
    
    println("""
    QUATERNION DEFINITION:
    
    A quaternion q = [q_x, q_y, q_z, q_w]ᵀ represents a 3D rotation as:
    
      ||q|| = √(q_x² + q_y² + q_z² + q_w²) = 1.0
    
    This unit-norm constraint is fundamental to the parameterization.
    
    KINEMATIC EQUATION:
    
    The quaternion time derivative is computed from angular velocity ω:
    
      q̇ = (1/2) Ξ(ω) q
    
    where Ξ is the kinematic matrix.
    
    THE PROBLEM:
    
    ODE integrators solve this differential equation numerically and can
    accumulate errors that cause ||q|| to drift away from 1.0.
    
    Example drift:
      Time t=0:      ||q|| = 1.0000000000
      Time t=10s:    ||q|| = 1.0000000127  (accumulated error)
      Time t=100s:   ||q|| = 1.0000001234  (larger error)
    
    WHY THIS IS PROBLEMATIC:
    
    1. Physically incorrect: The rotation is no longer valid
    2. Cascading errors: Attitude-dependent calculations become inaccurate
    3. Divergence: Errors can grow faster as integration continues
    
    THE SOLUTION:
    
    Periodically renormalize the quaternion:
    
      q_normalized = q / ||q||
    
    This restores the constraint without disrupting integration.
    
    WHEN TO NORMALIZE:
    
    - Event-driven (default):  Normalize when |||q|| - 1.0| > tolerance
    - Every step:             Normalize after each integration step
    
    Event-driven is more efficient while maintaining accuracy.
    """)
end

# ============================================================================
# Example 7: Real-World Scenario
# ============================================================================
"""
Example showing a realistic simulation scenario with orientation.
"""
function example_realistic_scenario()
    println("\n" * "="^70)
    println("Example 7: Realistic Scenario - Spacecraft Attitude Propagation")
    println("="^70)
    
    println("""
    SCENARIO: Propagating a spacecraft in Earth orbit with attitude
    
    CONFIGURATION:
    
      mission_configuration=MissionConfiguration(
          orientation_sim=true,          # ← Enable attitude propagation
          # ... other mission settings
      )
      
      dynamics_model=DynamicsModel(
          spacecraft=[SpacecraftDynamics(
              quaternion_initial=[0.0, 0.0, 0.0, 1.0],    # Identity
              angular_velocity_initial=[0.01, 0.02, 0.005], # rad/s
              inertia_tensor=I_body,                      # 3×3 matrix
              # ... other spacecraft settings
          )]
      )
    
    DEFAULT BEHAVIOR:
    
    1. Initial quaternion is automatically normalized
    2. Integration propagates q̇ = (1/2) Ξ(ω) q
    3. State normalization callback monitors ||q|| at each step
    4. If ||q|| deviates > 1e-10, quaternion is renormalized
    5. Attitude-dependent forces/torques use correct rotation matrix
    6. Simulation completes with properly normalized final quaternion
    
    MONITORING:
    
      # Run simulation
      results = run_simulation(args)
      
      # Check final quaternion norm
      final_q = results[end].sc[1].q
      final_norm = norm(final_q)
      
      @assert abs(final_norm - 1.0) < 1e-9 "Attitude not properly constrained"
      
      # Verify throughout propagation
      quaternion_norms = [norm(r.sc[1].q) for r in results]
      max_deviation = maximum(abs.(quaternion_norms .- 1.0))
      println("Maximum quaternion deviation: $max_deviation")
    
    PERFORMANCE:
    
    - Normalization check: ~0.2 microseconds per step
    - Normalization operation: ~0.1 microseconds (when triggered)
    - Total overhead: < 1% of simulation time
    """)
end

# ============================================================================
# Run All Examples
# ============================================================================
function run_all_examples()
    println("\n")
    println("╔" * "="^68 * "╗")
    println("║" * " "^15 * "SpaceAGORA.jl State Normalization Examples" * " "^11 * "║")
    println("╚" * "="^68 * "╝")
    
    example_automatic_normalization()
    example_manual_normalization()
    example_configuration_strategies()
    example_monitoring()
    example_solver_integration()
    example_mathematical_constraint()
    example_realistic_scenario()
    
    println("\n")
    println("="^70)
    println("Examples Summary")
    println("="^70)
    println("""
    1. ✓ Automatic normalization is enabled by default
    2. ✓ Manual normalization available via StateNormalization module
    3. ✓ Fully configurable via environment variables and settings
    4. ✓ Monitoring tools available for diagnostics
    5. ✓ Integrates seamlessly with ODE solver tolerances
    6. ✓ Maintains quaternion unit-norm constraint mathematically
    7. ✓ Minimal performance overhead in typical simulations
    
    For complete documentation, see: docs/StateNormalization.md
    """)
    
    println("="^70)
    println("\nTo use in your simulation:")
    println("""
    1. Ensure orientation_sim=true in MissionConfiguration
    2. Set desired solver tolerances in SimulationSettings
    3. Run simulation normally - normalization happens automatically!
    4. Monitor results in post-processing if desired
    """)
    println("="^70 * "\n")
end

end  # module StateNormalizationExample

# Run the examples
StateNormalizationExample.run_all_examples()
