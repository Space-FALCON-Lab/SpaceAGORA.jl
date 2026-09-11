# State Normalization in SpaceAGORA.jl

> **Historical document.** The standalone state-normalization module described
> below predates the typed engine, which is SI-native and applies quaternion
> unit-norm projection through the callback system
> (`src/simulation/callbacks/callbacks.jl`) using the utilities in
> `src/core/numerics/quaternion_utils.jl`. The type and function names below
> are retained for reference to the legacy design.

## Overview

The state normalization system in SpaceAGORA.jl provides robust handling of integration state constraints, particularly for quaternion-based attitude representation. This module ensures that rotational states maintain their geometric constraints (e.g., quaternion unit-norm) throughout numerical integration.

## Key Components

### 1. State Normalization Module

**Location:** legacy module (retired); quaternion projection now lives in the callback system

The module provides a comprehensive framework for state normalization with the following key types and functions:

#### Configuration Types

- **`QuaternionNormalizationStrategy`** - Configuration for quaternion normalization behavior
  - `enabled::Bool` - Enable/disable quaternion normalization
  - `tolerance::Float64` - Deviation from unit norm triggering normalization
  - `normalize_every_step::Bool` - Normalize at every step vs. on-demand

- **`StateNormalizationConfig`** - Overall state normalization configuration
  - `quaternion_strategy::QuaternionNormalizationStrategy`
  - `verbose::Bool` - Print normalization events
  - `track_statistics::Bool` - Track normalization statistics

#### Key Functions

```julia
normalize_quaternion!(q::AbstractVector{<:Real}) -> Float64
```
Normalize a quaternion vector to unit norm in-place. Returns the original norm.

```julia
normalize_state_component!(component::AbstractVector{<:Real}) -> Float64
```
Generic normalization for any vector state component.

```julia
normalize_full_state!(state, config::StateNormalizationConfig)
```
Normalize all applicable components of the state according to configuration.

```julia
create_quaternion_normalization_event(config::QuaternionNormalizationStrategy)
```
Create a `DiscreteCallback` for use with DifferentialEquations.jl solvers.

### 2. Integration with Simulation Framework

The state normalization callbacks are automatically integrated into the simulation framework via [src/simulation/callbacks/callbacks.jl](../src/simulation/callbacks/callbacks.jl):

- **`get_state_normalization_callback()`** - Creates the normalization callback
- **Automatic integration** - Callbacks are added to `CallbackSet` in `get_callbacks()`
- **Orientation simulation detection** - Only active when `orientation_sim=true`

## State Vector Structure

The integration state vector contains (per spacecraft):

```
State = [
    pos::SVector{3},           # Position [m]
    vel::SVector{3},           # Velocity [m/s]
    mass::Float64,             # Mass [kg]
    q::SVector{4},             # Quaternion [q_x, q_y, q_z, q_w]  <- normalized
    ω::SVector{3},             # Angular velocity [rad/s]
    heat_loads::Vector{Float64} # Heat per body [J/kg]
]
```

## Quaternion Constraints

A quaternion $\mathbf{q} = [q_x, q_y, q_z, q_w]^T$ must satisfy:

$$\|\mathbf{q}\| = \sqrt{q_x^2 + q_y^2 + q_z^2 + q_w^2} = 1.0$$

### Why Normalization is Required

1. **Numerical Integration Drift**: ODE integrators accumulate small errors when propagating quaternion kinematics
2. **Attitude Representation Validity**: Non-unit quaternions don't represent valid rotations
3. **Physical Correctness**: Violates the fundamental constraint of the attitude parameterization
4. **Error Propagation Prevention**: Unnormalized quaternions can cause misrepresentation in subsequent calculations

### Quaternion Kinematics

Quaternion time derivatives are computed via:

$$\dot{\mathbf{q}} = \frac{1}{2}\boldsymbol{\Xi}(\boldsymbol{\omega}) \mathbf{q}$$

where $\boldsymbol{\Xi}$ is the kinematic matrix relating angular velocity $\boldsymbol{\omega}$ to quaternion rates.

See [src/core/numerics/quaternion_utils.jl](../src/core/numerics/quaternion_utils.jl) for the mathematical details.

## Configuration & Usage

### Basic Usage (Automatic)

When `orientation_sim=true` in your `SimulationConfiguration`, state normalization is **automatically enabled** with default settings:

```julia
args = SimulationConfiguration(
    mission_configuration=MissionConfiguration(...),
    dynamics_model=DynamicsModel(...),
    environment_model=EnvironmentModel(...),
    # ... other settings
)
# Quaternion normalization is automatically active!
```

### Advanced Configuration

Control normalization behavior via **environment variables**:

#### Normalization Frequency

```bash
# Normalize at every integration step (more conservative, slightly slower)
export SPACEAGORA_NORMALIZE_EVERY_STEP=1

# Normalize only when tolerance exceeded (faster, default)
export SPACEAGORA_NORMALIZE_EVERY_STEP=0
```

#### Normalization Tolerance

Configure via `SimulationConfiguration.simulation_settings`:

```julia
settings = SimulationSettings(
    a_tol_quaternion=1e-11,  # Absolute tolerance for quaternion deviation
    # ... other settings
)
```

#### Integration Tolerances

Set solver tolerances for better quaternion accuracy:

```julia
reltol_quaternion = 1e-9
abstol_quaternion = 1e-11
```

These are automatically set in the solver configuration based on `a_tol_quaternion`.

## Normalization Mechanisms

### Automatic Event-Based Normalization

The simulation framework includes a `DiscreteCallback` that:

1. **Detects** quaternion deviations from unit norm
2. **Triggers** when $|\|\mathbf{q}\| - 1.0| > $ `a_tol_quaternion`
3. **Normalizes** by rescaling: $\mathbf{q}_{normalized} = \mathbf{q} / \|\mathbf{q}\|$
4. **Handles degeneracy** by setting to identity quaternion if norm ≈ 0

### Manual Normalization

For custom workflows, use the exported functions:

```julia
using StateNormalization

# Normalize a single quaternion
q = SVector{4, Float64}(0.1, 0.2, 0.3, 0.9)
q_normalized, original_norm = normalize_quaternion!(q)

# Check state deviation metrics
metrics = get_state_deviation_metrics(state)
println("Quaternion norm deviation: ", metrics["quaternion_norm_deviation"])
```

## Solver Integration

The callbacks are automatically added to the `CallbackSet` passed to the ODE solver:

```julia
# In run_simulation.jl, callbacks are built as:
callbacks = SimulationModel.get_callbacks(
    num_sats,
    effectors,
    args;
    saved_values=saved_values,
    save_fields=save_fields
)

# Includes state normalization callback if orientation_sim=true
# Used automatically in: ODEProblem(..., callback=callbacks)
```

## Monitoring & Diagnostics

### Enable Verbose Output

```bash
export SPACEAGORA_DEBUG=1
```

This will print normalization events and statistics.

### Check Statistics

The normalization system can track:

- Number of normalization events
- Maximum quaternion norm deviation
- Timing of normalization calls

Enable via `StateNormalizationConfig`:

```julia
norm_config = StateNormalizationConfig(
    quaternion_enabled=true,
    quaternion_tolerance=1e-10,
    verbose=true,
    track_statistics=true
)
```

## Best Practices

### 1. Tolerance Selection

```julia
# Strict constraints (high-precision attitude propagation)
a_tol_quaternion = 1e-12
reltol_quaternion = 1e-10

# Balanced (typical spacecraft simulations)
a_tol_quaternion = 1e-10
reltol_quaternion = 1e-9

# Relaxed (rapid prototyping, lower priority on attitude)
a_tol_quaternion = 1e-8
reltol_quaternion = 1e-6
```

### 2. Step Size Limits

Reduce max step sizes when attitude is critical:

```julia
args.simulation_settings.dt_max = 0.1  # seconds, during high-precision phases
```

### 3. Periodic Verification

Check quaternion norms in post-processing:

```julia
# After simulation, verify results
for result in results
    q_norm = norm(result.sc[1].q)
    @assert abs(q_norm - 1.0) < 1e-9 "Quaternion not normalized in final state"
end
```

### 4. Initialization

Initial quaternions are automatically normalized in [src/simulation/SimulationElements.jl:1399](src/simulation/SimulationElements.jl#L1399):

```julia
normalize!(in_cond[next_index:next_index+3])
```

Ensure initial conditions provide near-normalized quaternions to minimize solving work.

## Performance Impact

- **Normalization check**: < 1% overhead (simple norm computation)
- **Normalization operation**: < 0.1% overhead (only when needed)
- **Step reduction**: Using smaller step sizes for accuracy may be the dominant cost

**Recommendation**: Use `normalize_every_step=false` (default) and let the event-driven callback normalize only when necessary.

## Advanced Topics

### Custom Normalization Strategies

Extend the module to add specialized normalization:

```julia
struct MyCustomNormalization <: NormalizationStrategy
    # ... fields
end

function create_custom_callback(strategy::MyCustomNormalization)
    # ... implementation
end
```

### Quaternion Verification

Use `get_state_deviation_metrics()` to diagnose quaternion issues:

```julia
metrics = get_state_deviation_metrics(state)
if metrics["quaternion_norm_deviation"] > 1e-8
    @warn "Quaternion significantly denormalized at t=$(state.t)"
end
```

### Integration with Other Constraints

To add angle/angular velocity constraints, extend `StateNormalizationConfig` with additional strategy types.

## References

- Quaternion utilities: [src/core/numerics/quaternion_utils.jl](../src/core/numerics/quaternion_utils.jl)
- Callback system: [src/simulation/callbacks/callbacks.jl](../src/simulation/callbacks/callbacks.jl)
- Simulation setup: [src/simulation/run_simulation.jl](src/simulation/run_simulation.jl)
- State initialization: [src/simulation/SimulationElements.jl](src/simulation/SimulationElements.jl)
- Configuration: [src/core/state/simulation_configuration.jl](../src/core/state/simulation_configuration.jl)

## Troubleshooting

### Quaternion Still Not Normalized

**Issue**: Quaternion norm drifts despite normalization enabled

**Solutions**:
1. Reduce solver tolerance: `reltol=1e-10, abstol=1e-12`
2. Reduce max step size: `dt_max=0.1`
3. Enable `SPACEAGORA_NORMALIZE_EVERY_STEP=1` for conservative approach
4. Check `a_tol_quaternion` is tight enough

### Performance Degradation

**Issue**: Normalization callbacks slow down simulation

**Solutions**:
1. Use `normalize_every_step=false` (event-driven, default)
2. Relax `a_tol_quaternion` if acceptable for application
3. Profile using `SPACEAGORA_DEBUG=1` to find bottlenecks

### Attitude Errors

**Issue**: Attitude propagation accuracy degraded

**Solutions**:
1. Tighten quaternion tolerances
2. Reduce integration step size during critical phases
3. Use `orientation_sim=true` with adequate solver precision
4. Verify initial quaternion is normalized

## Summary

The state normalization system in SpaceAGORA.jl provides:

✓ Automatic quaternion normalization during integration
✓ Configurable tolerance and frequency
✓ Minimal performance impact
✓ Seamless integration with existing simulation framework
✓ Comprehensive diagnostics and monitoring tools

For most users, the default settings provide correct behavior with negligible overhead.
