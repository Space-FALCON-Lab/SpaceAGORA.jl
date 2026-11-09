# Using Sensors

Start here to attach an existing sensor, read its measurements and understand
when they update. A sensor measures the vehicle; a navigation model interprets
that measurement; a controller uses the estimate to compute a command.

## Run the gyro and estimator example

Use Julia 1.12 and start in the repository root. These terminal commands work
on Windows, Linux and macOS and need no GRAM or SPICE data:

```text
julia --startup-file=no --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/setup.jl
julia --startup-file=no --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/estimator.jl
```

The script mounts an exact gyro on a spinning cube, passes its measurement
through an estimator to a rate damper, and prints the final rate beside a hand
calculation. It checks the result and runs estimator/controller conformance
checks, including a checkpoint comparison. Read the complete implementation in
[Feedback Controllers and Estimators](gnc_feedback.md#Worked-example-B:-an-estimator-and-a-link).

This is a deliberately simple sensor and controller check. Passing it does not
validate an orbit-determination filter, a flight gyro or a noisy estimator.

## Mounting and frames

The example attaches `RateGyro(period_s=1.0)` with
`add_sensor!(spacecraft, root, RateGyro(period_s=1.0))`. The sensor belongs to
that spacecraft and is mounted on its root link. Its three measurements are
angular rates in rad/s in the mounting link's frame. In this example the root
frame is aligned with the body frame, so the estimator can pass the values to
the body-rate controller directly. A sensor on a rotated link needs the
appropriate frame conversion; do not treat every mounted measurement as a
body-frame vector.

`RateGyro` defaults to zero white noise, zero bias and zero bias random walk.
Its optional `angle_random_walk` is in rad/√s, `initial_bias` in rad/s and
`rate_random_walk` in rad/s/√s. Set these from a justified sensor model before
interpreting a noisy run. The [Sensors API](../generated/sensors_api.md) gives
the constructor arguments for the built-in kinds.

## Sampling, estimation and control are separate clocks

`period_s=1.0` samples the gyro at the start and every second thereafter. The
sensor period is a physical sampling interval, not the integrator's step size.
The estimator and controller have their own periods. At a shared time the
sensor samples first, navigation updates next, and control reads the estimate.
Between samples, the latest measurement is held.

A navigation model reads `sensor_measurement(spacecraft, k)` and
`sensor_measurement_time(spacecraft, k)`, where `k` is the sensor's position in
that spacecraft's sensor list. The time is seconds since the run began; use it
to detect a stale measurement. It is distinct from the latest controller tick.

Change one period at a time in a private copy of the example. A gyro slower
than the controller supplies repeated measurements at some control ticks; a
faster gyro may take samples between estimator ticks. The original hand
calculation assumes the exact gyro samples at the control ticks, so its test
is not an acceptance criterion for a different schedule or noisy sensor.

## Inspect the output

To inspect the existing example interactively, run this terminal command after
setup. The `-i` leaves you at the `julia>` prompt when its checks finish:

```text
julia --startup-file=no --project=examples/gnc_algorithm_env -i examples/gnc_algorithm_env/estimator.jl
```

Then, in the Julia REPL, request a table from a short run:

```julia
result = run_simulation(estimator_configuration(13.0); return_results=true)
result.table[:, [:time, :sc1_sensor_1, :sc1_sensor_2, :sc1_sensor_3]]
```

`time` is in seconds; the three sensor columns are the gyro's x, y and z rates
in rad/s. Columns flatten all sensors' components in mounting order: the
suffix is a component number, not necessarily a sensor number. A second gyro
would occupy the next three columns. The example uses `results=false`, so this
call returns the table in memory without writing a CSV.

Output rows follow `MissionConfiguration.data_rate`, independently of sensor
sampling. A coarse output grid can miss intermediate measurements. The table
contains the latest measurements at its saved times, not a separate log of
every sensor sample. See [Output Columns](output_columns.md#Sensors) for the
other sensor kinds and [Simulation Outputs](outputs.md) for writing tables.

## Reproducibility and next steps

Keep seeds, spacecraft IDs and sensor mounting order fixed when repeating a
run. Those values determine the sensor random streams. A checkpoint preserves
the built-in sensors' random state, held measurements and gyro bias. Set seeds
explicitly when comparing Monte Carlo samples; the default does not vary them
for you.

Before extending the example, record the revision, mounting frame, measurement
units, all three periods (sensor, estimator and controller), output cadence and
seeds. Explain which measurements each controller tick can see. Then compare
an uninterrupted run with a resumed one under identical settings.

To implement a new sensor kind, continue to
[Add a sensor](../extensibility.md#Add-a-sensor). To add a filter or controller,
continue through the [controller walkthrough](gnc_feedback.md#Follow-this-path-first).
