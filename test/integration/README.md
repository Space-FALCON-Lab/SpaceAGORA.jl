This directory currently hosts the legacy runtime-heavy test harness in `runtests.jl`.

Planned sub-buckets:

- `simulation/`
- `persistence/`
- `cli/`
- `examples/`
- `telemetry/`

Those folders are being introduced first so the suite can be split incrementally without changing the default package test command.
