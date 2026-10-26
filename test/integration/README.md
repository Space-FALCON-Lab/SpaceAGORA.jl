`runtests.jl` loads the package, prepares the shared test scope and runs the
numbered suites in `test/suites/`. It is the integration harness used by
`test/runtests.jl`. The `examples/` entrypoint runs the example regression suite.

Keep orchestration here and shared test support in `test/helpers/`. Create a
domain folder when it has tests to own; do not reserve folders for future work.
See [Test Layout](../README.md) for the suite boundaries and run commands.
