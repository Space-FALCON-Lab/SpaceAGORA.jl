# Every entry in Project.toml [deps] must be loaded by the package itself.
# Aqua's stale-dependency check catches the class of leftover that #92 and this
# gate's own PR removed by hand (Rotations, ControlSystemsBase, ..., Quaternions):
# packages that stay in the dependency list after their last use is deleted.
#
# No dependency is exempt. Packages that only the tests, the documentation or the
# examples use belong to their environments (test/Project.toml, docs/Project.toml
# and examples/Project.toml), not to [deps]: either use a dependency in src/ or
# remove it.
#
# Aqua comes from the test environment: run this gate with `julia --project=test`,
# after `julia --startup-file=no test/setup.jl`.
Base.find_package("Aqua") === nothing &&
    error("Run this gate from the test environment: `julia --startup-file=no test/setup.jl`, then `julia --project=test test/gates/ci_stale_deps_gate.jl`.")
using Aqua
using SpaceAGORA

Aqua.test_stale_deps(SpaceAGORA)
println("stale_deps_gate_ok")
