# Every entry in Project.toml [deps] must be loaded by the package itself.
# Aqua's stale-dependency check catches the class of leftover that #92 and this
# gate's own PR removed by hand (Rotations, ControlSystemsBase, ..., Quaternions):
# packages that stay in the dependency list after their last use is deleted.
#
# The ignore list is deliberate and short. These three are used by the test
# harness, the examples and the documentation build, not by the package, and
# they live in [deps] only because the repository has no separate test/docs
# project yet (an audit item of its own). Do not add package code to this list;
# either use the dependency in src/ or remove it.
using Aqua
using SpaceAGORA

Aqua.test_stale_deps(SpaceAGORA; ignore=[:Aqua, :Documenter, :Plots])
println("stale_deps_gate_ok")
