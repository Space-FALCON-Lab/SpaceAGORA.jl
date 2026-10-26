# Standalone reference-system test orchestration. Raw-load the canonical pure
# numeric owners into the includer before the frame leaf. Do not load SpaceAGORA
# or RuntimeServices here: cold lock-error and includer-local args probes depend
# on the sandbox remaining independent of the package bootstrap.
include(joinpath(@__DIR__, "..", "..", "src", "foundation", "core", "numerics", "geodesy.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "foundation", "core", "numerics", "orbital_geometry.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "foundation", "core", "interfaces", "reference_system.jl"))
