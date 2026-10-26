# The suites resolve default data paths, such as the SPICE kernels under
# data/GRAMSuite.jl, from the repository root, where the direct commands start;
# Pkg.test starts the test process in test/.
cd(dirname(@__DIR__))
include(joinpath(@__DIR__, "integration", "runtests.jl"))
