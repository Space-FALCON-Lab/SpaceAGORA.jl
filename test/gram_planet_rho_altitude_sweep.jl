const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "gram_planet_rho_altitude_sweep.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_sweep()
end
