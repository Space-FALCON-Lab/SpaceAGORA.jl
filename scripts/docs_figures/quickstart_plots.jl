# Runs the quickstart example into a temporary folder and copies two of its saved plots
# into docs/src/assets/figures/. Terminal, from the repository root:
#   julia --project=examples scripts/docs_figures/quickstart_plots.jl
root = normpath(joinpath(@__DIR__, "..", ".."))
out = joinpath(root, "docs", "src", "assets", "figures")
mkpath(out)
results = mktempdir()
example = joinpath(root, "examples", "AGORA_Basic_Quickstart.jl")
run(addenv(`$(Base.julia_cmd()) --startup-file=no --project=$(joinpath(root, "examples")) $example`,
           "SPACEAGORA_CLI_OUTPUT_DIR" => results))
for name in ("quickstart_altitude_speed.png", "quickstart_inertial_trajectory.png")
    cp(joinpath(results, "plots", name), joinpath(out, name); force=true)
    println("wrote ", joinpath("docs", "src", "assets", "figures", name))
end
