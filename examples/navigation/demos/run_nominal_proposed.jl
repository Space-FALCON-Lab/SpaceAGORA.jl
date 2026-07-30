user_arguments = copy(ARGS)
demo_parent = joinpath(tempdir(), "spaceagora_navigation", "demos")
mkpath(demo_parent)
demo_output_root = mktempdir(demo_parent; prefix="proposed_", cleanup=false)

empty!(ARGS)
append!(ARGS, [
    "--mode=single",
    "--mission-time=600",
    "--targets=300",
])
append!(ARGS, user_arguments)
append!(ARGS, [
    "--mode=single",
    "--case=proposed",
    "--output=$(demo_output_root)"
])

include(joinpath(@__DIR__, "..", "..", "run_navigation.jl"))

results_directory = joinpath(demo_output_root, "proposed")
plots_directory = joinpath(demo_output_root, "plots")
animation_path = joinpath(demo_output_root, "animation", "proposed.gif")
repo_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
analysis_script = joinpath(repo_root, "examples", "analyze_navigation.jl")

println()
println("Demo completed")
println("  results: $(results_directory)")
println()
println("Copy and paste to open the numerical results:")
println("code \"$(results_directory)\"")
println()
println("Copy and paste to generate the plots:")
println(
    "julia --project=\"$(repo_root)\" \"$(analysis_script)\" " *
    "--mode=plots --case=proposed --input=\"$(results_directory)\" " *
    "--output=\"$(plots_directory)\""
)
println()
println("After generation, copy and paste to view them:")
println(
    "code \"$(joinpath(plots_directory, "target_scenario.png"))\" " *
    "\"$(joinpath(plots_directory, "observer_initial_configuration_reference_orbits.png"))\" " *
    "\"$(joinpath(plots_directory, "rmse_all_tracks.png"))\" " *
    "\"$(joinpath(plots_directory, "target_errors"))\"/*.png"
)
println()
println("Copy and paste to generate the animation:")
println(
    "julia --project=\"$(repo_root)\" \"$(analysis_script)\" " *
    "--mode=animation --case=proposed --input=\"$(results_directory)\" " *
    "--output=\"$(animation_path)\" --interp-substeps=2 " *
    "--hide-target-trails=true --hide-est-trails=true " *
    "--hide-visibility-spheres=true"
)
println()
println("After rendering, copy and paste to view it:")
println("code \"$(animation_path)\"")
