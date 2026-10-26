# The study's working example: run the case and print a summary. The workflow
# runs it after the tests. Keep results out of the repository (see .gitignore),
# and record the SpaceAGORA commit of SpaceAGORASource.toml with each one.
using StarterName

solution = run_case()
println("run_case: t = ", last(solution.t), " s, final state finite: ", all(isfinite, solution.u[end]))
