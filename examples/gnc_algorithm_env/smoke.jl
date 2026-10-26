# Runs the worked examples of "Writing a GNC Algorithm" in one process, each in a
# module of its own, as CI does. Each script also runs alone.
const EXAMPLES = ("feedback_controller.jl", "estimator.jl", "rpo_planner.jl", "burn_planner.jl", "joint_planner.jl")

for script in EXAMPLES
    println("== ", script)
    Base.include(Module(Symbol("Example_", first(splitext(script)))), joinpath(@__DIR__, script))
end
