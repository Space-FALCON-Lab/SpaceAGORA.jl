# Runs the scripts of the two extension walkthroughs in one process, each in a
# module of its own, as CI does. Each script also runs alone.
const EXAMPLES = ("my_effector.jl", "my_actuator.jl")

for script in EXAMPLES
    println("== ", script)
    Base.include(Module(Symbol("Walkthrough_", first(splitext(script)))), joinpath(@__DIR__, script))
end
