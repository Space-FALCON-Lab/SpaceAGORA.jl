# startup.jl – Julia interactive startup template for SpaceAGORA.jl
#
# Copy this file to ~/.julia/config/startup.jl (or append its contents to an
# existing startup.jl).  It loads Revise.jl on every interactive Julia session
# so that changes to source files are hot-reloaded without restarting Julia,
# eliminating the costly package precompilation step during development.
#
# One-time setup:
#   1. Install Revise into your global Julia environment:
#        julia -e 'using Pkg; Pkg.add("Revise")'
#   2. Copy or symlink this file:
#        cp startup.jl ~/.julia/config/startup.jl
#
# After that, start Julia normally (with --project=.SpaceAGORA) and any edits
# you make to .jl files will be picked up automatically the next time you call
# a changed function.

atreplinit() do repl
    try
        @eval using Revise
    catch e
        @warn "Could not load Revise" exception = e
    end
end
