# startup.jl – Julia interactive startup template for SpaceAGORA.jl
#
# Copy this file to ~/.julia/config/startup.jl (or append its contents to an
# existing startup.jl) to make interactive work in this project more convenient.
#
# One-time setup:
#   1. Copy or symlink this file:
#        cp startup.jl ~/.julia/config/startup.jl
#
# This template intentionally avoids project-specific package loads so it stays
# usable even in minimal Julia environments.

atreplinit() do repl
    return nothing
end
