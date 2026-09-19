# Keep the package's pinned environment first; add script-only dependencies second.
# Usage: julia --project=. scripts/dev/run.jl viewer_demos/build_cygnss_ics.jl
import Pkg
isempty(ARGS) && error("Usage: julia --project=. scripts/dev/run.jl <script relative to scripts/dev> [arguments]")
script = normpath(joinpath(@__DIR__, popfirst!(ARGS)))
isfile(script) || error("No development script at $script")
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io=devnull)
empty!(LOAD_PATH)
append!(LOAD_PATH, ["@", @__DIR__, "@stdlib"])
# Match direct script execution so scripts guarded by PROGRAM_FILE enter main().
original_program_file = Base.PROGRAM_FILE
Base.PROGRAM_FILE = script
try
    include(script)
finally
    Base.PROGRAM_FILE = original_program_file
end
