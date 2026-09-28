using SpaceAGORA

# Load the optional public wrapper as a package so Julia activates its adapter
# extension. Other CLI commands do not require an atmospheric backend.
if length(ARGS) >= 2 && ARGS[1] == "assets" && ARGS[2] == "check" && length(ARGS) > 2
    Base.find_package("GRAMSuite") === nothing && error(
        "Preset validation needs the public GRAMSuite package. Prepare " *
        "examples/odyssey_surrogate_env/setup.jl and use that project; native GRAM is unnecessary.")
    @eval using GRAMSuite
end

exit(SpaceAGORA.run_cli(copy(ARGS)))
