# Full-feature tests explicitly install/load the reviewed renderer. Core-only
# installation probes never include this helper; missing support must not skip.
if !isdefined(@__MODULE__, :SpaceAGORAVisualizer)
    environment = get(ENV, "SPACEAGORA_VISUALIZER_TEST_ENV", "")
    if isempty(environment)
        isdefined(@__MODULE__, :VisualizerInstallation) ||
            include(joinpath(@__DIR__, "..", "..", "scripts", "setup_visualizer.jl"))
        environment = mktempdir(; prefix="spaceagora-visualizer-tests-")
        VisualizerInstallation.setup(environment)
        ENV["SPACEAGORA_VISUALIZER_TEST_ENV"] = environment
    end
    isfile(joinpath(environment, "Project.toml")) || error("Missing prepared visualizer test environment: $environment")
    environment in LOAD_PATH || push!(LOAD_PATH, environment)
end
using SpaceAGORAVisualizer
