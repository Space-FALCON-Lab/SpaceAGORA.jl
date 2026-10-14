module SpaceAGORAVisualizerExt

using SpaceAGORA
import SpaceAGORAVisualizer

contract_version() = SpaceAGORAVisualizer.CONTRACT_VERSION
require_version(version) = SpaceAGORAVisualizer.require_version(version)
viewer_directory() = SpaceAGORAVisualizer.viewer_directory()
function viewer_import_map(directory=nothing)
    directory === nothing && return SpaceAGORAVisualizer.viewer_import_map()
    return SpaceAGORAVisualizer.viewer_import_map(directory)
end
function render_viewer_html(payload; viewer_dir=nothing, kwargs...)
    viewer_dir === nothing && return SpaceAGORAVisualizer.render_viewer_html(payload; kwargs...)
    return SpaceAGORAVisualizer.render_viewer_html(payload; viewer_dir=viewer_dir, kwargs...)
end

function SpaceAGORAVisualizer.visualize(run::SpaceAGORA.MissionRun; kwargs...)
    run.payload_path === nothing && throw(ArgumentError(
        "This mission has no saved viewer payload. Run with save_payload=true or export a payload from its saved scene and results."))
    return SpaceAGORAVisualizer.visualize(run.payload_path; kwargs...)
end

end # module SpaceAGORAVisualizerExt
