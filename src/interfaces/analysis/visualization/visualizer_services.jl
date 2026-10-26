module VisualizerServices

const CONTRACT_VERSION = 1

"""HTML viewing was requested without the optional SpaceAGORAVisualizer package loaded."""
struct VisualizerUnavailableError <: Exception end
Base.showerror(io::IO, ::VisualizerUnavailableError) = print(io,
    "HTML viewing requires the optional SpaceAGORAVisualizer package. Install it in the same environment and load it with `using SpaceAGORAVisualizer`. Saved scene and viewer JSON export work without it. See the visualizer installation guide.")

_extension() = Base.get_extension(parentmodule(@__MODULE__), :SpaceAGORAVisualizerExt)
available() = (extension = _extension(); extension !== nothing && extension.contract_version() == CONTRACT_VERSION)

function require_available()
    extension = _extension()
    extension === nothing && throw(VisualizerUnavailableError())
    extension.require_version(CONTRACT_VERSION)
    return extension
end

viewer_directory() = require_available().viewer_directory()
viewer_import_map(directory=nothing) = require_available().viewer_import_map(directory)
render_viewer_html(payload; kwargs...) = require_available().render_viewer_html(payload; kwargs...)

end # module VisualizerServices
