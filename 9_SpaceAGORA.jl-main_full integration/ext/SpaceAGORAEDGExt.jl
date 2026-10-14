module SpaceAGORAEDGExt
using SpaceAGORA
import SpaceAGORAEDG
require_version(version) = SpaceAGORAEDG.require_version(version)
function call_kernel(name::Symbol, args...; kwargs...)
    name in SpaceAGORA.EDGServices.COMPATIBILITY_FUNCTIONS || throw(ArgumentError("Unsupported EDG kernel $name."))
    return getproperty(SpaceAGORAEDG, name)(args...; kwargs...)
end
function __init__()
    require_version(SpaceAGORA.EDGServices.CONTRACT_VERSION)
end
end # module SpaceAGORAEDGExt
