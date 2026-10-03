"""Availability contract for the optional HYPR companion."""
module HYPRSupport
export HYPRUnavailableError, hypr_available
"""HYPR was selected before loading its optional companion package."""
struct HYPRUnavailableError <: Exception end
Base.showerror(io::IO, ::HYPRUnavailableError) = print(io,
    "HYPR execution requires the optional SpaceAGORAHYPR package. Install it in your project and run `using SpaceAGORAHYPR` before selecting HYPR. See the optional HYPR installation guide.")
const _loaded = Ref(false)
"""Whether the optional HYPR implementation has been loaded in this process."""
hypr_available() = _loaded[]
activate!() = (_loaded[] = true; nothing)
function require_hypr()
    hypr_available() || throw(HYPRUnavailableError())
    return nothing
end
function unavailable(f, args)
    require_hypr()
    throw(MethodError(f, args))
end
end
