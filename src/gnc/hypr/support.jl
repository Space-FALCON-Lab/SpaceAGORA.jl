"""Availability contract for the optional HYPR companion."""
module HYPRSupport
export HYPRUnavailableError, HYPRCompatibilityError, hypr_available
"""HYPR was selected before loading its optional companion package."""
struct HYPRUnavailableError <: Exception end
Base.showerror(io::IO, ::HYPRUnavailableError) = print(io,
    "HYPR execution requires the optional SpaceAGORAHYPR package. Install it in your project and run `using SpaceAGORAHYPR` before selecting HYPR. See the optional HYPR installation guide.")
"""An incompatible optional planner implementation was requested in this process."""
struct HYPRCompatibilityError <: Exception
    message::String
end
Base.showerror(io::IO, err::HYPRCompatibilityError) = print(io, err.message)
const _loaded = Ref(false)
const _provider = Ref{Union{Nothing,Tuple{Symbol,VersionNumber}}}(nothing)
const _conflicted = Ref(false)
function check_provider(provider::Symbol, version::VersionNumber)
    _conflicted[] && throw(HYPRCompatibilityError("HYPR loading previously conflicted. Start a fresh process with the supported package pair."))
    current = _provider[]
    if current !== nothing && current != (provider, version)
        _conflicted[] = true
        _loaded[] = false
        throw(HYPRCompatibilityError("Another HYPR implementation is already active. Start a fresh process with only the supported package pair."))
    end
    return nothing
end
function activate!(provider::Symbol, version::VersionNumber)
    check_provider(provider, version)
    _provider[] = (provider, version)
    _loaded[] = true
    return nothing
end
"""Whether the optional HYPR implementation has been loaded in this process."""
hypr_available() = _loaded[] && !_conflicted[]
# The implementation-bearing 0.1 companion used the no-argument entry point.
# Refuse it even when loaded alone; it cannot share this versioned contract.
function activate!()
    _conflicted[] = true
    _loaded[] = false
    throw(HYPRCompatibilityError("SpaceAGORAHYPR 0.1 contains an unsupported implementation. Use the 0.2 compatibility package with HYPR and start a fresh process."))
end
function require_hypr()
    _conflicted[] && throw(HYPRCompatibilityError("HYPR loading conflicted. Start a fresh process with the supported package pair."))
    hypr_available() || throw(HYPRUnavailableError())
    return nothing
end
function unavailable(f, args)
    require_hypr()
    throw(MethodError(f, args))
end
end
