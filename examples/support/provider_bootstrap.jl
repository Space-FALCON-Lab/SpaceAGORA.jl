module ExampleProviderBootstrap

"""Load plotting-compatible core first, then the explicitly selected companion."""
function load!(target::Module, provider::Symbol)
    provider in (:HYPR, :SpaceAGORAEDG) || throw(ArgumentError("Unsupported example provider: $provider"))
    json = get(Base.loaded_modules, Base.PkgId(Base.UUID("682c06a0-de6a-54ab-a142-c8b1cf79cde6"), "JSON"), nothing)
    json === nothing || Base.pkgversion(json) < v"1" || throw(ArgumentError(
        "Plotting examples need JSON 0.21, but JSON 1 is already loaded. Start a fresh Julia process in the selected companion environment."))
    provider_project = Base.active_project()
    original_load_path = copy(LOAD_PATH)
    # common.jl selects examples before SpaceAGORA (and JSON) is imported.
    Base.include(target, joinpath(@__DIR__, "..", "common.jl"))
    try
        if provider_project !== nothing
            directory = dirname(provider_project)
            directory in LOAD_PATH || push!(LOAD_PATH, directory)
        end
        Core.eval(target, Expr(:using, Expr(:., provider)))
    finally
        # The provider is now loaded; keep the caller's load-path isolation.
        empty!(LOAD_PATH)
        append!(LOAD_PATH, original_load_path)
    end
    return nothing
end

end
