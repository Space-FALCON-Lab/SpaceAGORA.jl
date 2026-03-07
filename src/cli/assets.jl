"""
    AssetCheckItem

Single asset-root status entry produced by `check_assets()`.
"""
Base.@kwdef struct AssetCheckItem
    name::String
    scope::String
    path::String
    available::Bool
    required::Bool
    detail::String = ""
end

"""
    AssetCheckReport

Repository-local asset availability report returned by `check_assets()`.
"""
Base.@kwdef struct AssetCheckReport
    repo_root::String
    items::Vector{AssetCheckItem}
end

@inline function _asset_item(name::String, scope::String, path::String; available::Bool, required::Bool, detail::String="")
    return AssetCheckItem(
        name=name,
        scope=scope,
        path=path,
        available=available,
        required=required,
        detail=detail
    )
end

"""
    check_assets(; repo_root=REPO_ROOT) -> AssetCheckReport

Inspect the standard SpaceAGORA asset roots and return a typed availability report.
"""
function check_assets(; repo_root::String=REPO_ROOT)::AssetCheckReport
    gram_root = joinpath(repo_root, "data", "GRAMSuite.jl", "GRAM Suite 2.0")
    spice_dir = joinpath(gram_root, "SPICE")
    gravity_dir = joinpath(repo_root, "data", "Gravity_harmonics_data")
    topography_dir = joinpath(repo_root, "data", "Topography_harmonics_data")
    surrogate_dir = joinpath(repo_root, "data", "GRAM_surrogate")

    items = AssetCheckItem[]
    push!(items, _asset_item(
        "no_gram_mode",
        "baseline onboarding",
        repo_root;
        available=true,
        required=true,
        detail="Built-in fallback atmosphere models and simple ephemerides are available from the repository."
    ))
    push!(items, _asset_item(
        "gram_root",
        "high-fidelity atmosphere",
        gram_root;
        available=isdir(gram_root),
        required=false,
        detail="Licensed GRAM installation root. Required for GRAM-backed atmospheric studies."
    ))
    push!(items, _asset_item(
        "spice_directory",
        "high-fidelity ephemerides",
        spice_dir;
        available=isdir(spice_dir),
        required=false,
        detail="Machine-local SPICE kernels used by high-fidelity frame, SRP, and N-body workflows."
    ))
    push!(items, _asset_item(
        "gravity_harmonics_directory",
        "gravity harmonics",
        gravity_dir;
        available=isdir(gravity_dir),
        required=false,
        detail="Repository data directory for gravity harmonics CSV inputs."
    ))
    push!(items, _asset_item(
        "topography_harmonics_directory",
        "topography harmonics",
        topography_dir;
        available=isdir(topography_dir),
        required=false,
        detail="Repository data directory for topography harmonics inputs."
    ))
    push!(items, _asset_item(
        "gram_surrogate_directory",
        "GRAM surrogate grids",
        surrogate_dir;
        available=isdir(surrogate_dir),
        required=false,
        detail="Optional surrogate/static-grid bundle for GRAM-accelerated studies."
    ))

    return AssetCheckReport(repo_root=repo_root, items=items)
end

"""
    render_asset_report(report; io=stdout)

Render a human-readable summary of an `AssetCheckReport`.
"""
function render_asset_report(report::AssetCheckReport; io::IO=stdout)
    println(io, "SpaceAGORA asset check")
    println(io, "repo_root=$(report.repo_root)")
    for item in report.items
        status = item.available ? "available" : (item.required ? "missing-required" : "missing-optional")
        println(io, "- $(item.name): $(status)")
        println(io, "  scope: $(item.scope)")
        println(io, "  path: $(item.path)")
        println(io, "  detail: $(item.detail)")
    end
    return nothing
end
