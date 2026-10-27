# Canonical migration of the retired SpaceAGORAHYPR loader.
module LegacyEnvironmentMigration
using TOML
function retire_loader_manifest_data!(data)
    entries = get(data, "deps", data)
    for (name, records) in entries, record in records
        dependencies = get(record, "deps", String[])
        depends_on_loader = dependencies isa AbstractDict ?
            haskey(dependencies, "SpaceAGORAHYPR") : "SpaceAGORAHYPR" in dependencies
        depends_on_loader && error("$name still depends on the retired SpaceAGORAHYPR loader. Migrate that package to direct SpaceAGORA and HYPR imports before rerunning setup.")
        weak = get(record, "weakdeps", nothing)
        if weak isa AbstractVector && "SpaceAGORAHYPR" in weak
            # Keep the retired extension binding even if a stale manifest has
            # already lost its package record; it must remain optional.
            record["weakdeps"] = Dict(dependency =>
                dependency == "SpaceAGORAHYPR" ? "4096c904-79d0-4d9a-9f5d-f1d6c7a2f1b4" : only(entries[dependency])["uuid"]
                for dependency in weak)
        end
    end
    pop!(entries, "SpaceAGORAHYPR", nothing)
    return nothing
end

function retire_loader_manifest!(path)
    data = TOML.parsefile(path)
    entries = get(data, "deps", data)
    before = deepcopy(data)
    retire_loader_manifest_data!(data)
    data == before && return nothing
    open(io -> TOML.print(io, data; sorted=true), path, "w")
    return nothing
end

function retire_loader_project!(project)
    for section in ("deps", "compat", "sources")
        haskey(project, section) && pop!(project[section], "SpaceAGORAHYPR", nothing)
    end
    return project
end

# Read and validate both records before callers write either one. The original
# environment is never changed by this planning operation.
function migrated_records(project_path, manifest_path)
    project = isfile(project_path) ? TOML.parsefile(project_path) : Dict{String,Any}()
    manifest = isfile(manifest_path) ? TOML.parsefile(manifest_path) : nothing
    manifest === nothing || retire_loader_manifest_data!(manifest)
    retire_loader_project!(project)
    return (; project, manifest)
end

end
