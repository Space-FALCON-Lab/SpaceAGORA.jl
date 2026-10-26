"""Check (default) or regenerate shared companion support from its canonical core owner.

    julia --startup-file=no scripts/sync_companion_support.jl [--write] <companion>...

No repository is downloaded, no package is loaded and no git operation is run.
"""
module CompanionSupportSync
using SHA, TOML
const SOURCE_ROOT = joinpath(@__DIR__, "companion_support")

function expected_files(repository; source_root=SOURCE_ROOT)
    package = TOML.parsefile(joinpath(repository, "Project.toml"))["name"]
    specification = TOML.parsefile(joinpath(source_root, "manifest.toml"))
    specification["format"] == 1 || error("Unsupported companion support manifest format")
    selected = filter(item -> package in item["packages"], specification["files"])
    isempty(selected) && throw(ArgumentError("No shared support registered for $package"))
    result = Dict{String, String}()
    for item in selected
        source, destination = item["source"], item["destination"]
        # The registry is tracked core input; still keep all reads/writes in its roots.
        for relative in (source, destination)
            isabspath(relative) || ".." in splitpath(relative) ?
                throw(ArgumentError("Support paths must be relative and stay inside their root")) : nothing
        end
        body = read(joinpath(source_root, source), String)
        digest = bytes2hex(sha256(body))
        result[destination] = "# Generated from SpaceAGORA.jl scripts/companion_support/$source.\n" *
            "# Source SHA-256: $digest; regenerate with scripts/sync_companion_support.jl.\n" * body
    end
    return result
end

function sync_support(repository; write=false, source_root=SOURCE_ROOT)
    expected = expected_files(repository; source_root)
    mismatches = String[]
    for (relative, content) in sort!(collect(expected))
        path = joinpath(repository, relative)
        # Refuse symlinks in the destination path, including missing-file parents.
        current = path
        while current != dirname(repository) && current != dirname(current)
            islink(current) && throw(ArgumentError("Refusing a symlinked support destination: $current"))
            current == repository && break
            current = dirname(current)
        end
        isfile(path) && read(path, String) == content && continue
        push!(mismatches, relative)
        if write
            mkpath(dirname(path))
            Base.write(path, content)
        end
    end
    return mismatches
end

function main(args=ARGS)
    write = !isempty(args) && first(args) == "--write"
    repositories = write ? args[2:end] : args
    isempty(repositories) && throw(ArgumentError("Supply one or more companion checkout paths"))
    mismatches = String[]
    for repository in repositories
        root = realpath(repository)
        append!(mismatches, ["$root: $p" for p in sync_support(root; write)])
    end
    if !write && !isempty(mismatches)
        error("Shared companion support differs from core:\n" * join(mismatches, "\n"))
    end
    println(write ? "Regenerated $(length(mismatches)) shared support files." : "Shared companion support matches its canonical source.")
end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && CompanionSupportSync.main()
