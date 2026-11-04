using Dates, SHA, TOML

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const LIBRARY_NAME = "9_SpaceAGORA.jl-main_full integration"
const SNAPSHOT = joinpath(@__DIR__, "results", "source_snapshot_pre_absolute_RN")
const SOURCE_PATHS = [
    joinpath(LIBRARY_NAME, name) for name in (
        "src", "ext", "test", "packages", "Project.toml", "Manifest.toml",
        "Artifacts.toml", "CondaPkg.toml", "startup.jl", "LICENSE", "examples/common.jl",
        "data/assets_manifest.toml", "data/odyssey_surrogate_assets.toml",
        "data/surrogate_presets.toml", "data/textures/manifest.toml")
]
push!(SOURCE_PATHS, "Simulation/Code")

filehash(path) = open(io -> bytes2hex(sha256(io)), path)

function selected_hashes(root)
    hashes = Dict{String,String}()
    for relative in SOURCE_PATHS
        path = joinpath(root, relative)
        if isdir(path)
            for (directory, _, files) in walkdir(path), name in files
                file = joinpath(directory, name)
                islink(file) && error("Snapshot needs explicit handling for symlink: $file")
                hashes[relpath(file, root)] = filehash(file)
            end
        else
            isfile(path) || error("Missing snapshot input: $path")
            hashes[relative] = filehash(path)
        end
    end
    hashes
end

function main()
    ispath(SNAPSHOT) && error("Refusing to overwrite source snapshot: $SNAPSHOT")
    before = selected_hashes(PROJECT_ROOT)
    destination = joinpath(SNAPSHOT, "project")
    mkpath(destination)
    for relative in SOURCE_PATHS
        target = joinpath(destination, relative)
        mkpath(dirname(target))
        cp(joinpath(PROJECT_ROOT, relative), target)
    end
    before == selected_hashes(destination) || error("Copied snapshot failed checksum verification")
    before == selected_hashes(PROJECT_ROOT) || error("Sources changed during snapshot creation")
    manifest = Dict("created_utc" => string(now(UTC)), "julia_version" => string(VERSION),
        "source_project_root" => PROJECT_ROOT, "scope" => SOURCE_PATHS,
        "provenance" => "Filesystem snapshot; library is not a Git repository. Includes current local source changes.",
        "external_dependencies" => "Julia package depot and binary/data assets are not copied; versions and asset manifests are preserved.",
        "sha256" => before)
    open(joinpath(SNAPSHOT, "manifest.toml"), "w") do io
        TOML.print(io, manifest; sorted=true)
    end
    archive = joinpath(SNAPSHOT, "sources.tar.gz")
    run(`tar -czf $archive -C $SNAPSHOT project manifest.toml`)
    run(pipeline(`tar -tzf $archive`; stdout=devnull))
    open(joinpath(SNAPSHOT, "sources.tar.gz.sha256"), "w") do io
        println(io, filehash(archive), "  sources.tar.gz")
    end
    println("Verified source snapshot: $SNAPSHOT")
    println("Files: $(length(before)); archive bytes: $(filesize(archive))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
