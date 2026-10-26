using Test
using TOML

# The asset check must not call a Git LFS pointer, an empty directory or an
# empty payload "usable". The asset source is loaded into a sandbox module so
# this fixture runs without compiling the package; the file is located by name
# so the fixture follows it across source moves.

function _asset_source(repo::String)
    hits = String[]
    for (root, _, files) in walkdir(joinpath(repo, "src"))
        "assets.jl" in files && push!(hits, joinpath(root, "assets.jl"))
    end
    length(hits) == 1 || error("Expected exactly one assets.jl under src; found $(length(hits)).")
    return only(hits)
end

const _LFS_POINTER = "version https://git-lfs.github.com/spec/v1\noid sha256:" * "0"^64 * "\nsize 149545985\n"

function _write(path::String, text::String)
    mkpath(dirname(path))
    write(path, text)
    return path
end

# The sandbox is built at file top level so the methods it defines are visible
# to the test set below (a later top-level expression) without invokelatest.
const _ASSET_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _ASSET_FIXTURE = mktempdir()
const _ASSET_PARENT = Module(gensym(:AssetCatalogFixture))
Base.include(_ASSET_PARENT, joinpath(_ASSET_REPO, "src", "foundation", "assets", "asset_catalog.jl"))
Core.eval(_ASSET_PARENT, :(module SpaceAGORACLI end))
const _ASSET_SANDBOX = _ASSET_PARENT.SpaceAGORACLI
Core.eval(_ASSET_SANDBOX, :(const REPO_ROOT = $_ASSET_FIXTURE))
Base.include(_ASSET_SANDBOX, _asset_source(_ASSET_REPO))

@testset "Asset recovery guidance across shells" begin
    catalog = _ASSET_PARENT.AssetCatalog
    # Windows has several common shells. Preserve the owner and Git pattern
    # as literal data on every test host, without pretending to run Windows.
    owner = raw"C:\asset checkout's\repo"
    pattern = raw"grid data\owner's\*"
    hint = catalog._lfs_recovery_hint(owner, pattern; windows=true)
    @test occursin("documented asset route", hint)
    @test occursin("owning Git checkout: " * repr(owner), hint)
    @test occursin("Git LFS include pattern: \"grid data/owner's/*\"", hint)
    @test !occursin("git -C", hint)
    @test !occursin("lfs pull", hint)
    @test !occursin("'\"'\"'", hint)
    @test catalog._lfs_recovery_hint("owner", "grids/*") ==
          catalog._lfs_recovery_hint("owner", "grids/*"; windows=Sys.iswindows())

    # The existing POSIX command must still preserve unusual paths exactly.
    if !Sys.iswindows()
        unix_owner = "/asset checkout's/repo"
        unix_pattern = "grid data/owner's/*"
        unix_hint = catalog._lfs_recovery_hint(unix_owner, unix_pattern; windows=false)
        command = only(eachmatch(r"`([^`]+)`", unix_hint)).captures[1]
        arguments = readlines(Cmd(["/bin/sh", "-c", "set -- " * command * "; printf '%s\\n' \"\$@\""]))
        @test arguments == ["git", "-C", unix_owner, "lfs", "pull", "--include=" * unix_pattern]
    end
end

@testset "Asset readiness states" begin
    fixture = _ASSET_FIXTURE
    sandbox = _ASSET_SANDBOX
    try
        mkpath(joinpath(fixture, ".git"))
        # Directory entries.
        mkpath(joinpath(fixture, "empty_dir"))
        _write(joinpath(fixture, "dotfile_dir", ".DS_Store"), "x")
        _write(joinpath(fixture, "hidden_dir", ".cache", "payload.bin"), "ignored")
        _write(joinpath(fixture, "pointer_dir", "a.jls"), _LFS_POINTER)
        _write(joinpath(fixture, "pointer_dir", "b.jls"), _LFS_POINTER)
        _write(joinpath(fixture, "mixed_dir", "a.jls"), _LFS_POINTER)
        _write(joinpath(fixture, "mixed_dir", "nested", "b.jls"), "payload bytes")
        _write(joinpath(fixture, "zero_dir", "a.csv"), "")
        _write(joinpath(fixture, "zero_dir", "b.csv"), "")
        _write(joinpath(fixture, "usable_dir", "nested", "truth.csv"), "t,x\n0,1\n")
        _write(joinpath(fixture, "usable_dir", ".hidden"), "ignored")
        _write(joinpath(fixture, "usable_dir", ".cache", "pointer.jls"), _LFS_POINTER)
        # File entries.
        _write(joinpath(fixture, "empty.feather"), "")
        _write(joinpath(fixture, "pointer.jls"), _LFS_POINTER)
        _write(joinpath(fixture, "usable.stl"), "solid cube\nendsolid cube\n")
        # A short real payload must not be mistaken for a pointer.
        _write(joinpath(fixture, "short.toml"), "v = 1\n")

        entry(name, kind, path; required=false) = Dict(
            "name" => name, "scope" => "test", "kind" => kind, "relative_path" => path,
            "required" => required, "licensing" => "open")
        manifest = joinpath(fixture, "assets_manifest.toml")
        open(manifest, "w") do io
            TOML.print(io, Dict("asset" => [
                entry("builtin_entry", "builtin", "."; required=true),
                entry("missing_dir", "directory", "missing_dir"; required=true),
                entry("empty_dir", "directory", "empty_dir"),
                entry("dotfile_dir", "directory", "dotfile_dir"),
                entry("hidden_dir", "directory", "hidden_dir"),
                entry("pointer_dir", "directory", "pointer_dir"),
                entry("mixed_dir", "directory", "mixed_dir"),
                entry("zero_dir", "directory", "zero_dir"),
                entry("usable_dir", "directory", "usable_dir"; required=true),
                entry("missing_file", "file", "missing.feather"),
                entry("empty_file", "file", "empty.feather"),
                entry("pointer_file", "file", "pointer.jls"),
                entry("usable_file", "file", "usable.stl"; required=true),
                entry("short_file", "file", "short.toml"),
            ]))
        end

        report = sandbox.check_assets(; repo_root=fixture, manifest_path=manifest)
        service_report = _ASSET_PARENT.AssetCatalog.check_assets(; repo_root=fixture, manifest_path=manifest)
        @test getproperty.(service_report.items, :state) == getproperty.(report.items, :state)
        @test getproperty.(service_report.items, :reason) == getproperty.(report.items, :reason)
        states = Dict(item.name => item.state for item in report.items)
        reasons = Dict(item.name => item.reason for item in report.items)
        @test states["builtin_entry"] === :usable
        @test states["missing_dir"] === :missing
        @test states["empty_dir"] === :empty
        @test states["dotfile_dir"] === :empty
        @test states["hidden_dir"] === :empty
        @test states["pointer_dir"] === :pointer_only
        @test states["mixed_dir"] === :pointer_only
        @test states["zero_dir"] === :invalid
        @test states["usable_dir"] === :usable
        @test states["missing_file"] === :missing
        @test states["empty_file"] === :invalid
        @test states["pointer_file"] === :pointer_only
        @test states["usable_file"] === :usable
        @test states["short_file"] === :usable

        # `available` agrees with the state, and every non-usable entry explains itself.
        for item in report.items
            @test item.available == (item.state === :usable)
            @test isempty(item.reason) == (item.state === :usable)
        end
        @test occursin("2 of 2 files are undownloaded Git LFS pointers", reasons["pointer_dir"])
        @test occursin("1 of 2 files are undownloaded Git LFS pointers alongside 1 downloaded payload file", reasons["mixed_dir"])
        if Sys.iswindows()
            @test occursin("owning Git checkout: " * repr(realpath(fixture)), reasons["pointer_file"])
            @test occursin("Git LFS include pattern: \"pointer.jls\"", reasons["pointer_file"])
            @test occursin("Git LFS include pattern: \"pointer_dir/*\"", reasons["pointer_dir"])
        else
            @test occursin("git -C '$(realpath(fixture))' lfs pull --include='pointer.jls'", reasons["pointer_file"])
            @test occursin("git -C '$(realpath(fixture))' lfs pull --include='pointer_dir/*'", reasons["pointer_dir"])
        end
        @test occursin("all 2 files are empty", reasons["zero_dir"])

        # Rendering names the state, adds the reason and ends with a summary.
        text = sprint(io -> sandbox.render_asset_report(report; io=io))
        @test occursin("- missing_dir: missing-required", text)
        @test occursin("- pointer_dir: pointer-only-optional", text)
        @test occursin("- zero_dir: invalid-optional", text)
        @test occursin("- usable_dir: usable", text)
        @test !occursin(": available", text)
        @test occursin("  reason: the directory contains no files", text)
        @test occursin("summary: 4 of 14 entries usable; 1 required entries not usable", text)

        # Exit status follows required entries only.
        @test sandbox.asset_check_exit_status(report) == 1
        open(manifest, "w") do io
            TOML.print(io, Dict("asset" => [
                entry("builtin_entry", "builtin", "."; required=true),
                entry("usable_dir", "directory", "usable_dir"; required=true),
                entry("pointer_dir", "directory", "pointer_dir"),
            ]))
        end
        clean = sandbox.check_assets(; repo_root=fixture, manifest_path=manifest)
        @test sandbox.asset_check_exit_status(clean) == 0
        @test any(item -> item.name == "pointer_dir" && item.state === :pointer_only, clean.items)

        # Unknown kinds and unknown states are rejected, not guessed.
        open(manifest, "w") do io
            TOML.print(io, Dict("asset" => [entry("odd", "archive", "usable.stl")]))
        end
        @test_throws ArgumentError sandbox.check_assets(; repo_root=fixture, manifest_path=manifest)
        @test_throws ArgumentError sandbox._asset_item("x", "test", fixture; state=:present, required=false)
        @test_throws ArgumentError sandbox.AssetCheckItem("x", "test", fixture, false, true, "", :present, "")
        @test_throws ArgumentError sandbox.AssetCheckItem(name="x", scope="test", path=fixture, available=true, required=true, state=:missing)

        # Keyword construction without a state still works for existing callers.
        legacy = sandbox.AssetCheckItem(name="x", scope="test", path=fixture, available=false, required=true)
        @test legacy.state === :missing
        @test isempty(legacy.reason)
        for available in (false, true)
            positional = sandbox.AssetCheckItem("old", "test", fixture, available, true, "legacy detail")
            @test positional.available == available
            @test positional.state === (available ? :usable : :missing)
            @test positional.required
            @test positional.detail == "legacy detail"
            @test isempty(positional.reason)
        end

        # A nested checkout/submodule owns its pointers, even when the command
        # inspects the parent repository. Its .git marker is a file, as it is in
        # linked worktrees; the apostrophe also exercises shell argument quoting.
        vendor = joinpath(fixture, "vendor checkout's")
        _write(joinpath(vendor, ".git"), "gitdir: ../metadata\n")
        _write(joinpath(vendor, "grids", "pointer.jls"), _LFS_POINTER)
        open(manifest, "w") do io
            TOML.print(io, Dict("asset" => [entry("vendor", "directory", "vendor checkout's/grids")]))
        end
        vendor_item = only(sandbox.check_assets(; repo_root=fixture, manifest_path=manifest).items)
        @test vendor_item.state === :pointer_only
        vendor_command = match(r"`([^`]+)`", vendor_item.reason)
        if Sys.iswindows()
            @test vendor_command === nothing
            @test occursin("owning Git checkout: " * repr(realpath(vendor)), vendor_item.reason)
            @test occursin("Git LFS include pattern: \"grids/*\"", vendor_item.reason)
        else
            @test vendor_command !== nothing
        end
        if vendor_command !== nothing && !Sys.iswindows()
            # Parse, but never execute, the suggested git command. Positional
            # shell arguments establish the exact checkout and include pattern.
            parse_command = "set -- " * vendor_command.captures[1] * "; printf '%s\\n' \"\$@\""
            arguments = readlines(Cmd(Cmd(["/bin/sh", "-c", parse_command]); dir=dirname(fixture)))
            @test arguments == ["git", "-C", realpath(vendor), "lfs", "pull", "--include=grids/*"]
        end

        mktempdir() do detached
            _write(joinpath(detached, "pointer.jls"), _LFS_POINTER)
            detached_manifest = joinpath(detached, "assets_manifest.toml")
            open(detached_manifest, "w") do io
                TOML.print(io, Dict("asset" => [entry("detached", "file", "pointer.jls")]))
            end
            detached_item = only(sandbox.check_assets(; repo_root=detached, manifest_path=detached_manifest).items)
            @test detached_item.state === :pointer_only
            @test occursin("no Git checkout found", detached_item.reason)
            @test !occursin("git -C", detached_item.reason)
        end

        # File symlinks classify their target; nested directory symlinks do
        # not add payloads or pointers from outside the inspected directory.
        if !Sys.iswindows()
            links = joinpath(fixture, "links")
            mkpath(links)
            symlink(joinpath(fixture, "usable.stl"), joinpath(links, "payload.stl"))
            symlink(joinpath(fixture, "pointer.jls"), joinpath(links, "pointer.jls"))
            symlink(joinpath(vendor, "grids", "pointer.jls"), joinpath(links, "vendor_pointer.jls"))
            symlink(joinpath(fixture, "absent"), joinpath(links, "broken"))
            mkpath(joinpath(links, "only_directory_link"))
            symlink(joinpath(fixture, "usable_dir"), joinpath(links, "only_directory_link", "payloads"))
            symlink(joinpath(fixture, "pointer_dir"), joinpath(fixture, "usable_dir", "pointer_link"))
            open(manifest, "w") do io
                TOML.print(io, Dict("asset" => [
                    entry("linked_payload", "file", "links/payload.stl"),
                    entry("linked_pointer", "file", "links/pointer.jls"),
                    entry("broken_link", "file", "links/broken"),
                    entry("directory_link_only", "directory", "links/only_directory_link"),
                    entry("usable_with_directory_link", "directory", "usable_dir"),
                    entry("cross_checkout_links", "directory", "links"),
                ]))
            end
            linked = Dict(item.name => item for item in sandbox.check_assets(; repo_root=fixture, manifest_path=manifest).items)
            @test linked["linked_payload"].state === :usable
            @test linked["linked_pointer"].state === :pointer_only
            @test linked["broken_link"].state === :missing
            @test linked["directory_link_only"].state === :empty
            @test linked["usable_with_directory_link"].state === :usable
            @test linked["cross_checkout_links"].state === :pointer_only
            linked_commands = [m.captures[1] for m in eachmatch(r"`([^`]+)`", linked["cross_checkout_links"].reason)]
            @test length(linked_commands) == 2
            linked_arguments = [readlines(Cmd(["/bin/sh", "-c", "set -- " * command * "; printf '%s\\n' \"\$@\""]))
                                for command in linked_commands]
            @test ["git", "-C", realpath(fixture), "lfs", "pull", "--include=pointer.jls"] in linked_arguments
            @test ["git", "-C", realpath(vendor), "lfs", "pull", "--include=grids/pointer.jls"] in linked_arguments
        end
    finally
        rm(fixture; recursive=true, force=true)
    end
end

# Only the unrelated environment-model imports are stubbed. The asset logic,
# CLI dispatcher, and both process entrypoints below are unchanged production
# files. Preserving their relative locations also exercises repo-root discovery.
function _asset_entrypoint_fixture(fixture::String)
    asset_source = _asset_source(_ASSET_REPO)
    cli_dir = dirname(asset_source)
    catalog_source = joinpath(_ASSET_REPO, "src", "foundation", "assets", "asset_catalog.jl")
    for source in (catalog_source, asset_source, joinpath(cli_dir, "spaceagora_cli.jl"),
                   joinpath(cli_dir, "main.jl"), joinpath(_ASSET_REPO, "scripts", "assets", "check_assets.jl"))
        destination = joinpath(fixture, relpath(source, _ASSET_REPO))
        mkpath(dirname(destination))
        cp(source, destination)
    end
    cli_include = joinpath(relpath(cli_dir, joinpath(_ASSET_REPO, "src")), "spaceagora_cli.jl")
    _write(joinpath(fixture, "src", "SpaceAGORA.jl"), """
        module SpaceAGORA
        module SimulationModel
        module EnvironmentModels
        function available_surrogate_presets end
        function resolve_surrogate_preset end
        function surrogate_preset_model end
        function atmosphere_provenance end
        end
        end
        include("foundation/assets/asset_catalog.jl")
        include($(repr(cli_include)))
        using .SpaceAGORACLI: check_assets, render_asset_report, run_cli
        end
        """)
    open(joinpath(fixture, "Project.toml"), "w") do io
        TOML.print(io, Dict("name" => "SpaceAGORA",
                           "uuid" => "afbfb69f-5c0b-4832-b760-43725dff8540",
                           "deps" => Dict("TOML" => string(Base.PkgId(TOML).uuid))))
    end
    return joinpath(fixture, relpath(joinpath(cli_dir, "main.jl"), _ASSET_REPO)),
           joinpath(fixture, "scripts", "assets", "check_assets.jl")
end

@testset "Asset check process exit status" begin
    mktempdir() do fixture
        cli_main, asset_script = _asset_entrypoint_fixture(fixture)
        _write(joinpath(fixture, "pointer.jls"), _LFS_POINTER)
        _write(joinpath(fixture, "invalid.bin"), "")
        _write(joinpath(fixture, "payload.bin"), "payload bytes")
        mkpath(joinpath(fixture, "empty_dir"))
        mkpath(joinpath(fixture, "data"))
        mkpath(joinpath(fixture, ".git"))
        entry(name, kind, path, required) = Dict(
            "name" => name, "scope" => "test", "kind" => kind,
            "relative_path" => path, "required" => required, "licensing" => "open")
        unavailable = [
            ("missing", "directory", "absent", "missing"),
            ("pointer", "file", "pointer.jls", "pointer-only"),
            ("empty", "directory", "empty_dir", "empty"),
            ("invalid", "file", "invalid.bin", "invalid"),
        ]
        scenarios = [(name=name, entries=[entry(name, kind, path, true)], code=1,
                      labels=["- $(name): $(state)-required"])
                     for (name, kind, path, state) in unavailable]
        push!(scenarios, (name="optional entries",
                          entries=vcat([entry("builtin", "builtin", ".", true)],
                                       [entry(name, kind, path, false) for (name, kind, path, _) in unavailable]),
                          code=0, labels=["- $(name): $(state)-optional" for (name, _, _, state) in unavailable]))
        push!(scenarios, (name="required payload", entries=[entry("payload", "file", "payload.bin", true)],
                          code=0, labels=["- payload: usable"]))
        for scenario in scenarios
            open(joinpath(fixture, "data", "assets_manifest.toml"), "w") do io
                TOML.print(io, Dict("asset" => scenario.entries))
            end
            for (name, script, args) in (("CLI", cli_main, ["assets", "check"]),
                                        ("asset script", asset_script, String[]))
                @testset "$(name): $(scenario.name)" begin
                    output = IOBuffer()
                    command = Cmd(`$(Base.julia_cmd()) --startup-file=no --history-file=no --compiled-modules=no --project=$fixture $script $args`;
                                  dir=dirname(fixture))
                    command = addenv(command, "JULIA_LOAD_PATH" => "@:@stdlib", "JULIA_PKG_PRECOMPILE_AUTO" => "0")
                    process = run(pipeline(ignorestatus(command); stdout=output, stderr=output))
                    text = String(take!(output))
                    @test process.exitcode == scenario.code
                    for label in scenario.labels
                        @test occursin(label, text)
                    end
                    @test occursin("SpaceAGORA asset check", text)
                    @test occursin("summary: ", text)
                    @test !occursin("ERROR:", text)
                end
            end
        end
    end
end
