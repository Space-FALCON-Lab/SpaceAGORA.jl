using Test
using Serialization
using TOML

# Exercise the reusable service and the CLI facade without loading the simulator.
module AssetCatalogBoundaryFixture
    const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
    include(joinpath(ROOT, "src", "foundation", "assets", "asset_catalog.jl"))
    module SpaceAGORACLI
        const REPO_ROOT = parentmodule(@__MODULE__).ROOT
        include(joinpath(REPO_ROOT, "src", "interfaces", "cli", "assets.jl"))
    end
end

@testset "Asset catalog service preserves the CLI contract" begin
    catalog = AssetCatalogBoundaryFixture.AssetCatalog
    cli = AssetCatalogBoundaryFixture.SpaceAGORACLI
    @test parentmodule(cli.AssetCheckItem) === cli
    @test parentmodule(cli.AssetCheckReport) === cli
    @test parentmodule(cli.AssetManifestEntry) === cli
    @test parentmodule(cli.check_assets) === cli
    @test parentmodule(cli.load_asset_manifest) === cli
    @test catalog.REPO_ROOT == AssetCatalogBoundaryFixture.ROOT

    mktempdir() do root
        mkpath(joinpath(root, "data", "present"))
        write(joinpath(root, "data", "present.txt"), "fixture")
        write(joinpath(root, "data", "pointer.jls"), "version https://git-lfs.github.com/spec/v1\noid sha256:" * "0"^64 * "\nsize 1\n")
        write(joinpath(root, "data", "zero.bin"), "")
        manifest = joinpath(root, "manifest.toml")
        write(manifest, """
        [[asset]]
        name = "builtin"
        scope = "baseline"
        kind = "builtin"
        [[asset]]
        name = "directory"
        scope = "optional"
        relative_path = "data/present"
        licensing = "open"
        detail = "sample directory"
        [[asset]]
        name = "file"
        scope = "baseline"
        relative_path = "data/present.txt"
        kind = "file"
        required = true
        [[asset]]
        name = "missing"
        scope = "optional"
        relative_path = "data/absent"
        kind = "directory"
        required = true
        [[asset]]
        name = "optional"
        scope = "optional"
        relative_path = "data/absent.txt"
        kind = "file"
        [[asset]]
        name = "pointer"
        scope = "optional"
        relative_path = "data/pointer.jls"
        kind = "file"
        [[asset]]
        name = "zero"
        scope = "optional"
        relative_path = "data/zero.bin"
        kind = "file"
        """)
        before = sort(readdir(root))
        raw = catalog.load_asset_manifest(; repo_root=root, manifest_path=manifest)
        entries = cli.load_asset_manifest(; repo_root=root, manifest_path=manifest)
        @test length(raw) == length(entries) == 7
        @test raw[1] == (name="builtin", scope="baseline", relative_path=".",
            kind="builtin", required=false, licensing="unspecified", detail="")
        @test raw[2].kind == "directory"
        for (record, entry) in zip(raw, entries)
            @test record == NamedTuple{fieldnames(typeof(entry))}(Tuple(getfield(entry, f) for f in fieldnames(typeof(entry))))
            @test cli._manifest_entry_path(root, entry) == catalog.manifest_entry_path(root, record)
            @test cli._manifest_entry_available(root, entry) == catalog.manifest_entry_available(root, record)
            @test cli._manifest_entry_readiness(root, entry) == catalog.manifest_entry_readiness(root, record)
        end
        service_report = catalog.check_assets(; repo_root=root, manifest_path=manifest)
        report = mktempdir() do elsewhere
            cd(elsewhere) do
                cli.check_assets(; repo_root=root, manifest_path=manifest)
            end
        end
        @test report isa cli.AssetCheckReport
        @test report.repo_root == root
        @test getproperty.(report.items, :available) == [true, false, true, false, false, false, false]
        @test getproperty.(report.items, :required) == [false, false, true, true, false, false, false]
        @test getproperty.(report.items, :state) == [:usable, :empty, :usable, :missing, :missing, :pointer_only, :invalid]
        @test getproperty.(report.items, :reason) == getproperty.(service_report.items, :reason)
        @test cli.asset_check_exit_status(report) == 1
        @test report.items[1].path == root
        @test report.items[2].detail == "sample directory [licensing=open]"
        @test report.items[1].detail == "licensing=unspecified"
        for (record, item) in zip(service_report.items, report.items)
            @test record.state === item.state
            @test record.reason == item.reason
            @test record == NamedTuple{fieldnames(typeof(item))}(Tuple(getfield(item, f) for f in fieldnames(typeof(item))))
        end
        @test sort(readdir(root)) == before
        @test !ispath(joinpath(root, "data", "absent"))
        @test !ispath(joinpath(root, "data", "absent.txt"))
        buffer = IOBuffer()
        serialize(buffer, report)
        seekstart(buffer)
        restored = deserialize(buffer)
        @test typeof(restored) === cli.AssetCheckReport
        @test restored.repo_root == report.repo_root
        @test restored.items == report.items
        rendered = sprint(io -> cli.render_asset_report(report; io))
        @test startswith(rendered, "SpaceAGORA asset check\nrepo_root=$root\n")
        @test occursin("- builtin: usable\n", rendered)
        @test occursin("- directory: empty-optional\n", rendered)
        @test occursin("- file: usable\n", rendered)
        @test occursin("- pointer: pointer-only-optional\n", rendered)
        @test occursin("- zero: invalid-optional\n", rendered)
        @test occursin("- missing: missing-required\n", rendered)
        @test occursin("- optional: missing-optional\n", rendered)
        @test occursin("  detail: sample directory [licensing=open]\n", rendered)
        @test occursin("summary: 2 of 7 entries usable; 1 required entries not usable", rendered)

        write(manifest, "")
        @test isempty(catalog.check_assets(; repo_root=root, manifest_path=manifest).items)
        @test isempty(cli.check_assets(; repo_root=root, manifest_path=manifest).items)
        write(manifest, "[[asset]]\nname = \"bad\"\nscope = \"test\"\nkind = \"unknown\"\n")
        @test_throws ArgumentError catalog.check_assets(; repo_root=root, manifest_path=manifest)
        @test_throws ArgumentError cli.check_assets(; repo_root=root, manifest_path=manifest)
        write(manifest, "[[asset]]\nname = \"bad\"\n")
        @test_throws KeyError cli.load_asset_manifest(; repo_root=root, manifest_path=manifest)
    end
end
