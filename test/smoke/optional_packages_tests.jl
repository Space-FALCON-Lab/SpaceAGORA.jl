# Bounded packaging tests: local Git fixtures only; no Pkg, native library or core import.
using Test, TOML, SHA
include(joinpath(@__DIR__, "..", "..", "scripts", "optional_packages.jl"))

const OPTIONAL_FIXTURE_PACKAGES = (
    hypr=(name="HYPR", path="packages/robotics/HYPR.jl",
          uuid="342c8629-b435-4bda-9d24-004c735fb297", version="0.1.1"),
    mujoco=(name="SpaceAGORAMuJoCo", path="packages/robotics/SpaceAGORAMuJoCo.jl",
            uuid="a46e04e1-adb7-43c3-a03b-26575bf4fb67", version="0.1.0"),
    edg=(name="SpaceAGORAEDG", path="packages/gnc/SpaceAGORAEDG.jl",
         uuid="a782e8ee-96d7-4dd6-8231-bb64649e6ad4", version="0.1.0"),
)

# Automatic maintenance after a commit or fetch runs in the background and can
# repack loose objects while `fixture_snapshot` walks `.git`, so it is disabled.
function fixture_git(args...)
    addenv(Cmd(["git", "-c", "init.defaultBranch=main", "-c", "protocol.file.allow=always",
                "-c", "user.name=Optional package fixture", "-c", "user.email=fixture@example.invalid",
                "-c", "core.autocrlf=false", "-c", "maintenance.auto=false", "-c", "gc.auto=0",
                String.(args)...]),
           "GIT_DIR" => nothing, "GIT_WORK_TREE" => nothing, "GIT_INDEX_FILE" => nothing,
           "GIT_COMMON_DIR" => nothing, "GIT_OBJECT_DIRECTORY" => nothing,
           "GIT_ALTERNATE_OBJECT_DIRECTORIES" => nothing, "GIT_CONFIG_NOSYSTEM" => "1",
           "GIT_CONFIG_GLOBAL" => (Sys.iswindows() ? "NUL" : "/dev/null"))
end
fixture_run(args...) = run(pipeline(fixture_git(args...); stdout=devnull, stderr=devnull))
fixture_read(args...) = readchomp(fixture_git(args...))

function fixture_toml(path, data)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
end

function optional_fixture(f)
    mktempdir(; prefix="spaceagora-optional-fixture-") do temporary
        temporary = realpath(temporary)
        root = joinpath(temporary, "core")
        mkpath(root)
        fixture_run("init", "--quiet", root)
        fixture_toml(joinpath(root, "Project.toml"), Dict(
            "name" => "SpaceAGORA", "uuid" => "afbfb69f-5c0b-4832-b760-43725dff8540",
            "version" => "0.2.0", "deps" => Dict{String,String}()))
        write(joinpath(root, "Manifest.toml"), "# Core manifest sentinel; no installation is allowed.\n")
        mkpath(joinpath(root, "scripts"))
        write(joinpath(root, "scripts", "setup_hypr.jl"), "error(\"Fixture installer must never execute\")\n")
        registry = Dict{String,Any}()
        origins = Dict{String,String}()
        revisions = Dict{String,String}()
        next_revisions = Dict{String,String}()
        modules = IOBuffer()
        for key in (:hypr, :mujoco, :edg)
            package = getproperty(OPTIONAL_FIXTURE_PACKAGES, key)
            label = String(key)
            origin = joinpath(temporary, label * "-origin")
            mkpath(origin)
            fixture_run("init", "--quiet", origin)
            fixture_toml(joinpath(origin, "Project.toml"), Dict(
                "name" => package.name, "uuid" => package.uuid, "version" => package.version))
            mkpath(joinpath(origin, "scripts"))
            write(joinpath(origin, "scripts", "setup_env.jl"), "error(\"Fixture installer must never execute\")\n")
            fixture_run("-C", origin, "add", ".")
            fixture_run("-C", origin, "commit", "--quiet", "-m", "Pinned fixture")
            revision = fixture_read("-C", origin, "rev-parse", "HEAD")
            write(joinpath(origin, "later.txt"), "A later commit must not be selected implicitly.\n")
            fixture_run("-C", origin, "add", "later.txt")
            fixture_run("-C", origin, "commit", "--quiet", "-m", "Unpinned fixture revision")
            origins[label] = origin
            revisions[label] = revision
            next_revisions[label] = fixture_read("-C", origin, "rev-parse", "HEAD")
            registry[label] = Dict("name" => package.name, "path" => package.path,
                                   "uuid" => package.uuid, "version" => package.version, "url" => origin)
            println(modules, "[submodule \"", package.name, "\"]\n\tpath = ", package.path,
                    "\n\turl = ", origin, "\n\tfetchRecurseSubmodules = false")
            fixture_run("-C", root, "update-index", "--add", "--cacheinfo", "160000,$revision,$(package.path)")
        end
        # This unrelated submodule deliberately has no usable remote. Selected
        # checkout must never initialize it or contact its URL.
        println(modules, "[submodule \"GRAMSuite.jl\"]\n\tpath = data/GRAMSuite.jl\n\turl = ",
                joinpath(temporary, "absent-gram-origin"))
        fixture_run("-C", root, "update-index", "--add", "--cacheinfo",
                    "160000,$(revisions["hypr"]),data/GRAMSuite.jl")
        write(joinpath(root, ".gitmodules"), String(take!(modules)))
        fixture_toml(joinpath(root, "packages", "OptionalPackages.toml"), registry)
        fixture_toml(joinpath(root, "packages", "HYPRSource.toml"), Dict(
            "url" => origins["hypr"], "rev" => revisions["hypr"],
            "uuid" => OPTIONAL_FIXTURE_PACKAGES.hypr.uuid, "version" => "0.1.1"))
        fixture_run("-C", root, "add", ".gitmodules", "Project.toml", "Manifest.toml", "scripts",
                    "packages/OptionalPackages.toml", "packages/HYPRSource.toml")
        fixture_run("-C", root, "commit", "--quiet", "-m", "Fixture core and optional gitlinks")
        # Local file transport is enabled only for these disposable repositories.
        withenv("GIT_CONFIG_COUNT" => "1", "GIT_CONFIG_KEY_0" => "protocol.file.allow",
                "GIT_CONFIG_VALUE_0" => "always", "GIT_CONFIG_NOSYSTEM" => "1",
                "GIT_CONFIG_GLOBAL" => (Sys.iswindows() ? "NUL" : "/dev/null")) do
            f((; root, temporary, origins, revisions, next_revisions))
        end
    end
end

function fixture_snapshot(root)
    entries = Dict{String,Any}()
    for (directory, dirs, files) in walkdir(root; follow_symlinks=false)
        for name in vcat(dirs, files)
            path = joinpath(directory, name)
            entries[relpath(path, root)] = islink(path) ? (:link, readlink(path)) :
                isdir(path) ? (:directory,) : (:file, bytes2hex(sha256(read(path))))
        end
    end
    entries
end

function assert_refused_unchanged(fixture, name)
    before = fixture_snapshot(fixture.root)
    @test OptionalPackages.status(name; root=fixture.root).state == :invalid
    @test_throws Exception OptionalPackages.check(name; root=fixture.root)
    @test_throws Exception OptionalPackages.checkout(name; root=fixture.root)
    @test fixture_snapshot(fixture.root) == before
end

@testset "Optional package metadata and source archives" begin
    optional_fixture() do fixture
        for name in ("hypr", "mujoco", "edg")
            expected = getproperty(OPTIONAL_FIXTURE_PACKAGES, Symbol(name))
            resolved = OptionalPackages.spec(name; root=fixture.root)
            @test resolved.name == expected.name
            @test resolved.path == expected.path
            @test string(resolved.uuid) == expected.uuid
            @test string(resolved.version) == expected.version
            @test resolved.url == fixture.origins[name]
            @test resolved.revision == fixture.revisions[name]
            before = fixture_snapshot(fixture.root)
            @test OptionalPackages.status(name; root=fixture.root).state == :absent
            @test_throws Exception OptionalPackages.check(name; root=fixture.root)
            @test fixture_snapshot(fixture.root) == before
        end
        for unknown in ("all", "../hypr", "")
            @test_throws Exception OptionalPackages.spec(unknown; root=fixture.root)
        end
        archive = joinpath(fixture.temporary, "source-zip")
        mkpath(joinpath(archive, "packages"))
        cp(joinpath(fixture.root, "packages", "OptionalPackages.toml"),
           joinpath(archive, "packages", "OptionalPackages.toml"))
        before = fixture_snapshot(archive)
        @test_throws Exception OptionalPackages.spec("hypr"; root=archive)
        @test fixture_snapshot(archive) == before
    end
end

@testset "Only the selected pinned repository is initialized" begin
    for name in ("hypr", "mujoco", "edg")
        optional_fixture() do fixture
            package = getproperty(OPTIONAL_FIXTURE_PACKAGES, Symbol(name))
            other = name == "hypr" ? "mujoco" : "hypr"
            path = joinpath(fixture.root, package.path)
            @test OptionalPackages.checkout(name; root=fixture.root) == path
            @test isfile(joinpath(path, ".git"))
            @test fixture_read("-C", path, "rev-parse", "HEAD") == fixture.revisions[name]
            @test OptionalPackages.status(name; root=fixture.root).state == :ready
            @test OptionalPackages.check(name; root=fixture.root) == path
            @test OptionalPackages.status(other; root=fixture.root).state == :absent
            @test !ispath(joinpath(fixture.root, "data", "GRAMSuite.jl", "Project.toml"))
            @test !ispath(joinpath(fixture.root, ".git", "modules", "GRAMSuite.jl"))
            before = fixture_snapshot(fixture.root)
            @test OptionalPackages.checkout(name; root=fixture.root) == path
            @test fixture_snapshot(fixture.root) == before
        end
    end
end

@testset "Existing checkouts and local work are never overwritten" begin
    for change in (:tracked, :staged, :untracked, :head, :origin, :uuid, :name, :version)
        optional_fixture() do fixture
            path = OptionalPackages.checkout("mujoco"; root=fixture.root)
            project_path = joinpath(path, "Project.toml")
            if change == :tracked || change == :staged
                open(project_path, "a") do io
                    println(io, "# User's local change")
                end
                change == :staged && fixture_run("-C", path, "add", "Project.toml")
            elseif change == :untracked
                write(joinpath(path, "user-notes.txt"), "Keep my work.\n")
            elseif change == :head
                fixture_run("-C", path, "fetch", "--quiet", "origin", fixture.next_revisions["mujoco"])
                fixture_run("-C", path, "checkout", "--quiet", "--detach", fixture.next_revisions["mujoco"])
            elseif change == :origin
                fixture_run("-C", path, "remote", "set-url", "origin", fixture.origins["hypr"])
            else
                # Commit the wrong metadata and update the fixture gitlink so the
                # package identity itself, rather than dirt or HEAD, must fail.
                project = TOML.parsefile(project_path)
                project[String(change)] = change == :uuid ? "00000000-0000-0000-0000-000000000000" :
                    change == :name ? "UnrelatedPackage" : "9.9.9"
                fixture_toml(project_path, project)
                fixture_run("-C", path, "add", "Project.toml")
                fixture_run("-C", path, "commit", "--quiet", "-m", "Wrong package metadata")
                revision = fixture_read("-C", path, "rev-parse", "HEAD")
                fixture_run("-C", fixture.root, "update-index", "--cacheinfo",
                    "160000,$revision,$(OPTIONAL_FIXTURE_PACKAGES.mujoco.path)")
            end
            assert_refused_unchanged(fixture, "mujoco")
        end
    end
end

@testset "Pin and Git metadata must agree" begin
    for change in (:hypr_pin, :module_url, :module_path, :not_gitlink)
        optional_fixture() do fixture
            if change == :hypr_pin
                pin_path = joinpath(fixture.root, "packages", "HYPRSource.toml")
                pin = TOML.parsefile(pin_path)
                pin["rev"] = fixture.next_revisions["hypr"]
                fixture_toml(pin_path, pin)
            elseif change == :module_url
                fixture_run("-C", fixture.root, "config", "--file", ".gitmodules",
                            "submodule.HYPR.url", fixture.origins["mujoco"])
            elseif change == :module_path
                fixture_run("-C", fixture.root, "config", "--file", ".gitmodules",
                            "submodule.HYPR.path", "packages/robotics/AnotherPackage.jl")
            else
                fixture_run("-C", fixture.root, "update-index", "--force-remove", OPTIONAL_FIXTURE_PACKAGES.hypr.path)
            end
            before = fixture_snapshot(fixture.root)
            @test_throws Exception OptionalPackages.spec("hypr"; root=fixture.root)
            @test fixture_snapshot(fixture.root) == before
        end
    end
end

@testset "Unmanaged and symlinked targets are refused" begin
    for hazard in (:nongit, :file, :symlink, :parent_symlink, :git_symlink)
        optional_fixture() do fixture
            path = joinpath(fixture.root, OPTIONAL_FIXTURE_PACKAGES.mujoco.path)
            mkpath(dirname(path))
            outside = joinpath(fixture.temporary, "outside")
            mkpath(outside)
            write(joinpath(outside, "keep.txt"), "Do not alter files outside this source tree.\n")
            if hazard == :nongit
                mkpath(path)
                write(joinpath(path, "keep.txt"), "Not a package checkout.\n")
            elseif hazard == :file
                write(path, "A user's file occupies the target.\n")
            elseif hazard == :symlink
                symlink(outside, path; dir_target=true)
            elseif hazard == :parent_symlink
                rm(dirname(path)) # newly created and empty fixture directory only
                symlink(outside, dirname(path); dir_target=true)
            else
                OptionalPackages.checkout("mujoco"; root=fixture.root)
                gitfile = joinpath(path, ".git")
                saved = joinpath(outside, "gitfile")
                mv(gitfile, saved)
                symlink(saved, gitfile)
            end
            outside_before = fixture_snapshot(outside)
            assert_refused_unchanged(fixture, "mujoco")
            @test fixture_snapshot(outside) == outside_before
        end
    end
    optional_fixture() do fixture
        fixture_run("-C", fixture.root, "config", "submodule.SpaceAGORAMuJoCo.url", fixture.origins["hypr"])
        before = fixture_snapshot(fixture.root)
        @test_throws Exception OptionalPackages.checkout("mujoco"; root=fixture.root)
        @test fixture_snapshot(fixture.root) == before
    end
end

@testset "Installation command routing and destination safety" begin
    optional_fixture() do fixture
        hypr = OptionalPackages.checkout("hypr"; root=fixture.root)
        mujoco = OptionalPackages.checkout("mujoco"; root=fixture.root)
        edg = OptionalPackages.checkout("edg"; root=fixture.root)
        root_alias = joinpath(fixture.temporary, "core-alias")
        symlink(fixture.root, root_alias; dir_target=true)
        before = fixture_snapshot(fixture.root)
        for name in ("hypr", "mujoco", "edg")
            for environment in (fixture.root, fixture.root * "/", joinpath(fixture.root, "."),
                                joinpath(fixture.root, "generated", "environment"), hypr, mujoco, edg,
                                root_alias, joinpath(root_alias, "not-created", "environment"))
                @test_throws Exception OptionalPackages.installation_command(name, environment; root=fixture.root)
            end
            environment = joinpath(fixture.temporary, "core-environments", name * " environment")
            command = OptionalPackages.installation_command(name, environment; root=root_alias)
            @test command isa Cmd
            @test environment in command.exec
            @test "--startup-file=no" in command.exec
            @test "--project=@stdlib" in command.exec
            @test !ispath(environment)
            if name == "hypr"
                @test joinpath(fixture.root, "scripts", "setup_hypr.jl") in command.exec
                @test fixture.root in command.exec
                @test hypr in command.exec
                @test any(argument -> occursin("HYPRInstallation.setup", argument) &&
                            occursin("local_path=ARGS[4]", argument), command.exec)
            else
                @test joinpath(name == "edg" ? edg : mujoco, "scripts", "setup_env.jl") in command.exec
                @test "SPACEAGORA_CORE_PATH=$(fixture.root)" in command.env
                @test "SPACEAGORA_CORE_REV=" in command.env
            end
        end
        @test fixture_snapshot(fixture.root) == before
    end
    # Invalid destinations and missing prerequisites fail before acquisition,
    # installer execution, environment creation or active-project changes.
    optional_fixture() do fixture
        active = Base.active_project()
        before = fixture_snapshot(fixture.root)
        @test_throws Exception OptionalPackages.setup("hypr", fixture.root; root=fixture.root)
        @test_throws Exception OptionalPackages.setup("mujoco", joinpath(fixture.root, "environment"); root=fixture.root)
        @test fixture_snapshot(fixture.root) == before
        @test Base.active_project() == active
        rm(joinpath(fixture.root, "Manifest.toml"))
        before = fixture_snapshot(fixture.root)
        environment = joinpath(fixture.temporary, "external-environment")
        @test_throws Exception OptionalPackages.setup("hypr", environment; root=fixture.root)
        @test !ispath(environment)
        @test fixture_snapshot(fixture.root) == before
        @test Base.active_project() == active
    end
end

@testset "Optional package domains and unpublished EDG are fail closed" begin
    optional_fixture() do fixture
        registry_path = joinpath(fixture.root, "packages", "OptionalPackages.toml")
        original = TOML.parsefile(registry_path)
        for (name, invalid) in (("edg", "packages/robotics/SpaceAGORAEDG.jl"),
                                ("edg", "packages/gnc/../SpaceAGORAEDG.jl"),
                                ("edg", "packages/gnc/nested/SpaceAGORAEDG.jl"),
                                ("mujoco", "packages/gnc/SpaceAGORAMuJoCo.jl"))
            modified = deepcopy(original)
            modified[name]["path"] = invalid
            fixture_toml(registry_path, modified)
            before = fixture_snapshot(fixture.root)
            @test_throws ArgumentError OptionalPackages.spec(name; root=fixture.root)
            @test fixture_snapshot(fixture.root) == before
        end
        fixture_toml(registry_path, original)
        fixture_run("-C", fixture.root, "update-index", "--force-remove", OPTIONAL_FIXTURE_PACKAGES.edg.path)
        before = fixture_snapshot(fixture.root)
        @test_throws ArgumentError OptionalPackages.checkout("edg"; root=fixture.root)
        @test fixture_snapshot(fixture.root) == before
    end
end
