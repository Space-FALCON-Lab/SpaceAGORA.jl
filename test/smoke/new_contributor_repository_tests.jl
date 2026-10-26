# Self-test of scripts/new_contributor_repository.jl, with the standard library only.
# Core CI runs it in the contributor-kit job, before that job tests each generated
# starter through the contributor workflow.
using Dates, Test, TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "scripts", "new_contributor_repository.jl"))
include(joinpath(ROOT, "scripts", "contributor_environment.jl"))
const G = NewContributorRepository
const STARTERS = joinpath(ROOT, "templates", "starters")

git(dir, args...) = run(`git -c maintenance.auto=false -c gc.auto=0 -c user.name=t -c user.email=t@example.com -C $dir $args`)
refused(f, fragment) = try f(); false catch e; e isa ArgumentError && occursin(fragment, e.msg) end
files_under(dir) = sort!([relpath(joinpath(folder, file), dir) for (folder, _, files) in walkdir(dir) for file in files])
clean_checkout() = isempty(readchomp(`git -C $ROOT status --porcelain=v1 --untracked-files=no --ignore-submodules=all`))

@testset "new_contributor_repository.jl" begin
    @testset "every template file parses with its placeholders" begin
        for kind in G.KINDS
            dir = joinpath(STARTERS, kind)
            @test isfile(joinpath(dir, "src", "StarterName.jl"))
            for rel in files_under(dir)
                text = read(joinpath(dir, rel), String)
                if endswith(rel, ".jl")
                    @test Meta.parseall(text; filename=rel) isa Expr
                elseif endswith(rel, ".toml")
                    @test TOML.tryparse(text) isa Dict
                end
            end
            pin = TOML.parsefile(joinpath(dir, "SpaceAGORASource.toml"))
            @test pin["rev"] == G.REV_TOKEN && pin["version"] == "0.0.0"
            # A package or a companion ships core's MIT license with the year and holder left open.
            if kind == "project"
                @test !isfile(joinpath(dir, "LICENSE"))
            else
                core_license = readlines(joinpath(ROOT, "LICENSE"))
                license = readlines(joinpath(dir, "LICENSE"))
                @test license[3] == G.LICENSE_TOKEN
                @test license[[1:2; 4:end]] == core_license[[1:2; 4:end]]
            end
            @test !haskey(get(TOML.parsefile(joinpath(dir, "Project.toml")), "sources", Dict()), "SpaceAGORA")
        end
    end

    # The kit's clone commands keep SpaceAGORA's LFS data out of the clone for good: a
    # variable on the clone alone would not cover the checkout of the pinned commit.
    # They skip old file versions too, and the double quotes reach Git as `*` in bash,
    # PowerShell and Command Prompt alike.
    @testset "clone commands exclude SpaceAGORA's LFS data and old file versions" begin
        pages = [joinpath(ROOT, "docs", "src", "user", "research_repository.md");
                 [joinpath(STARTERS, kind, "README.md") for kind in G.KINDS]]
        for page in pages
            clones = filter(line -> occursin("git clone", line) && occursin("SpaceAGORA.jl.git", line), readlines(page))
            @test !isempty(clones)
            @test all(line -> occursin("git clone --filter=blob:none -c lfs.fetchexclude=\"*\" ", line), clones)
            @test !occursin("GIT_LFS_SKIP_SMUDGE", read(page, String))
        end
    end

    @testset "names, series and targets" begin
        @test G.check_name("MyPlanner") == "MyPlanner"
        for bad in ("myPlanner", "My-Planner", "My Planner", "SpaceAGORA", "StarterName", "MyPlanner.jl", "")
            @test refused(() -> G.check_name(bad), "--name")
        end
        # A name that contains the placeholder is refused before anything is written.
        @test refused(() -> G.check_name("MyStarterName"), "the starters' placeholder")
        mktempdir() do parent
            target = joinpath(parent, "MyStarterName")
            @test refused(() -> G.generate("package", target; name="MyStarterName", require_clean=false),
                "the starters' placeholder")
            @test !ispath(target)
        end
        @test G.compat_series(v"0.2.0") == "0.2" && G.compat_series(v"0.3.4") == "0.3" && G.compat_series(v"1.4.0") == "1"
        mktempdir() do dir
            write(joinpath(dir, "keep.txt"), "")
            @test refused(() -> G.check_target(dir, ROOT), "is not empty")
        end
        @test refused(() -> G.check_target(joinpath(ROOT, "starter-probe"), ROOT), "inside the SpaceAGORA checkout")
        @test refused(() -> G.generate("library", mktempdir(); name="Abc", require_clean=false), "Unknown kind")
        @test G.check_holder("  Space-FALCON Lab ") == "Space-FALCON Lab"
        @test refused(() -> G.check_holder(" "), "--holder must name")
        @test refused(() -> G.check_holder("A\nB"), "single line")
        mktempdir() do parent
            target = joinpath(parent, "Abc")
            @test refused(() -> G.generate("package", target; name="Abc", holder="", require_clean=false), "--holder")
            @test !ispath(target)
        end
    end

    @testset "a checkout with tracked changes is refused" begin
        mktempdir() do root
            write(joinpath(root, "Project.toml"), "name = \"SpaceAGORA\"\nuuid = \"afbfb69f-5c0b-4832-b760-43725dff8540\"\nversion = \"0.2.0\"\n")
            facade = joinpath(root, G.FACADE_SOURCE)
            mkpath(dirname(facade))
            write(facade, "const CONTRACT_VERSION = v\"2.0.0\"\n")
            git(root, "init", "-q")
            git(root, "add", ".")
            git(root, "commit", "-q", "-m", "fixture")
            core = G.core_identity(root)
            @test occursin(r"^[0-9a-f]{40}$", core.rev) && core.version == v"0.2.0" && core.contract == v"2.0.0"
            write(facade, "const CONTRACT_VERSION = v\"2.1.0\"\n")
            @test refused(() -> G.core_identity(root), "has tracked changes")
            @test G.core_identity(root; require_clean=false).contract == v"2.1.0"
        end
    end

    @testset "the starter comes from the commit, not from untracked files" begin
        mktempdir() do root
            write(joinpath(root, "Project.toml"), "name = \"SpaceAGORA\"\nuuid = \"afbfb69f-5c0b-4832-b760-43725dff8540\"\nversion = \"0.2.0\"\n")
            facade = joinpath(root, G.FACADE_SOURCE)
            mkpath(dirname(facade))
            write(facade, "const CONTRACT_VERSION = v\"2.0.0\"\n")
            template = joinpath(root, "templates", "starters", "package")
            mkpath(dirname(template))
            cp(joinpath(STARTERS, "package"), template)
            git(root, "init", "-q")
            git(root, "add", ".")
            git(root, "commit", "-q", "-m", "fixture")
            folder = joinpath("templates", "starters", "package")
            committed = sort!([relpath(path, folder) for path in
                split(readchomp(`git -C $root ls-tree -r --name-only HEAD -- $folder`), '\n')])
            write(joinpath(template, "untracked.jl"), "error(\"not part of the commit\")\n")
            write(joinpath(template, "src", "untracked_helper.jl"), "x = 1\n")
            @test isempty(readchomp(`git -C $root status --porcelain=v1 --untracked-files=no`))
            mktempdir() do parent
                target = G.generate("package", joinpath(parent, "LinePlanner"); name="LinePlanner", root)
                @test sort!(map(path -> replace(path, "LinePlanner" => "StarterName"), files_under(target))) == committed
                @test !isfile(joinpath(target, "untracked.jl")) && !isfile(joinpath(target, "src", "untracked_helper.jl"))
            end
            # Work on the templates themselves uses the working tree.
            @test "untracked.jl" in first.(G.template_files(root, "package"; committed=false))
        end
    end

    @testset "each starter, filled in for this checkout" begin
        core = G.core_identity(ROOT; require_clean=false)
        uuids = String[]
        mktempdir() do parent
            for (kind, name) in (("project", "OrbitStudy"), ("package", "LinePlanner"), ("companion", "CircleCompanion"))
                target = G.generate(kind, joinpath(parent, name); name, require_clean=false)
                rel = files_under(target)
                @test "src/$(name).jl" in rel && !any(contains("StarterName"), rel)
                @test ("Manifest.toml" in rel) == (kind == "project")
                for path in rel
                    text = read(joinpath(target, path), String)
                    @test !occursin(G.NAME_TOKEN, text) && !occursin(G.UUID_TOKEN, text) && !occursin(G.REV_TOKEN, text)
                    endswith(path, ".jl") && @test Meta.parseall(text; filename=path) isa Expr
                end
                project = TOML.parsefile(joinpath(target, "Project.toml"))
                @test project["name"] == name
                push!(uuids, project["uuid"])
                @test project["deps"]["SpaceAGORA"] == core.uuid
                @test project["compat"]["SpaceAGORA"] == G.compat_series(core.version)
                pin = TOML.parsefile(joinpath(target, "SpaceAGORASource.toml"))
                @test pin["rev"] == core.rev && pin["uuid"] == core.uuid && pin["version"] == string(core.version)
                workflow = read(joinpath(target, ".github", "workflows", "ci.yml"), String)
                @test occursin("contributor-ci.yml@$(core.rev)", workflow)
                @test occursin("spaceagora-rev: \"$(core.rev)\"", workflow) && occursin("kind: $kind", workflow)
                if kind == "companion"
                    @test occursin("const CORE_CONTRACT_VERSION = v\"$(core.contract)\"", read(joinpath(target, "src", "$(name).jl"), String))
                end
                if kind != "project"
                    @test readlines(joinpath(target, "LICENSE"))[3] ==
                          "Copyright (c) $(Dates.year(Dates.today())) $(G.HOLDER_PLACEHOLDER)"
                end
                # The contributor workflow's pin check accepts it, at this commit when the checkout is clean.
                if clean_checkout()
                    checked = ContributorEnvironment.check_repository(kind, target; rev=core.rev)
                    @test checked.pin == "rev $(core.rev)"
                end
                # Projects and packages use root names and SpaceAGORA.Conformance only.
                kind == "companion" || @test isempty(ContributorEnvironment.scan(kind, target))
            end
        end
        @test allunique(uuids)
        # The holder and the year fill the license line.
        mktempdir() do parent
            target = G.generate("package", joinpath(parent, "HeldPlanner"); name="HeldPlanner", holder="Space-FALCON Lab",
                                year=2026, require_clean=false)
            @test readlines(joinpath(target, "LICENSE"))[3] == "Copyright (c) 2026 Space-FALCON Lab"
        end
    end
end
