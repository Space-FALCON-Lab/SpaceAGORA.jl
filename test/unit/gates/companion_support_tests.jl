using Test, TOML
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "sync_companion_support.jl"))

@testset "Canonical companion support detects drift" begin
    mktempdir() do temporary
        for package in ("HYPR", "SpaceAGORAEDG", "SpaceAGORAMuJoCo", "SpaceAGORA", "SpaceAGORAVisualizer")
            repository = mkpath(joinpath(temporary, package))
            write(joinpath(repository, "Project.toml"), "name = \"$package\"\n")
            expected = CompanionSupportSync.expected_files(repository)
            @test !isempty(expected)
            @test Set(CompanionSupportSync.sync_support(repository)) == Set(keys(expected))
            @test !ispath(joinpath(repository, "scripts")) # checking never edits
            CompanionSupportSync.sync_support(repository; write=true)
            @test isempty(CompanionSupportSync.sync_support(repository))
            file = first(sort!(collect(keys(expected))))
            path = joinpath(repository, file)
            write(path, read(path, String) * "\n# unsynchronized change\n")
            @test file in CompanionSupportSync.sync_support(repository)
            @test_throws ErrorException CompanionSupportSync.main([repository])
            CompanionSupportSync.sync_support(repository; write=true)
            @test isempty(CompanionSupportSync.sync_support(repository))
            rm(path)
            @test file in CompanionSupportSync.sync_support(repository)
            # A generation command must not follow a link into unrelated user data.
            sentinel = joinpath(temporary, "sentinel-$package")
            write(sentinel, "keep")
            symlink(sentinel, path)
            @test_throws ArgumentError CompanionSupportSync.sync_support(repository; write=true)
            @test read(sentinel, String) == "keep"
        end
    end
end

include(joinpath(@__DIR__, "..", "..", "..", "scripts", "companion_support", "installation.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "companion_support", "installation_contract_tests.jl"))
@testset "Canonical installer works before core is installed" begin
    mktempdir() do temporary
        for package in ("HYPR", "SpaceAGORAEDG", "SpaceAGORAMuJoCo")
            repository = realpath(mkpath(joinpath(temporary, package)))
            pin_file = package == "SpaceAGORAEDG" ? "CoreSource.toml" : "SpaceAGORASource.toml"
            marker = "." * lowercase(package) * "-environment.toml"
            core_section = package == "HYPR" ? "weakdeps" : "deps"
            write(joinpath(repository, "Project.toml"), """
                name = "$package"
                [$core_section]
                SpaceAGORA = "afbfb69f-5c0b-4832-b760-43725dff8540"
                """)
            open(joinpath(repository, pin_file), "w") do io
                TOML.print(io, Dict("url" => "https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git",
                    "uuid" => "afbfb69f-5c0b-4832-b760-43725dff8540", "rev" => repeat("a", 40)))
            end
            installer = (;
                PACKAGE_ROOT=repository,
                ENVIRONMENT_MARKER=marker,
                installation_plan=(directory; kwargs...) -> CompanionInstaller.installation_plan(repository, directory;
                    pin_file, marker_name=marker, kwargs...),
                setup=(plan; kwargs...) -> CompanionInstaller.setup(plan;
                    marker_name=marker, package_name=package, kwargs...),
            )
            CompanionInstallationContract.test_installation_contract(installer)
        end
    end
end

@testset "Core dependency identities and optional dependency isolation" begin
    mktempdir() do temporary
        repository = mkpath(joinpath(temporary, "companion"))
        destination = joinpath(temporary, "environment")
        project_path = joinpath(repository, "Project.toml")
        pin_path = joinpath(repository, "SpaceAGORASource.toml")
        core_uuid = "afbfb69f-5c0b-4832-b760-43725dff8540"
        old_loader_uuid = "4096c904-79d0-4d9a-9f5d-f1d6c7a2f1b4"
        pin = Dict("url" => "https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git",
            "uuid" => core_uuid, "rev" => repeat("a", 40))
        write_toml(path, data) = open(io -> TOML.print(io, data), path, "w")
        write_toml(pin_path, pin)
        plan() = CompanionInstaller.installation_plan(repository, destination;
            pin_file="SpaceAGORASource.toml", marker_name=".fixture-environment.toml",
            core_path="", core_revision="")
        project = Dict{String,Any}("name" => "HYPR", "deps" => Dict("Random" => "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"),
            "weakdeps" => Dict("SpaceAGORA" => core_uuid, "SpaceAGORAHYPR" => old_loader_uuid),
            "extras" => Dict("Test" => "8dfed614-e22c-5e08-85e1-65c5234f0b40"))
        write_toml(project_path, project)
        @test Set(keys(plan().dependencies)) == Set(["Random", "Test"])
        @test plan().core.revision == repeat("a", 40)
        @test !ispath(destination)
        # Both tables may declare the same optional extension identity.
        project["deps"]["SpaceAGORA"] = core_uuid
        write_toml(project_path, project)
        @test !haskey(plan().dependencies, "SpaceAGORA")
        @test !haskey(plan().dependencies, "SpaceAGORAHYPR")
        project["deps"]["SpaceAGORA"] = old_loader_uuid
        write_toml(project_path, project)
        @test_throws ArgumentError plan()
        delete!(project["deps"], "SpaceAGORA")
        project["weakdeps"]["SpaceAGORA"] = old_loader_uuid
        write_toml(project_path, project)
        @test_throws ArgumentError plan()
        delete!(project["weakdeps"], "SpaceAGORA")
        write_toml(project_path, project)
        @test_throws ArgumentError plan()
        project["weakdeps"]["SpaceAGORA"] = core_uuid
        write_toml(project_path, project)
        pin["uuid"] = old_loader_uuid
        write_toml(pin_path, pin)
        @test_throws ArgumentError plan()
        @test !ispath(destination)
    end
end
