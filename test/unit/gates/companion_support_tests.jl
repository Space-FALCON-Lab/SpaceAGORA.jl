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
        for package in ("SpaceAGORAEDG", "SpaceAGORAMuJoCo")
            repository = realpath(mkpath(joinpath(temporary, package)))
            pin_file = package == "SpaceAGORAEDG" ? "CoreSource.toml" : "SpaceAGORASource.toml"
            marker = package == "SpaceAGORAEDG" ? ".spaceagora-edg-environment.toml" : ".spaceagora-mujoco-environment.toml"
            write(joinpath(repository, "Project.toml"), """
                name = "$package"
                [deps]
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
