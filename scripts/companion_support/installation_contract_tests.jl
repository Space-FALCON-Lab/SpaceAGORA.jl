# Common behavioral contract for standalone installers. The assertions are
# independent of the path-validation implementation; no Pkg APIs are invoked.
module CompanionInstallationContract
using Test, TOML

function test_installation_contract(installation)
    @testset "Shared installation destination contract" begin
        mktempdir() do temporary
            temporary = realpath(temporary)
            core = mkpath(joinpath(temporary, "core"))
            write(joinpath(core, "Project.toml"), "name = \"SpaceAGORA\"\nuuid = \"afbfb69f-5c0b-4832-b760-43725dff8540\"\n")
            env = joinpath(temporary, "environment")
            plan(path=env) = installation.installation_plan(path; core_path=core, core_revision="")
            active = Base.active_project()
            @test plan().core.path == core
            @test !ispath(env)
            @test_throws ArgumentError plan(joinpath(core, "new", "environment"))
            @test_throws ArgumentError plan(joinpath(installation.PACKAGE_ROOT, "new", "environment"))
            alias = joinpath(temporary, "core-alias")
            symlink(core, alias; dir_target=true)
            @test_throws ArgumentError plan(joinpath(alias, "new", "environment"))
            dangling = joinpath(temporary, "dangling")
            symlink(joinpath(temporary, "absent"), dangling; dir_target=true)
            @test_throws ArgumentError plan(joinpath(dangling, "environment"))
            file = joinpath(temporary, "file")
            write(file, "keep")
            @test_throws ArgumentError plan(file)
            @test_throws ArgumentError plan(joinpath(file, "environment"))
            @test read(file, String) == "keep"

            prepared = plan()
            mkpath(env)
            sentinel = joinpath(env, "sentinel")
            write(sentinel, "user data")
            @test_throws ArgumentError plan()
            # Revalidate at setup time, before acquisition or package mutation.
            @test_throws ArgumentError installation.setup(prepared; precompile=false)
            @test read(sentinel, String) == "user data"
            @test !isfile(joinpath(env, "Project.toml"))
            marker = joinpath(env, installation.ENVIRONMENT_MARKER)
            function record(companion, source)
                open(marker, "w") do io
                    TOML.print(io, Dict("companion" => companion, "core" => source))
                end
            end
            record(installation.PACKAGE_ROOT, core)
            @test plan().destination == env
            @test read(sentinel, String) == "user data"
            record(installation.PACKAGE_ROOT, "another core")
            @test_throws ArgumentError plan()
            record(temporary, core)
            @test_throws ArgumentError plan()
            record(installation.PACKAGE_ROOT, core)
            original_marker = read(marker)
            external_marker = joinpath(temporary, "external-marker")
            mv(marker, external_marker)
            symlink(external_marker, marker)
            @test_throws ArgumentError plan()
            @test_throws ArgumentError installation.setup(prepared; precompile=false)
            @test read(external_marker) == original_marker
            @test Base.active_project() == active
            @test !ispath(joinpath(core, "new"))
            @test !ispath(joinpath(temporary, "absent"))
        end
    end
end
end
