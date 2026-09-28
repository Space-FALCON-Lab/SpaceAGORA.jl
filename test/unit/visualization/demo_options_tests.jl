module ViewerDemoOptionsTests
using Test
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "viewer_demos", "demo_options.jl"))
@testset "viewer demo options preserve existing outputs" begin
    mktempdir() do dir
        env=Dict("SPACEAGORA_VIEWER_DEMO_OUT"=>dir)
        opts=viewer_demo_options("sample",60;argv=String[],env)
        @test opts==(duration_s=60.0,output_dir=joinpath(dir,"sample"))
        @test !ispath(opts.output_dir)
        explicit=joinpath(dir,"explicit")
        @test viewer_demo_options("sample",60;argv=["--duration-s","2.5","--output-dir",explicit],env).duration_s==2.5
        @test viewer_demo_options("sample",60;argv=["--output-dir",explicit],env).output_dir==explicit
        for argv in (["--duration-s","NaN"],["--duration-s","Inf"],["--duration-s","0"],["--duration-s","-1"],
            ["--duration-s","junk"],["--duration-s"],["--unknown","1"],["--duration-s","1","--duration-s","2"],
            ["--output-dir",""],["--output-dir"])
            @test_throws ArgumentError viewer_demo_options("sample",60;argv,env)
        end
        mkpath(explicit)
        @test viewer_demo_options("sample",60;argv=["--output-dir",explicit],env).output_dir==explicit
        write(joinpath(explicit,"keep.txt"),"previous results")
        @test_throws ArgumentError viewer_demo_options("sample",60;argv=["--output-dir",explicit],env)
        @test read(joinpath(explicit,"keep.txt"),String)=="previous results"
        @test_throws ArgumentError viewer_demo_options("sample",60;argv=["--output-dir",joinpath(explicit,"keep.txt")],env)
        @test_throws ArgumentError viewer_demo_options("sample",Inf;argv=String[],env)
    end
end
end
