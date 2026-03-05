const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "scripts", "performance_paper_pipeline.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
