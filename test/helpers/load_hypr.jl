# Full-feature tests use an isolated installation. The default resolves the
# reviewed Git pin; an explicit SPACEAGORA_HYPR_PATH is development-only.
if !isdefined(@__MODULE__, :SpaceAGORAHYPR)
    include(joinpath(@__DIR__, "..", "..", "scripts", "setup_hypr.jl"))
    environment = get(ENV, "SPACEAGORA_HYPR_TEST_ENV", "")
    if isempty(environment)
        environment = mktempdir(; prefix="spaceagora-hypr-tests-")
        HYPRInstallation.setup(environment)
    end
    environment in LOAD_PATH || push!(LOAD_PATH, environment)
end
using SpaceAGORAHYPR
