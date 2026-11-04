# Existing EDG numerical suites always require the real provider. Missing
# installation is a test failure; core-only absence is tested separately.
let environment = get(ENV, "SPACEAGORA_EDG_TEST_ENV", "")
    if !isempty(environment)
        isfile(joinpath(environment, "Project.toml")) ||
            error("SPACEAGORA_EDG_TEST_ENV is not a prepared environment.")
        environment in LOAD_PATH || push!(LOAD_PATH, environment)
    end
end
using SpaceAGORA, SpaceAGORAEDG
SpaceAGORA.EDGServices.require_provider()

