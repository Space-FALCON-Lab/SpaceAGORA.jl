include("CampaignPaths.jl")
include(joinpath(CampaignPaths.BASE_CODE, "StudyRunner.jl"))

function run_worker(args; allowed_policies=(:naive_next_entering, :positive_along_track))
    5 <= length(args) <= 8 || error(
        "Usage: run_case.jl OUTPUT_ROOT TARGET_KM INCLINATION_DEG HELPERS POLICY [ORBITS] [POWER_W] [HELPER_KM]")
    root, altitude, inclination, helpers, policy_name = args[1:5]
    policy = Symbol(policy_name)
    policy in allowed_policies || error("Unsupported policy for this output: $policy")
    parse(Float64, altitude) in (1000, 1050, 1100) ||
        error("This comparison requires target altitude 1000, 1050, or 1100 km")
    orbits = length(args) >= 6 ? parse(Float64, args[6]) : 100.0
    power = length(args) >= 7 ? parse(Float64, args[7]) :
        (policy === :no_laser ? 0.0 : 10_000.0)
    helper_altitude = length(args) >= 8 ? parse(Float64, args[8]) : 1000.0
    case = OracleStudy.StudyCase(parse(Float64, altitude), parse(Float64, inclination),
        parse(Int, helpers), policy; orbits, power_W=power, helper_altitude_km=helper_altitude)
    OracleStudy.run_case(case; output_root=abspath(root),
        provenance=CampaignPaths.campaign_sources(OracleStudy.StudyConfig.source_hashes))
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_worker(ARGS)
end
