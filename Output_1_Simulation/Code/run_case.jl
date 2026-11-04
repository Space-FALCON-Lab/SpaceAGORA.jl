include("StudyRunner.jl")
using .OracleStudy

function main(args)
    5 <= length(args) <= 8 || error(
        "Usage: julia run_case.jl OUTPUT_ROOT TARGET_KM INCLINATION_DEG HELPERS POLICY [ORBITS=100] [POWER_W=10000] [HELPER_KM=1000]")
    root, altitude, inclination, helpers, policy = args[1:5]
    orbits = length(args) >= 6 ? parse(Float64, args[6]) : 100.0
    power = length(args) >= 7 ? parse(Float64, args[7]) : 10_000.0
    helper_altitude = length(args) >= 8 ? parse(Float64, args[8]) : 1000.0
    case = StudyCase(parse(Float64, altitude), parse(Float64, inclination),
        parse(Int, helpers), Symbol(policy); orbits, power_W=power,
        helper_altitude_km=helper_altitude)
    run_case(case; output_root=abspath(root))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
