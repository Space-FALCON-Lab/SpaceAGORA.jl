# Run the targeted inclination x argument-of-periapsis grid for the shallow
# and deep periapsis regimes, complementing run_inclination_argp_supplement.jl.
#
# This fills out the non-nominal panels in slice_atlas_omega_vs_inclination_reference_apo:
# - shallow and deep periapsis
# - one representative apoapsis per body
# - mass scale 1
# - full_environment with nominal density only
# - full inclination x omega grid
#
# Usage:
#   julia --project=. benchmarks/studies/aerobraking_perturbation_mc/run_inclination_argp_periapsis_supplement.jl --procs 8

include(joinpath(@__DIR__, "study.jl"))

using .AerobrakingPerturbationMC

const SUPPLEMENT_APOAPSIS_ALTITUDES_KM = Dict{Symbol, Vector{Float64}}(
    :mars => [5_000.0],
    :venus => [10_000.0],
    :earth => [36_000.0],
    :titan => [10_000.0],
)

function _parse_supplement_args(args::Vector{String})
    procs = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_PROCS", "0"))
    output_dir = joinpath(AerobrakingPerturbationMC.REPO_ROOT, "output", "aerobraking_perturbation_mc_inclination_argp_periapsis_supplement")
    resume = false
    resume_dir = ""
    worker_max_cases = AerobrakingPerturbationMC.DEFAULT_WORKER_MAX_CASES
    worker_max_rss_gb = AerobrakingPerturbationMC.DEFAULT_WORKER_MAX_RSS_GB
    case_timeout_min = AerobrakingPerturbationMC.DEFAULT_CASE_TIMEOUT_MIN
    aero_solver = AerobrakingPerturbationMC.DEFAULT_AERO_SOLVER
    aero_stiff_tol_scale = AerobrakingPerturbationMC.DEFAULT_AERO_STIFF_TOL_SCALE

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--procs"
            i += 1; procs = parse(Int, args[i])
        elseif arg == "--output-dir"
            i += 1; output_dir = abspath(args[i])
        elseif arg == "--resume"
            resume = true
        elseif arg == "--resume-dir"
            i += 1; resume_dir = abspath(args[i])
        elseif arg == "--worker-max-cases"
            i += 1; worker_max_cases = parse(Int, args[i])
        elseif arg == "--worker-max-rss-gb"
            i += 1; worker_max_rss_gb = parse(Float64, args[i])
        elseif arg == "--case-timeout-min"
            i += 1; case_timeout_min = parse(Float64, args[i])
        elseif arg == "--aero-solver"
            i += 1; aero_solver = Symbol(args[i])
        elseif arg == "--aero-stiff-tol-scale"
            i += 1; aero_stiff_tol_scale = parse(Float64, args[i])
        elseif arg in ("-h", "--help")
            println("""
            Usage:
              julia --project=. benchmarks/studies/aerobraking_perturbation_mc/run_inclination_argp_periapsis_supplement.jl [options]

            Options:
              --procs N
              --output-dir PATH
              --resume
              --resume-dir PATH
              --worker-max-cases N
              --worker-max-rss-gb GB
              --case-timeout-min MINUTES
              --aero-solver MODE
              --aero-stiff-tol-scale X
            """)
            exit(0)
        else
            throw(ArgumentError("Unknown argument '$arg'. Use --help for options."))
        end
        i += 1
    end

    return (;
        procs,
        output_dir,
        resume,
        resume_dir,
        worker_max_cases,
        worker_max_rss_gb,
        case_timeout_min,
        aero_solver,
        aero_stiff_tol_scale,
    )
end

function main(args=ARGS)
    opts = _parse_supplement_args(collect(args))
    AerobrakingPerturbationMC.configure_gram_trajectory_density!()
    AerobrakingPerturbationMC.ensure_aero_perturb_workers!(opts.procs)

    spec = AerobrakingPerturbationMC.StudySpec(
        planets=[:mars, :venus, :earth, :titan],
        periapsis_regimes=[:shallow, :deep],
        apoapsis_altitudes_km=deepcopy(SUPPLEMENT_APOAPSIS_ALTITUDES_KM),
        dynamics_cases=[:full_environment],
        density_scales=Dict(:nominal => 1.0),
        spacecraft_mass_scales=[1.0],
        inclinations_deg=copy(AerobrakingPerturbationMC.DEFAULT_INCLINATIONS_DEG),
        argp_degs=copy(AerobrakingPerturbationMC.DEFAULT_ARGP_DEGS),
        spacecraft_grid=:inclination_argp,
        nominal_mass_scale=1.0,
        nominal_inclination_deg=AerobrakingPerturbationMC.NOMINAL_INCLINATION_DEG,
        nominal_argp_deg=AerobrakingPerturbationMC.NOMINAL_ARGP_DEG,
        norbits=1,
        procs=opts.procs,
        worker_max_cases=opts.worker_max_cases,
        worker_max_rss_gb=opts.worker_max_rss_gb,
        aero_solver=opts.aero_solver,
        aero_stiff_tol_scale=opts.aero_stiff_tol_scale,
        case_timeout_min=opts.case_timeout_min,
        output_dir=opts.output_dir,
        resume=opts.resume,
        resume_dir=opts.resume_dir,
    )
    return AerobrakingPerturbationMC.run_study(spec)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
