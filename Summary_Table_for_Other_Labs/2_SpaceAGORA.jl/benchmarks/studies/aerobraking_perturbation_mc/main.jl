include(joinpath(@__DIR__, "study.jl"))

using Distributed
using .AerobrakingPerturbationMC

function main(args=ARGS)
    spec = AerobrakingPerturbationMC.spec_from_args(collect(args))
    AerobrakingPerturbationMC.configure_gram_trajectory_density!()

    AerobrakingPerturbationMC.ensure_aero_perturb_workers!(spec.procs)

    return AerobrakingPerturbationMC.run_study(spec)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
