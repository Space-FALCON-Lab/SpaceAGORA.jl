include(joinpath(@__DIR__, "study.jl"))

using Distributed
using .AerobrakingPerturbationMC

function main(args=ARGS)
    spec = AerobrakingPerturbationMC.spec_from_args(collect(args))
    AerobrakingPerturbationMC.configure_gram_trajectory_density!()

    if spec.procs > 0 && nworkers() < spec.procs
        study_dir = @__DIR__
        project = AerobrakingPerturbationMC.REPO_ROOT
        println("[aero-perturb] starting $(spec.procs - nworkers()) worker process(es); this can take a bit on first run")
        flush(stdout)
        addprocs(spec.procs - nworkers(); exeflags="--project=$(project) --compiled-modules=existing")
        for worker in workers()
            println("[aero-perturb] initializing worker=$worker")
            flush(stdout)
            ex = quote
                import Pkg
                Pkg.activate($project; io=devnull)
                include(joinpath($study_dir, "study.jl"))
                nothing
            end
            remotecall_wait(Core.eval, worker, Main, ex)
            println("[aero-perturb] worker=$worker ready")
            flush(stdout)
        end
    end

    return AerobrakingPerturbationMC.run_study(spec)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
