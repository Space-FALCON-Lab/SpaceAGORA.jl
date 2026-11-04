include(joinpath(@__DIR__, "..", "Simulation", "Code", "StudyConfig.jl"))
using .StudyConfig
using TOML

const EVIDENCE = joinpath(@__DIR__, "results", "production_pipeline_compact")
const TEST_SCRIPT = joinpath(@__DIR__, "production_pipeline_tests.jl")
const VERIFY_SCRIPT = joinpath(@__DIR__, "..", "Simulation", "Code", "verify_output.jl")

function child(command, label)
    open(joinpath(EVIDENCE, label * ".log"), "w") do io
        run(pipeline(command; stdout=io, stderr=io))
    end
end

function main()
    mkpath(EVIDENCE)
    julia = Base.julia_cmd()
    outputs = Dict{String,String}()
    for (label, h, inc, helpers, policy, orbits, mode) in (
        ("final_sma_recorded", 850, 1, 50, :gve_sma, 2, "recorded"),
        ("final_sma_control", 850, 1, 50, :gve_sma, 2, "control"),
        ("final_inc_recorded", 1050, 1, 50, :gve_inc, 2, "recorded"),
        ("final_inc_control", 1050, 1, 50, :gve_inc, 2, "control"),
        ("final_on", 1050, 1, 1, :gve_sma, 0.02, "recorded"),
        ("final_zero", 1050, 1, 1, :gve_sma, 0.02, "zero"),
        ("final_coplanar_inc", 1050, 0, 1, :gve_inc, 0.02, "recorded"))
        root = joinpath(EVIDENCE, label)
        child(`$julia --startup-file=no --threads=1 $TEST_SCRIPT pilot $root $h $inc $helpers $policy $orbits $mode`, label)
        case = StudyCase(h, inc, helpers, policy; orbits, power_W=mode == "zero" ? 0.0 : 10_000.0)
        directory = case_directory(root, case)
        outputs[label] = directory
        if mode != "control"
            child(`$julia --startup-file=no --threads=1 $VERIFY_SCRIPT $directory`, label * "_verify")
        end
    end
    for name in ("sma", "inc")
        recorded, control = outputs["final_$(name)_recorded"], outputs["final_$(name)_control"]
        child(`$julia --startup-file=no --threads=1 $TEST_SCRIPT compare $recorded $control`,
            "final_$(name)_noninterference")
    end
    atomic_toml(joinpath(EVIDENCE, "short_pilot_outputs.toml"), outputs)
    println("All short production pilots and recorder-control comparisons passed.")
end

main()
