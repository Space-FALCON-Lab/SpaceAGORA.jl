# Adds the heat-limited altitude columns (h_qhigh_m, h_qlow_m) to the previous study's 50-deg passive passes
# (no simulation; fit evaluation only). Usage: julia --project=. run_static_refband.jl <ref_dir> <out_dir>
include(joinpath(@__DIR__, "common.jl"))
function refband_main(ref, out)
    study = load_study(); cfg = make_study_config(study; raan_deg=0.0, omega_deg=0.0); fit = fitted_atmosphere_parameters(cfg)
    for tag in readdir(joinpath(ref, "grid"))
        df = CSV.read(joinpath(ref, "grid", tag, "passes_3b.csv"), DataFrame)
        b = [band_altitudes(fit, v, study) for v in df.vrel_mps]
        df.h_qhigh_m = first.(b); df.h_qlow_m = last.(b)
        mkpath(joinpath(out, "ref50", tag)); CSV.write(joinpath(out, "ref50", tag, "passes_3b.csv"), df)
    end
end
refband_main(abspath(ARGS[1]), abspath(ARGS[2]))
