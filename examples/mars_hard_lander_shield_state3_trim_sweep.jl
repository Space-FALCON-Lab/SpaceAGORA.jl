include(joinpath(dirname(@__DIR__), "analysis", "mars_hard_lander_prelim", "MarsHardLanderPrelim.jl"))

using .MarsHardLanderPrelim
using CSV
using DataFrames

const OUTPUT_ROOT = "/Users/josephine/Research/Space-FALCON/_Proposals/[FY27]SURP/SHIELD_State3_Trim_Study"

geom = MarsHardLanderPrelim.shield_published_surrogate_geometry()
adapter_cfg = MarsHardLanderPrelim.default_config(
    output_root=OUTPUT_ROOT,
    atmosphere_mode=:gram,
    run_secondary_sweeps=false,
    include_reverse_bangbang_sweep=false,
    enable_crossrange_proxy=false,
    enable_target_range_guidance=false,
    enable_monte_carlo=false,
    generate_plots=false,
    h_switch_grid_m=collect(40e3:-2e3:0.0),
    h_jettison_grid_m=[0.0],
    representative_switches_m=[40e3, 32e3, 24e3, 16e3, 8e3, 0.0],
    geometry=geom,
)

adapter = MarsHardLanderPrelim.build_atmosphere_adapter(adapter_cfg)
rows = DataFrame()

for frac in (0.02, 0.05, 0.10)
    out = joinpath(OUTPUT_ROOT, "frac_" * replace(string(frac), "." => "p"))
    cfg = MarsHardLanderPrelim.default_config(
        output_root=out,
        atmosphere_mode=:gram,
        run_secondary_sweeps=false,
        include_reverse_bangbang_sweep=false,
        enable_crossrange_proxy=false,
        enable_target_range_guidance=false,
        enable_monte_carlo=false,
        generate_plots=false,
        h_switch_grid_m=collect(40e3:-2e3:0.0),
        h_jettison_grid_m=[0.0],
        representative_switches_m=[40e3, 32e3, 24e3, 16e3, 8e3, 0.0],
        geometry=geom,
    )
    result = MarsHardLanderPrelim.run_shield_state3_trim_study(
        cfg,
        adapter;
        base_entry_mass_kg=120.0,
        trim_panel_area_fraction_of_deployed=frac,
        differential_fraction=0.4,
        cant_deg=45.0,
    )
    tmp = result.rows
    tmp[!, :trim_fraction] .= frac
    tmp[!, :shield_areal_density_proxy_kg_m2] .= result.shield_areal_density_proxy_kg_m2
    append!(rows, tmp; cols=:union)
end

mkpath(OUTPUT_ROOT)
CSV.write(joinpath(OUTPUT_ROOT, "shield_state3_trim_sweep.csv"), rows)
