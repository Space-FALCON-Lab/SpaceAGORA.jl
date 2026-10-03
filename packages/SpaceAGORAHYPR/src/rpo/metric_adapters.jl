# HYPR keeps the existing signatures; shared kernels consume numeric inputs only.
"""Forward existing HYPR cost-normalization inputs to the shared metric kernel."""
function rpo_path_cost_normalization_refs(points, cfg::RPOPSOConfig)
    return rpo_path_cost_normalization_refs(points;
        cost_ref_distance_m=cfg.cost_ref_distance_m, sample_ds_m=cfg.sample_ds_m,
        tf_s=cfg.tf_s, mass_kg=cfg.mass_kg, isp_s=cfg.isp_s, g0_mps2=cfg.g0_mps2)
end

"""Forward existing HYPR finite-difference fuel inputs to the shared metric kernel."""
function rpo_fuel_proxy_from_samples(samples, cfg::RPOPSOConfig)
    return rpo_fuel_proxy_from_samples(samples;
        tf_s=cfg.tf_s, mass_kg=cfg.mass_kg, isp_s=cfg.isp_s, g0_mps2=cfg.g0_mps2)
end
