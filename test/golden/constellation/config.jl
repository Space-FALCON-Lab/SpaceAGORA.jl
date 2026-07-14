function build_constellation_config()
    sc_a = make_spacecraft(ra_alt_m=520e3, rp_alt_m=430e3, ν_deg=150.0)
    sc_b = make_spacecraft(ra_alt_m=700e3, rp_alt_m=650e3, ν_deg=110.0)
    sc_c = make_spacecraft(ra_alt_m=560e3, rp_alt_m=560e3, i_deg=53.0, ν_deg=60.0)

    return build_config_multi(
        spacecraft=[sc_a, sc_b, sc_c],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=3600.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
end
