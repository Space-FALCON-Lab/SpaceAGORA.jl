function _spaceagora_precompile_args()
    SM = SimulationModel
    planet = SM.Mars()
    spacecraft = TelemetryVerification.make_three_body_spacecraft(
        bus_dims=(1.2, 1.1, 0.9),
        panel_dims=(0.01, 0.8, 0.4),
        bus_mass=150.0,
        panel_mass_each=2.0,
        panel_offset_y=0.7,
        ic=SM.InitialCondition(
            ra=planet.Rp_e + 220e3,
            rp=planet.Rp_e + 150e3,
            i=28.0,
            ω=10.0,
            Ω=15.0,
            ν=165.0
        ),
        prop_mass=15.0,
        id=1
    )

    return TelemetryVerification.make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=5.0,
        initial_time=SM.InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredGravityModel(),),
        density_model=SM.ExponentialAtmosphereModel(planet),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=false,
        keplerian=true,
        EI_km=140.0,
        verbose=false,
        results=false,
        results_directory=joinpath(tempdir(), "spaceagora_precompile")
    )
end

const _SPACEAGORA_PRECOMPILE_ENV = Dict("SPACEAGORA_PARALLEL_PROFILE" => "R2", "SPACEAGORA_SAVE_BUNDLE" => "0", "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0")

function _run_spaceagora_precompile_workload(; workspace::AbstractString=tempdir())
    parse_parallel_profile("R2")
    engine_config = simulation_engine_config_from_env(_SPACEAGORA_PRECOMPILE_ENV)
    args = _spaceagora_precompile_args()
    mktempdir(workspace) do tmp
        cd(tmp) do
            run_simulation(engine_config, args; return_solution=true)
        end
    end
end

@setup_workload begin
    @compile_workload _run_spaceagora_precompile_workload()
    # MANDATORY whenever a workload above touches a SPICE-backed planet, and
    # cheap insurance when none does. `_FURNISHED_KERNELS` and the planet
    # instance caches are module-level `const`s, so anything furnished HERE is
    # serialised into the pkgimage; at run time `_furnsh_once` would then skip
    # the furnish for a kernel CSPICE never actually loaded, and every lookup
    # needing it fails against an empty pool (`utc2et` losing the leapseconds
    # kernel is the first symptom). Same hazard `_reset_furnished_kernels!`
    # documents for `kclear()`, reached by a different route. Verified by
    # observation, not theory: a workload constructing `Earth(...)` here left
    # every subsequent process unable to resolve a UTC epoch until the pkgimage
    # was rebuilt.
    SimulationModel.Planets._reset_furnished_kernels!()
end
