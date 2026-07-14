@testset "Run Simulation Metadata Return" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        keplerian=true
    )

    metadata = mktempdir() do tmp
        cd(tmp) do
            run_simulation(args; return_solution=true, return_solver_metadata=true)
        end
    end

    @test metadata isa NamedTuple
    @test hasproperty(metadata, :solution)
    @test hasproperty(metadata, :solver_mode)
    @test hasproperty(metadata, :solver_trace)
    @test hasproperty(metadata, :parallel_policy)
    @test hasproperty(metadata, :spice_counters)
end

@testset "Run Simulation Metadata Return (Gravity Backbone)" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        keplerian=true
    )

    metadata = mktempdir() do tmp
        cd(tmp) do
            withenv(
                "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
                "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "2.0"
            ) do
                run_simulation(args; return_solution=true, return_solver_metadata=true)
            end
        end
    end

    @test metadata isa NamedTuple
    @test metadata.solver_mode == "gravity_backbone_split"
    @test _is_gravity_backbone_state(metadata.solution.u[end])
    @test isapprox(
        _state_mass_kg(metadata.solution.u[end], args, 1),
        args.dynamics_model.spacecraft[1].dry_mass + args.dynamics_model.spacecraft[1].prop_mass;
        atol=0.0,
        rtol=0.0
    )
    @test isempty(_state_heat_loads(metadata.solution.u[end], args, 1)) == false
end

@testset "RHS Completeness: Mass Derivative" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )

    u = build_initial_conditions(args)
    du = copy(u)
    du.sc[1].mass = 789.0
    p = ODEParams(n_sats=1, args=args)
    spacecraft_dynamics!(du, u, p, 0.0)
    @test du.sc[1].mass == 0.0

    du_inactive = copy(u)
    du_inactive.sc[1].mass = 123.0
    p_inactive = ODEParams(n_sats=1, args=args, is_active=[false])
    spacecraft_dynamics!(du_inactive, u, p_inactive, 0.0)
    @test du_inactive.sc[1].mass == 0.0
end







