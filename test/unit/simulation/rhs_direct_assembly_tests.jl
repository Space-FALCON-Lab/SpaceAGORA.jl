using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel
using ComponentArrays

const RHSDA_SE = SpaceAGORA.SimulationEngine

function _rhs_direct_assembly_config(; n_sats::Int=3)
    planet = Earth()
    spacecraft = SpacecraftModel[]
    for i in 1:n_sats
        root = Link(root=true, m=500.0, ref_area=12.0)
        ic = InitialCondition(
            ra=planet.Rp_e + 550_000.0 + 100.0 * i,
            rp=planet.Rp_e + 550_000.0 + 100.0 * i,
            i=53.0,
            ω=0.0,
            Ω=10.0,
            ν=360.0 * (i - 1) / n_sats,
        )
        push!(
            spacecraft,
            SpacecraftModel(
                Joint[],
                [root],
                root,
                true,
                500.0,
                0.0,
                root.inertia,
                0,
                0,
                ic,
                i,
            ),
        )
    end
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false,
            save_csv=false,
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=120.0,
            orientation_sim=false,
            num_steps_to_save=20,
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=300.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false,
            ephemerides_model=SimpleEphemeridesModel(),
        ),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-9,
            abstol_orbit=1e-9,
            dt_max_orbit=2.0,
        ),
    )
end

@testset "RHS direct final assembly writes translational derivatives" begin
    args = _rhs_direct_assembly_config(n_sats=3)
    u = RHSDA_SE.build_initial_conditions(args)
    du = zero(u)
    p = ODEParams(n_sats=3, args=args)
    p.is_active[2] = false

    totals = zeros(Float64, 6, 3)
    totals[1:3, 1] .= (6.0, -9.0, 12.0)
    totals[1:3, 2] .= (10.0, 20.0, 30.0)
    totals[1:3, 3] .= (-15.0, 18.0, -21.0)

    env = withenv("SPACEAGORA_RHS_FINAL_ASSEMBLY_DIRECT_LAYOUT" => "1") do
        RHSDA_SE._snapshot_rhs_plan_env_config()
    end
    @test env.final_assembly_direct_layout
    @test RHSDA_SE._try_assign_flat_translational_rhs_direct_layout!(
        du,
        u,
        p,
        totals,
        env,
        :full,
        2,
    )

    for sat_idx in (1, 3)
        @test collect(du.sc[sat_idx].pos) == collect(u.sc[sat_idx].vel)
        @test collect(du.sc[sat_idx].vel) ≈ collect(totals[1:3, sat_idx]) ./ u.sc[sat_idx].mass
        @test du.sc[sat_idx].mass == 0.0
        @test all(iszero, du.sc[sat_idx].heat_loads)
    end
    @test all(iszero, ComponentArrays.getdata(du.sc[2]))
end

@testset "RHS direct final assembly declines unsupported cases" begin
    args = _rhs_direct_assembly_config(n_sats=2)
    u = RHSDA_SE.build_initial_conditions(args)
    du = zero(u)
    p = ODEParams(n_sats=2, args=args)
    totals = zeros(Float64, 6, 2)

    disabled = withenv("SPACEAGORA_RHS_FINAL_ASSEMBLY_DIRECT_LAYOUT" => "0") do
        RHSDA_SE._snapshot_rhs_plan_env_config()
    end
    @test !RHSDA_SE._try_assign_flat_translational_rhs_direct_layout!(
        du,
        u,
        p,
        totals,
        disabled,
        :full,
        1,
    )

    enabled = withenv("SPACEAGORA_RHS_FINAL_ASSEMBLY_DIRECT_LAYOUT" => "1") do
        RHSDA_SE._snapshot_rhs_plan_env_config()
    end
    @test !RHSDA_SE._try_assign_flat_translational_rhs_direct_layout!(
        du,
        u,
        p,
        totals,
        enabled,
        :implicit,
        1,
    )
end
