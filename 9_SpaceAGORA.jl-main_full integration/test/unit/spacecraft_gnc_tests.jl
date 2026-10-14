using SpaceAGORA, Test, SpaceAGORA.SimulationModel
using StaticArrays

# Cheap indexed control effector: a constant inertial force on one spacecraft.
struct SgPushEffector <: SpaceAGORA.AbstractControlEffectorModel
    sat_idx::Int
    force_n::SVector{3, Float64}
end
SpaceAGORA.calcControlEffect!(::SgPushEffector, u, p::SpaceAGORA.SimulationModel.ODEParams, t::Float64, i::Int) = nothing
SpaceAGORA.calcControlForceTorque(m::SgPushEffector, u::AbstractVector, p::SpaceAGORA.SimulationModel.ODEParams, i::Int, t::Float64) =
    i == m.sat_idx ? (m.force_n, SVector{3, Float64}(0.0, 0.0, 0.0)) :
                     (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
SpaceAGORA.bind_spacecraft(m::SgPushEffector, sat_idx::Int) = SgPushEffector(sat_idx, m.force_n)

const SG_PLANET = SimulationModel.make_no_gram_planet(:earth)
const SG_PUSH = SVector{3, Float64}(0.0, 5.0, 0.0)

function sg_sat(id::Int, raan_deg::Float64; guidance=nothing, navigation=nothing, control=nothing)
    root = Link(root=true, m=100.0, ref_area=2.0)
    ic = InitialCondition(ra=SG_PLANET.Rp_e + 700e3, rp=SG_PLANET.Rp_e + 650e3, i=45.0, ω=0.0, Ω=raan_deg, ν=10.0)
    kw = control === nothing ? (;) : (; control)
    guidance === nothing || (kw = merge(kw, (; guidance)))
    navigation === nothing || (kw = merge(kw, (; navigation)))
    return SpacecraftModel(joints=Joint[], links=[root], root=root, prop_mass=0.0,
        inertia_tensor=root.inertia, initial_condition=ic, id=id; kw...)
end

sg_control(effectors...; rate=1.0) =
    ControlModel(control_effectors=Tuple(effectors), control_rates=fill(rate, length(effectors)))

function sg_config(sats; control=sg_control())
    planet = SG_PLANET
    env = EnvironmentModel(planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=300.0, orientation_sim=false, num_steps_to_save=1000),
        environment_model=env,
        dynamics_model=DynamicsModel(SpacecraftModel[sats...], (InverseSquaredJ2GravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=control,
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

final_state(args) = collect(run_simulation(args; return_solution=true).u[end])

@testset "Per-spacecraft GNC" begin
    @testset "11-argument positional constructor" begin
        root = Link(root=true, m=100.0, ref_area=2.0)
        sc = SpacecraftModel(Joint[], [root], root, true, root.m, 0.0, root.inertia, 0, 0, InitialCondition(), 7)
        @test sc.id == 7
        @test isempty(sc.guidance.guidance_effectors) && isempty(sc.navigation.navigation_effectors) &&
              isempty(sc.control.control_effectors)
    end

    @testset "equivalence with configuration-level declaration" begin
        push_cfg = sg_config([sg_sat(1, 0.0)]; control=sg_control(SgPushEffector(1, SG_PUSH)))
        push_sc = sg_config([sg_sat(1, 0.0; control=sg_control(SgPushEffector(1, SG_PUSH)))])
        base = sg_config([sg_sat(1, 0.0)])
        u_cfg = final_state(push_cfg)
        @test final_state(push_sc) == u_cfg
        @test u_cfg != final_state(base)
        # The caller's configuration is untouched, so a second run still applies the effector.
        @test length(push_sc.dynamics_model.spacecraft[1].control.control_effectors) == 1
        @test final_state(push_sc) == u_cfg
        # No per-spacecraft GNC: the flattening is the identity.
        @test SpaceAGORA.SimulationEngine._flatten_spacecraft_gnc(base) === base
        flat = SpaceAGORA.SimulationEngine._flatten_spacecraft_gnc(push_sc)
        @test flat.control_model.control_effectors isa Tuple{SgPushEffector}
        @test SpaceAGORA.SimulationEngine._flatten_spacecraft_gnc(flat) === flat
    end

    @testset "binds to the spacecraft's index" begin
        sats_decl = [sg_sat(1, 0.0), sg_sat(2, 40.0; control=sg_control(SgPushEffector(1, SG_PUSH)))]
        flat = SpaceAGORA.SimulationEngine._flatten_spacecraft_gnc(sg_config(sats_decl; control=sg_control()))
        @test [e.sat_idx for e in flat.control_model.control_effectors] == [2]
        u_decl = final_state(sg_config(sats_decl))
        u_cfg = final_state(sg_config([sg_sat(1, 0.0), sg_sat(2, 40.0)]; control=sg_control(SgPushEffector(2, SG_PUSH))))
        @test u_decl == u_cfg
        # Configuration-level effectors come first, then spacecraft effectors in spacecraft order.
        mixed = sg_config([sg_sat(1, 0.0; control=sg_control(SgPushEffector(9, SG_PUSH))), sg_sat(2, 40.0; control=sg_control(SgPushEffector(9, SG_PUSH)))];
                          control=sg_control(SgPushEffector(1, SG_PUSH); rate=2.0))
        flat = SpaceAGORA.SimulationEngine._flatten_spacecraft_gnc(mixed)
        @test [e.sat_idx for e in flat.control_model.control_effectors] == [1, 1, 2]
        @test flat.control_model.control_rates == [2.0, 1.0, 1.0]
    end

    @testset "built-in effectors" begin
        m = MagneticMomentumManagerModel(sat_idx=1)
        @test bind_spacecraft(m, 1) === m
        @test bind_spacecraft(m, 3).sat_idx == 3
        @test m.sat_idx == 1
        r = RobotArmControlEffector()
        @test bind_spacecraft(r, 2).spacecraft_idx == 2
        # Multi-vehicle and unindexed effectors cannot be declared on a spacecraft.
        @test_throws ArgumentError bind_spacecraft(RPOMPCControlModel(), 1)
        err = try bind_spacecraft(RPOMPCControlModel(), 1) catch e e end
        @test occursin("chaser", sprint(showerror, err)) && occursin("configuration level", sprint(showerror, err))
        thr = BaseThrusterModel(thrust=[0.0], direction=[0.0], Δv=[0.0], start_burn_time=[0.0], stop_burn_time=[0.0], Isp=[300.0])
        err = try bind_spacecraft(thr, 1) catch e e end
        @test err isa ArgumentError && occursin("every spacecraft", sprint(showerror, err))
        sc = sg_sat(1, 0.0; control=sg_control(RPOMPCControlModel()))
        @test_throws ArgumentError run_simulation(sg_config([sc]))
    end

    @testset "built-in binding methods" begin
        # Index-field effectors: new index, original untouched, same index returns the same object.
        m = MagneticMomentumManagerModel(sat_idx=1, mu_gain=0.02)
        mb = bind_spacecraft(m, 2)
        @test mb isa MagneticMomentumManagerModel && mb.sat_idx == 2 && m.sat_idx == 1 && mb !== m
        @test mb.mu_gain == 0.02
        @test bind_spacecraft(mb, 2) === mb
        r = RobotArmControlEffector(spacecraft_idx=1, joint_kp=3.0)
        rb = bind_spacecraft(r, 2)
        @test rb isa RobotArmControlEffector && rb.spacecraft_idx == 2 && r.spacecraft_idx == 1 && rb !== r
        @test rb.joint_kp == 3.0
        @test bind_spacecraft(rb, 2) === rb

        # Apollo guidance and control, state sized for two spacecraft.
        tg, ta = apollo11_descent_targets()
        radius = SG_PLANET.Rp_e - 2_000.0
        cfg = ApolloDescentConfig(reference_radius_m=radius, site_lat_deg=30.0, site_lon_deg=0.0, braking=tg, approach=ta)
        state = ApolloDescentState(2)
        g = ApolloDescentGuidanceModel(cfg, state; spacecraft_indices=(1,))
        c = ApolloDescentControlModel(ApolloDescentControlConfig(), cfg, state; spacecraft_indices=(1,))
        for (model, ctor_state) in ((g, state), (c, state))
            mbound = bind_spacecraft(model, 2)
            @test typeof(mbound) == typeof(model) && mbound !== model
            @test mbound.spacecraft_indices == (2,) && model.spacecraft_indices == (1,)
            @test mbound.state === ctor_state
            @test bind_spacecraft(mbound, 2) === mbound
            # State sized for two spacecraft cannot serve a third.
            err = try bind_spacecraft(model, 3) catch e e end
            @test err isa ArgumentError && occursin("spacecraft_indices", sprint(showerror, err))
        end
        @test bind_spacecraft(c, 2).actuators === c.actuators
    end

    @testset "constellation ensemble" begin
        sats = [sg_sat(11, 0.0), sg_sat(22, 40.0), sg_sat(33, 80.0; control=sg_control(SgPushEffector(1, SG_PUSH)))]
        res = run_constellation_ensemble(sg_config(sats); return_solution=true)   # no allow_gnc_effectors
        @test isempty(res.failed)
        u_member = collect(res.samples[3].value.u[end])
        u_solo = final_state(sg_config([sg_sat(33, 80.0)]; control=sg_control(SgPushEffector(1, SG_PUSH))))
        @test maximum(abs.(u_member .- u_solo)) <= 1e-9
        # The unaffected members match a run without GNC.
        u_plain = final_state(sg_config([sg_sat(11, 0.0)]))
        @test maximum(abs.(collect(res.samples[1].value.u[end]) .- u_plain)) <= 1e-9
        # Configuration-level effectors still need the flag.
        @test_throws ArgumentError run_constellation_ensemble(sg_config(sats; control=sg_control(SgPushEffector(1, SG_PUSH))))
    end
end

# Both scheduling lanes must survive flattening, including mutable-state isolation.
struct SgGuidanceProbe <: SpaceAGORA.AbstractGuidanceModel
    sat_idx::Int
    visits::Vector{Int}
end
struct SgNavigationProbe
    sat_idx::Int
    visits::Vector{Int}
end
SpaceAGORA.bind_spacecraft(m::SgGuidanceProbe, i::Int) = SgGuidanceProbe(i, m.visits)
SpaceAGORA.bind_spacecraft(m::SgNavigationProbe, i::Int) = SgNavigationProbe(i, m.visits)
function SpaceAGORA.SimulationModel.GuidanceHooks.calcGuidanceEffect!(m::SgGuidanceProbe, u, p, t::Float64, i::Int)
    i == m.sat_idx && push!(m.visits, i)
    return nothing
end
function SpaceAGORA.SimulationModel.NavigationHooks.calcNavigationEffect!(m::SgNavigationProbe, u, p, t::Float64, i::Int)
    i == m.sat_idx && push!(m.visits, i)
    return nothing
end

@testset "Per-spacecraft GNC scheduling and ownership" begin
    g = SgGuidanceProbe(9, Int[])
    n = SgNavigationProbe(9, Int[])
    gm = GuidanceModel((g,), [2.0]); nm = NavigationModel((n,), [3.0])
    args = sg_config([sg_sat(1, 0.0), sg_sat(2, 40.0; guidance=gm, navigation=nm)])
    flat = SpaceAGORA.SimulationEngine._flatten_spacecraft_gnc(args)
    @test flat.guidance_model.guidance_effectors isa Tuple{SgGuidanceProbe}
    @test flat.navigation_model.navigation_effectors isa Tuple{SgNavigationProbe}
    @test flat.guidance_model.guidance_effectors[1].sat_idx == 2
    @test flat.navigation_model.navigation_effectors[1].sat_idx == 2
    @test flat.guidance_model.guidance_rates == [2.0]
    @test flat.navigation_model.navigation_rates == [3.0]
    @test flat.guidance_model.guidance_rates !== gm.guidance_rates
    @test flat.navigation_model.navigation_rates !== nm.navigation_rates
    @test SpaceAGORA.SimulationEngine._flatten_spacecraft_gnc(flat) === flat
    @test (run_simulation(args); true)
    @test isempty(g.visits) && isempty(n.visits)
    @test (run_simulation(args; isolate_state=false); true)
    @test !isempty(g.visits) && all(==(2), g.visits)
    @test !isempty(n.visits) && all(==(2), n.visits)
    @test length(g.visits) > length(n.visits)
    @test g.sat_idx == 9 && n.sat_idx == 9
    @test args.dynamics_model.spacecraft[2].guidance === gm
    @test args.dynamics_model.spacecraft[2].navigation === nm

    # Built-in rebinds copy scalar fields but retain nested mutable members.
    arm = RobotArmControlEffector(spacecraft_idx=1)
    bound = bind_spacecraft(arm, 2)
    bound.updated_at_s = 5.0
    push!(bound.held.joint_torque_nm, 7.0)
    @test arm.updated_at_s == 0.0
    @test bound.held === arm.held
    @test arm.held.joint_torque_nm == [7.0]
end

@testset "Per-spacecraft GNC returned configuration" begin
    g = SgGuidanceProbe(9, Int[])
    n = SgNavigationProbe(9, Int[])
    args = sg_config([sg_sat(1, 0.0), sg_sat(2, 40.0;
        guidance=GuidanceModel((g,), [2.0]), navigation=NavigationModel((n,), [3.0]))])
    result = run_simulation(args; return_results=true, return_solution=true)
    @test result isa SimulationResults
    @test result.solution !== nothing
    @test size(result.table, 1) > 0
    @test isempty(result.files)
    ran_g = only(result.configuration.guidance_model.guidance_effectors)
    ran_n = only(result.configuration.navigation_model.navigation_effectors)
    @test ran_g !== g && ran_n !== n
    @test !isempty(ran_g.visits) && all(==(2), ran_g.visits)
    @test !isempty(ran_n.visits) && all(==(2), ran_n.visits)
    @test isempty(g.visits) && isempty(n.visits)
    @test all(sc -> isempty(sc.guidance.guidance_effectors) &&
        isempty(sc.navigation.navigation_effectors), result.configuration.dynamics_model.spacecraft)
    rerun = run_simulation(result.configuration; return_results=true)
    @test length(rerun.configuration.guidance_model.guidance_effectors) == 1
    @test length(rerun.configuration.navigation_model.navigation_effectors) == 1
    @test_throws ArgumentError run_simulation(args; return_results=true, return_solver_metadata=true)
end
