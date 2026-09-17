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

# Explicit axes can retain the same property names while changing storage order.
function _rhs_direct_assembly_reaxis(u; pos=1:3, vel=4:6, mass=7, heat_loads=8:8)
    n_sats = length(u.sc)
    stride = length(ComponentArrays.getdata(u.sc[1]))
    fields = (pos=pos, vel=vel, mass=mass, heat_loads=heat_loads)
    axis = ComponentArrays.Axis(sc=ComponentArrays.ViewAxis(
        1:length(u), ComponentArrays.PartitionedAxis(stride, fields)))
    return ComponentArray(copy(ComponentArrays.getdata(u)), axis)
end

@testset "RHS direct layout probe rejects noncanonical field offsets" begin
    u = RHSDA_SE.build_initial_conditions(_rhs_direct_assembly_config(n_sats=3))
    du = zero(u)
    @test RHSDA_SE._probe_rhs_final_assembly_direct_stride(u.sc, du.sc, u, du) == 8
    for bad in (
        _rhs_direct_assembly_reaxis(u; pos=4:6, vel=1:3),
        _rhs_direct_assembly_reaxis(u; mass=8, heat_loads=7:7),
        _rhs_direct_assembly_reaxis(u; pos=1:2:5, vel=2:2:6),
    )
        @test propertynames(bad.sc[1]) == propertynames(u.sc[1])
        @test RHSDA_SE._probe_rhs_final_assembly_direct_stride(bad.sc, du.sc, bad, du) == 0
        @test RHSDA_SE._probe_rhs_final_assembly_direct_stride(u.sc, bad.sc, u, bad) == 0
    end

    # First/last starts and total length can agree with stride 8 even when the
    # two intermediate satellites have different body counts (2 and 0).
    # Construct the child views directly to exercise that probe boundary.
    four = RHSDA_SE.build_initial_conditions(_rhs_direct_assembly_config(n_sats=4))
    four_du = zero(four)
    child_views = ComponentVector[]
    for (start, n_heat) in ((1, 1), (9, 2), (18, 0), (25, 1))
        width = 7 + n_heat
        push!(child_views, ComponentArray(
            view(ComponentArrays.getdata(four), start:(start + width - 1)),
            ComponentArrays.Axis(pos=1:3, vel=4:6, mass=7, heat_loads=8:width)))
    end
    @test RHSDA_SE._probe_rhs_final_assembly_direct_stride(
        child_views, four_du.sc, four, four_du) == 0
end

@testset "RHS direct layout cache follows both state and derivative layouts" begin
    args = _rhs_direct_assembly_config(n_sats=3)
    u = RHSDA_SE.build_initial_conditions(args)
    du = zero(u)
    p = ODEParams(n_sats=3, args=args)
    env = withenv("SPACEAGORA_RHS_FINAL_ASSEMBLY_DIRECT_LAYOUT" => "1") do
        RHSDA_SE._snapshot_rhs_plan_env_config()
    end
    totals = ones(Float64, 6, 3)
    assign!(out, state) = RHSDA_SE._try_assign_flat_translational_rhs_direct_layout!(
        out, state, p, totals, env, :full, 1)

    @test assign!(du, u)
    @test assign!(zero(u), copy(u)) # Same axes with fresh backing buffers.
    bad = _rhs_direct_assembly_reaxis(u; pos=4:6, vel=1:3)
    fill!(du, 123.0)
    @test !assign!(du, bad)
    @test all(==(123.0), ComponentArrays.getdata(du))
    @test assign!(du, u) # A cached rejection must not disable a valid layout.
    fill!(bad, 456.0)
    @test !assign!(bad, u)
    @test all(==(456.0), ComponentArrays.getdata(bad))
    @test assign!(du, u)

    # Array storage and its length are also part of the cache boundary.
    storage = zeros(Float64, 2 * length(u))
    strided = ComponentArray(view(storage, 1:2:length(storage)), ComponentArrays.getaxes(u))
    @test !assign!(strided, u)
    @test all(iszero, storage)
    @test assign!(du, u)
    shorter = copy(u)
    resize!(ComponentArrays.getdata(shorter), length(shorter) - 1)
    @test !assign!(du, shorter)
    @test assign!(du, u)
end

@testset "RHS direct final assembly matches the generic translational assembly" begin
    args = _rhs_direct_assembly_config(n_sats=3)
    u = RHSDA_SE.build_initial_conditions(args)
    p = ODEParams(n_sats=3, args=args)
    p.is_active[2] = false
    totals = reshape(collect(1.0:18.0), 6, 3)
    env = withenv("SPACEAGORA_RHS_FINAL_ASSEMBLY_DIRECT_LAYOUT" => "1") do
        RHSDA_SE._snapshot_rhs_plan_env_config()
    end
    for kind in (:full, :explicit, :slow), mass in (500.0, 0.0, eps(Float64), NaN, Inf, -500.0)
        u.sc[1].mass = mass
        direct = fill!(zero(u), 99.0)
        generic = fill!(zero(u), 99.0)
        for i in eachindex(p.is_active)
            if !p.is_active[i]
                generic.sc[i] .= 0.0
                continue
            end
            forces = view(totals, 1:3, i)
            if kind == :slow
                SpaceAGORA.SimulationModel.DynamicsTranslational.assign_slow_translational_rhs!(generic.sc[i], u.sc[i], forces)
            else
                SpaceAGORA.SimulationModel.DynamicsTranslational.assign_full_translational_rhs!(generic.sc[i], u.sc[i], forces, 0.0)
            end
            generic.sc[i].heat_loads .= 0.0
        end
        @test RHSDA_SE._try_assign_flat_translational_rhs_direct_layout!(
            direct, u, p, totals, env, kind, 2)
        @test isequal(ComponentArrays.getdata(direct), ComponentArrays.getdata(generic))
    end
end

@testset "Flat RHS preserves derivatives when direct assembly is enabled or declined" begin
    args = _rhs_direct_assembly_config(n_sats=3)
    canonical = RHSDA_SE.build_initial_conditions(args)
    noncanonical = _rhs_direct_assembly_reaxis(canonical; pos=4:6, vel=1:3)
    plan = (
        mode=:flat_constellation_effector_queue, allotment=1, scheduler=:static,
        dominant_axis=:effector, policy_applied=false,
        effector_decision=(use_threads=false, allotment=1, mode=:off, policy_applied=false),
    )
    for u in (canonical, noncanonical), kind in (:full, :explicit, :slow)
        outputs = ComponentVector[]
        for enabled in (false, true)
            p = ODEParams(n_sats=3, args=args)
            p.is_active[2] = false
            p.shared_buffers.rhs_env_config[] = withenv(
                "SPACEAGORA_RHS_FINAL_ASSEMBLY_DIRECT_LAYOUT" => (enabled ? "1" : "0"),
            ) do
                RHSDA_SE._snapshot_rhs_plan_env_config()
            end
            du = fill!(zero(u), 99.0)
            RHSDA_SE._spacecraft_dynamics_flat_constellation_effector_queue!(
                du, u, p, 0.0, plan; rhs_kind=kind,
                partition=kind == :explicit ? :explicit : nothing,
            )
            expected_status = !enabled ? Int8(0) : (u === canonical ? Int8(1) : Int8(-1))
            @test p.shared_buffers.rhs_final_assembly_direct_layout_status[] == expected_status
            @test all(iszero, ComponentArrays.getdata(du.sc[2]))
            push!(outputs, du)
        end
        @test isequal(ComponentArrays.getdata(outputs[1]), ComponentArrays.getdata(outputs[2]))
    end
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
