using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel
using ComponentArrays

const RHSDA_SE = SpaceAGORA.SimulationEngine

function _rhs_direct_assembly_config(; n_sats::Int=3, effectors::Tuple=(InverseSquaredGravityModel(),))
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
        dynamics_model=DynamicsModel(spacecraft, effectors),
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

# The harmonics-only flat route fuses the harmonics pre-pass, the reduction and
# the final assembly into one parallel region (each worker finishes the
# spacecraft range its slice covers). It must write exactly what the unfused
# composition writes -- pre-pass and reduction, then every spacecraft assembled
# from the totals -- for both assemblies and both pre-pass dispatches (persistent
# pool, spin barrier), with inactive spacecraft at the ends and in the gaps
# between slices.
@testset "Fused harmonics RHS matches the unfused route bit for bit" begin
    planet = Earth()
    gravity_file = joinpath(normpath(joinpath(@__DIR__, "..", "..", "..")),
        "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    model = GravitationalHarmonicsModel(20, 20, gravity_file, planet)
    n = 37
    args = _rhs_direct_assembly_config(n_sats=n, effectors=(model,))
    u = RHSDA_SE.build_initial_conditions(args)
    plan = (
        mode=:flat_constellation_effector_queue, allotment=max(2, Threads.nthreads()),
        scheduler=:static, dominant_axis=:flat_effector, policy_applied=false,
        effector_decision=(use_threads=false, allotment=1, mode=:off, policy_applied=false),
    )
    PP = SpaceAGORA.SimulationModel.ParallelPolicy
    for spin in (false, true), direct in (false, true), inactive in (Int[], [1, 2, 9, 10, 11, 36, 37], [5])
        fresh_p() = begin
            p = ODEParams(n_sats=n, args=args)
            p.is_active .= true
            p.is_active[inactive] .= false
            p.shared_buffers.rhs_env_config[] = withenv(
                "SPACEAGORA_RHS_FINAL_ASSEMBLY_DIRECT_LAYOUT" => (direct ? "1" : "0"),
                "SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER" => (spin ? "1" : "0"),
            ) do
                RHSDA_SE._snapshot_rhs_plan_env_config()
            end
            p
        end
        p_fused = fresh_p()
        fused = fill!(zero(u), 99.0)
        RHSDA_SE._spacecraft_dynamics_flat_constellation_effector_queue!(
            fused, u, p_fused, 30.0, plan; rhs_kind=:full)

        p_ref = fresh_p()
        ref = fill!(zero(u), 99.0)
        RHSDA_SE._accumulate_dynamic_effectors_flat_batch!(u.sc, p_ref, 30.0, (model,), plan)
        totals = p_ref.shared_buffers.rhs_flat_effector_totals[]
        env = p_ref.shared_buffers.rhs_env_config[]
        stride = RHSDA_SE._flat_translational_direct_layout_stride(ref, u, p_ref, env, :full)
        @test (stride > 0) == direct
        for i in 1:n
            RHSDA_SE._assemble_flat_satellite!(
                ref.sc, u.sc, ComponentArrays.getdata(ref), ComponentArrays.getdata(u),
                p_ref, 30.0, totals, args.dynamics_model.spacecraft,
                p_ref.shared_buffers.debug_control[], stride, :full, i)
        end
        @test !any(==(99.0), ComponentArrays.getdata(fused))
        @test isequal(ComponentArrays.getdata(fused), ComponentArrays.getdata(ref))
        @test isequal(p_fused.shared_buffers.rhs_flat_effector_totals[][:, 1:n], totals[:, 1:n])
    end
    # Stop the spin-barrier workers so they do not hold threads for later tests.
    PP._destroy_persistent_foreach_scope!(PP._active_policy_scope_id())
end

# The per-region width search halves the width while the region keeps within
# the stop ratio of its best time, then keeps the fastest width only if it beats
# the full allotment by the margin.
@testset "Per-region width search" begin
    T = SpaceAGORA.SimulationModel.RhsRegionTuner
    feed!(t, cost) = for _ in 1:RHSDA_SE._REGION_SAMPLES
        RHSDA_SE._region_tuner_observe!(t, Int64(cost(t.width)))
    end
    # Dispatch dominated: every halving is faster, so the search runs to serial.
    t = T(16)
    while t.searching
        feed!(t, w -> 1000 + 100 * w)
    end
    @test t.width == 1
    # Work dominated: halving doubles the time, so it stops after one probe.
    t = T(16)
    while t.searching
        feed!(t, w -> 160_000 ÷ w)
    end
    @test t.width == 16
    # Interior optimum of work/w + d*w at w = 4.
    t = T(32)
    while t.searching
        feed!(t, w -> 16_000 ÷ w + 1000 * w)
    end
    @test t.width == 4
    # A narrower width within the margin does not displace the allotment.
    t = T(8)
    while t.searching
        feed!(t, w -> w == 8 ? 1000 : 980)
    end
    @test t.width == 8
    # The search re-runs after the reprobe interval.
    for _ in 1:RHSDA_SE._REGION_REPROBE
        RHSDA_SE._region_tuner_observe!(t, Int64(1))
    end
    @test t.searching && t.width == 8
end

@testset "Per-region width search labels its observations" begin
    T = SpaceAGORA.SimulationModel.RhsRegionTuner
    fake(tuning) = (; shared_buffers=(; rhs_region_tuning=Ref(tuning), rhs_region_tuners=Dict{Symbol, T}()))
    p = fake(true)
    # The first samples run at the full allotment; once the search narrows the
    # region, its calls must not be recorded as full-width observations.
    for _ in 1:RHSDA_SE._REGION_SAMPLES
        RHSDA_SE._rhs_region(w -> nothing, p, :probe, 8)
        @test RHSDA_SE._rhs_region_ran_full(p, :probe)
    end
    RHSDA_SE._rhs_region(w -> nothing, p, :probe, 8)
    @test p.shared_buffers.rhs_region_tuners[:probe].last_width < 8
    @test !RHSDA_SE._rhs_region_ran_full(p, :probe)
    # Search off, or a region never run: every observation is at full width.
    @test RHSDA_SE._rhs_region_ran_full(fake(false), :probe)
    @test RHSDA_SE._rhs_region_ran_full(fake(true), :never_run)
    # An unrecognised setting is an error, not a silent "auto".
    withenv("SPACEAGORA_RHS_REGION_WIDTH" => "8") do
        @test_throws ArgumentError RHSDA_SE._rhs_region_width_enabled()
    end
end

# With the search on, the width of every region changes from call to call
# while it searches; the derivative must not.
@testset "Per-region width search leaves the RHS bit-identical" begin
    planet = Earth()
    gravity_file = joinpath(normpath(joinpath(@__DIR__, "..", "..", "..")),
        "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    model = GravitationalHarmonicsModel(20, 20, gravity_file, planet)
    n = 37
    args = _rhs_direct_assembly_config(n_sats=n, effectors=(model, InverseSquaredGravityModel()))
    u = RHSDA_SE.build_initial_conditions(args)
    plan = (
        mode=:flat_constellation_effector_queue, allotment=max(2, Threads.nthreads()),
        scheduler=:static, dominant_axis=:flat_effector, policy_applied=false,
        effector_decision=(use_threads=false, allotment=1, mode=:off, policy_applied=false),
    )
    for effectors in ((model,), (model, InverseSquaredGravityModel()))
        a = _rhs_direct_assembly_config(n_sats=n, effectors=effectors)
        fresh_p() = begin
            p = ODEParams(n_sats=n, args=a)
            p.is_active .= true
            p.shared_buffers.rhs_env_config[] = withenv(
                "SPACEAGORA_RHS_REDUCTION_MIN_SATS_PER_WORKER" => "1",
            ) do
                RHSDA_SE._snapshot_rhs_plan_env_config()
            end
            p
        end
        p_ref = fresh_p()
        ref = zero(u)
        withenv("SPACEAGORA_RHS_REDUCTION_MIN_SATS_PER_WORKER" => "1") do
            RHSDA_SE._spacecraft_dynamics_flat_constellation_effector_queue!(ref, u, p_ref, 30.0, plan; rhs_kind=:full)
        end
        p = fresh_p()
        p.shared_buffers.rhs_region_tuning[] = true
        widths = Set{Int}()
        for _ in 1:60
            du = fill!(zero(u), 99.0)
            withenv("SPACEAGORA_RHS_REDUCTION_MIN_SATS_PER_WORKER" => "1") do
                RHSDA_SE._spacecraft_dynamics_flat_constellation_effector_queue!(du, u, p, 30.0, plan; rhs_kind=:full)
            end
            @test isequal(ComponentArrays.getdata(du), ComponentArrays.getdata(ref))
            for t in values(p.shared_buffers.rhs_region_tuners)
                push!(widths, t.width)
            end
        end
        @test !isempty(p.shared_buffers.rhs_region_tuners)
        Threads.nthreads() > 1 && @test length(widths) > 1
    end
end

# SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER: 1/0 force the barrier, `auto` (the
# default) probes both arms during the solve and keeps the faster.
@testset "Spin barrier env modes" begin
    PP = SpaceAGORA.SimulationModel.ParallelPolicy
    for (raw, mode) in ((nothing, :auto), ("auto", :auto), (" AUTO ", :auto), ("1", :on), ("on", :on),
                        ("true", :on), ("0", :off), ("off", :off), ("false", :off))
        withenv("SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER" => raw) do
            @test PP.harmonics_batch_spin_barrier_mode() === mode
            # Only a forced `on` counts as enabled for the pre-solve paths.
            @test PP.harmonics_batch_spin_barrier_enabled() == (mode === :on)
            @test RHSDA_SE._snapshot_rhs_plan_env_config().harmonics_batch_spin_barrier === mode
        end
    end
    withenv("SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER" => "sometimes") do
        @test_throws ArgumentError PP.harmonics_batch_spin_barrier_mode()
    end
    # Resolution per call: forced modes ignore the tuner; auto follows its arm,
    # and is off where there is no tuner.
    T = SpaceAGORA.SimulationModel.SpinBarrierTuner
    env_for(raw) = withenv(() -> RHSDA_SE._snapshot_rhs_plan_env_config(),
        "SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER" => raw)
    for on in (false, true)
        t = T(); t.on = on
        p = (; shared_buffers=(; spin_barrier_tuner=t))
        @test RHSDA_SE._spin_barrier_active(p, env_for("1"))
        @test !RHSDA_SE._spin_barrier_active(p, env_for("0"))
        @test RHSDA_SE._spin_barrier_active(p, env_for("auto")) == on
    end
    @test !RHSDA_SE._spin_barrier_active(nothing, env_for("auto"))
    # The barrier lifts the flat route's floor; without it the floor applies.
    env = env_for("auto")
    @test RHSDA_SE._rhs_flat_default_admits(env, 64, 32, true)
    @test !RHSDA_SE._rhs_flat_default_admits(env, 64, 32, false)
end

@testset "Spin barrier auto: decision logic" begin
    T = SpaceAGORA.SimulationModel.SpinBarrierTuner
    SE = RHSDA_SE
    span = SE._SPIN_WINDOWS * SE._SPIN_SAMPLES
    fake() = (; shared_buffers=(; rhs_plan_step_cache=Ref{Any}(:cached)))
    observe!(p, t, ns) = SE._spin_barrier_observe!(p, t, Int64(ns))
    # Run calls costing cost(t) until the current probe decides; returns the call count.
    probe!(p, t, cost) = (n = 0; while t.probing; observe!(p, t, cost(t)); n += 1; end; n)
    arm(on_ns, off_ns) = t -> t.on ? on_ns : off_ns
    fresh() = (t = T(); SE._spin_barrier_tuner_reset!(t, true); t)

    # Starts off and times the incumbent over the full span before switching.
    p = fake(); t = fresh()
    @test !t.on && !t.incumbent && t.probing
    for _ in 1:(SE._SPIN_SETTLE + span - 1)
        observe!(p, t, 1000)
    end
    @test !t.on
    observe!(p, t, 1000)
    @test t.on && t.probing                                    # now timing the challenger
    @test p.shared_buffers.rhs_plan_step_cache[] === nothing   # cached plan dropped
    # A challenger ahead by more than the margin takes over after its full span.
    n = probe!(p, t, arm(500, 1000))
    @test n == SE._SPIN_SETTLE + span
    @test t.on && t.incumbent
    @test t.switches == 1 && t.probes == 1 && t.period == SE._SPIN_REPROBE

    # It reprobes after the period; if the barrier now loses, the solve goes back.
    for _ in 1:t.period
        observe!(p, t, 500)
    end
    @test t.probing
    probe!(p, t, arm(2000, 1000))
    @test !t.on && !t.incumbent && t.switches == 2 && t.period == SE._SPIN_REPROBE

    # A clearly slower challenger is cut after one window, and the period
    # doubles while the choice holds.
    p = fake(); t = fresh()
    n = probe!(p, t, arm(9000, 1000))
    @test n == 2 * SE._SPIN_SETTLE + span + SE._SPIN_SAMPLES
    @test !t.on && !t.incumbent && t.switches == 0
    @test t.period == 2 * SE._SPIN_REPROBE
    @test p.shared_buffers.rhs_plan_step_cache[] === nothing
    for k in 2:10
        for _ in 1:t.period
            observe!(p, t, 1000)
        end
        probe!(p, t, arm(9000, 1000))
    end
    @test t.period == SE._SPIN_REPROBE_MAX && !t.on

    # Faster by less than the margin: the incumbent stays.
    p = fake(); t = fresh()
    probe!(p, t, arm(980, 1000))
    @test !t.incumbent && !t.on

    # Means, not medians: a challenger fast on most calls but with rare long
    # stalls that make it slower overall does not take over.
    p = fake(); t = fresh()
    k = Ref(0)
    probe!(p, t, t -> t.on ? (k[] += 1; k[] % 64 == 0 ? 60_000 : 500) : 1000)
    @test !t.incumbent

    # Settle calls after a switch are not samples: a slow start on the new arm
    # (pool start-up) does not decide the probe.
    p = fake(); t = fresh()
    probe!(p, t, t -> t.skip > 0 ? 10^9 : (t.on ? 500 : 1000))
    @test t.incumbent
end

@testset "Spin barrier auto: observations labelled by arm" begin
    T = SpaceAGORA.SimulationModel.SpinBarrierTuner
    R = SpaceAGORA.SimulationModel.RhsRegionTuner
    t = T(); RHSDA_SE._spin_barrier_tuner_reset!(t, true)
    p = (; shared_buffers=(; spin_barrier_tuner=t, rhs_region_tuning=Ref(false),
        rhs_region_tuners=Dict{Symbol, R}(), rhs_plan_step_cache=Ref{Any}(nothing)))
    @test RHSDA_SE._rhs_region_ran_full(p, :rhs_harmonics_batch)        # incumbent arm
    t.on = true                                                           # challenger arm
    @test !RHSDA_SE._rhs_region_ran_full(p, :rhs_harmonics_batch)
    @test !RHSDA_SE._rhs_region_ran_full(p, :rhs_flat_queue)
    t.incumbent = true                                                    # settled on
    @test RHSDA_SE._rhs_region_ran_full(p, :rhs_harmonics_batch_spin)
    t.enabled = false; t.incumbent = false                                # forced modes: no probe
    @test RHSDA_SE._rhs_region_ran_full(p, :rhs_harmonics_batch_spin)
end

# The arm (dispatch, and with it the default route's floors) changes from call
# to call under auto; the derivative and the solve must not.
@testset "Spin barrier on/off/auto leave the RHS and the solve bit-identical" begin
    PP = SpaceAGORA.SimulationModel.ParallelPolicy
    planet = Earth()
    gravity_file = joinpath(normpath(joinpath(@__DIR__, "..", "..", "..")),
        "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    model = GravitationalHarmonicsModel(20, 20, gravity_file, planet)
    n = 37
    args = _rhs_direct_assembly_config(n_sats=n, effectors=(model,))
    u = RHSDA_SE.build_initial_conditions(args)
    fresh_p(raw) = begin
        p = ODEParams(n_sats=n, args=args)
        p.is_active .= true
        p.shared_buffers.rhs_env_config[] = withenv(() -> RHSDA_SE._snapshot_rhs_plan_env_config(),
            "SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER" => raw)
        p
    end
    rhs(p) = (du = fill!(zero(u), 99.0); RHSDA_SE.spacecraft_dynamics!(du, u, p, 30.0); ComponentArrays.getdata(du))
    ref = rhs(fresh_p("0"))
    @test !any(==(99.0), ref)
    @test isequal(rhs(fresh_p("1")), ref)
    p = fresh_p("auto")
    RHSDA_SE._spin_barrier_tuner_reset!(p.shared_buffers.spin_barrier_tuner, true)
    arms = Set{Bool}()
    calls = RHSDA_SE._SPIN_SETTLE + RHSDA_SE._SPIN_WINDOWS * RHSDA_SE._SPIN_SAMPLES + 40
    same = true
    for k in 1:calls
        push!(arms, p.shared_buffers.spin_barrier_tuner.on)
        same &= isequal(rhs(p), ref)
    end
    @test same
    @test arms == Set([false, true])
    PP._destroy_persistent_foreach_scope!(PP._active_policy_scope_id())

    finals = Dict{String, Vector{Float64}}()
    for raw in ("0", "1", "auto")
        sol = withenv("SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER" => raw) do
            SpaceAGORA.run_simulation(_rhs_direct_assembly_config(n_sats=n, effectors=(model,));
                return_solution=true)
        end
        s = sol isa Tuple ? last(sol) : sol
        finals[raw] = copy(ComponentArrays.getdata(s.u[end]))
    end
    @test isequal(finals["1"], finals["0"])
    @test isequal(finals["auto"], finals["0"])
end
