@testset "Thruster Edge Cases (AI Second-Pass)" begin
    @testset "BaseThrusterModel Validation" begin
        model_ok = BaseThrusterModel(
            thrust=[1.0, 2.0],
            direction=[0.0, π],
            Δv=[10.0, 20.0],
            start_burn_time=[0.0, 0.0],
            stop_burn_time=[100.0, 120.0],
            Isp=[300.0, 300.0]
        )
        @test length(model_ok.thrust) == 2
        @test length(model_ok.direction) == 2

        @test_throws ArgumentError BaseThrusterModel(
            thrust=[1.0, 2.0],
            direction=[0.0],
            Δv=[10.0, 20.0],
            start_burn_time=[0.0, 0.0],
            stop_burn_time=[100.0, 120.0],
            Isp=[300.0, 300.0]
        )
    end

    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=600.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    p = ODEParams{1}(args=args)
    u = build_initial_conditions(args).sc[1]

    _burn_plan(p_, sat_idx) = p_.shared_buffers.maneuver_burn_plans[sat_idx]
    _clear_burn_plan!(p_, sat_idx) = (p_.shared_buffers.maneuver_burn_plans[sat_idx] = PropulsiveBurnPlan())
    _expected_burn_duration(mass_kg, thrust_n, isp_s, delta_v_mps) =
        mass_kg * isp_s * 9.80665 * (1.0 - exp(-delta_v_mps / (isp_s * 9.80665))) / thrust_n

    @testset "calcControlForceTorque" begin
        model = make_base_thruster_model(thrust=2.0, direction=0.0, start_burn_time=50.0, stop_burn_time=150.0)

        force_pre, torque_pre = calcControlForceTorque(model, u, p, 1, 49.9)
        @test force_pre == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test torque_pre == SVector{3, Float64}(0.0, 0.0, 0.0)

        expected_dir = normalize(SVector{3, Float64}(u.vel))
        force_start, torque_start = calcControlForceTorque(model, u, p, 1, 50.0)
        @test isapprox(force_start, 2.0 .* expected_dir; atol=1e-12, rtol=0.0)
        @test torque_start == SVector{3, Float64}(0.0, 0.0, 0.0)

        force_stop, _ = calcControlForceTorque(model, u, p, 1, 150.0)
        @test isapprox(force_stop, 2.0 .* expected_dir; atol=1e-12, rtol=0.0)

        force_post, _ = calcControlForceTorque(model, u, p, 1, 150.1)
        @test force_post == SVector{3, Float64}(0.0, 0.0, 0.0)

        model.direction[1] = π
        force_retro, _ = calcControlForceTorque(model, u, p, 1, 100.0)
        @test dot(force_retro, expected_dir) < 0.0

        u_zero = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])
        force_zero, torque_zero = calcControlForceTorque(model, u_zero, p, 1, 100.0)
        @test force_zero == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test torque_zero == SVector{3, Float64}(0.0, 0.0, 0.0)

        force_oob, torque_oob = calcControlForceTorque(model, u, p, 2, 100.0)
        @test force_oob == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test torque_oob == SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    @testset "calcControlMassFlowRate" begin
        model = make_base_thruster_model(
            thrust=2.0,
            direction=0.0,
            start_burn_time=50.0,
            stop_burn_time=150.0,
            Isp=300.0
        )

        mdot_pre = calcControlMassFlowRate(model, u, p, 1, 49.9)
        @test mdot_pre == 0.0

        mdot_on = calcControlMassFlowRate(model, u, p, 1, 100.0)
        expected_mdot = -2.0 / (300.0 * 9.80665)
        @test isapprox(mdot_on, expected_mdot; atol=1e-14, rtol=0.0)

        mdot_post = calcControlMassFlowRate(model, u, p, 1, 150.1)
        @test mdot_post == 0.0

        model_bad_isp = make_base_thruster_model(
            thrust=2.0,
            direction=0.0,
            start_burn_time=50.0,
            stop_burn_time=150.0,
            Isp=0.0
        )
        @test calcControlMassFlowRate(model_bad_isp, u, p, 1, 100.0) == 0.0

        u_zero = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])
        @test calcControlMassFlowRate(model, u_zero, p, 1, 100.0) == 0.0
        @test calcControlMassFlowRate(model, u, p, 2, 100.0) == 0.0

        model_subtyped = TimedTangentialThrusterModel(1.0, +1.0, 10.0, 20.0)
        @test calcControlMassFlowRate(model_subtyped, u, p, 1, 15.0) == 0.0

        struct UntypedControlEffector end
        @test calcControlMassFlowRate(UntypedControlEffector(), u, p, 1, 100.0) == 0.0
    end

    @testset "calcControlEffect!" begin
        model = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        state = build_initial_conditions(args)
        _clear_burn_plan!(p, 1)
        calcControlEffect!(model, state, p, 100.0, 1)
        plan = _burn_plan(p, 1)
        @test plan.valid
        @test isfinite(plan.start_burn_s)
        @test isfinite(plan.stop_burn_s)
        @test plan.start_burn_s < plan.stop_burn_s
        @test isfinite(model.start_burn_time[1])
        @test isfinite(model.stop_burn_time[1])
        @test model.start_burn_time[1] == plan.start_burn_s
        @test model.stop_burn_time[1] == plan.stop_burn_s
        expected_burn_duration = _expected_burn_duration(state.sc[1].mass, model.thrust[1], model.Isp[1], model.Δv[1])
        @test isapprox(
            plan.stop_burn_s - plan.start_burn_s,
            expected_burn_duration;
            atol=1e-9,
            rtol=0.0
        )
        @test isapprox(plan.commanded_impulse_n_s, model.thrust[1] * expected_burn_duration; atol=1e-9, rtol=0.0)
        @test expected_burn_duration < state.sc[1].mass * model.Δv[1] / model.thrust[1]

        model_signed_prograde = make_base_thruster_model(thrust=2.0, direction=π, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        _clear_burn_plan!(p, 1)
        calcControlEffect!(model_signed_prograde, state, p, 100.0, 1)
        @test model_signed_prograde.Δv[1] == 20.0
        @test isapprox(model_signed_prograde.direction[1], 0.0; atol=1e-12, rtol=0.0)
        signed_prograde_plan = _burn_plan(p, 1)
        @test signed_prograde_plan.valid
        @test isfinite(signed_prograde_plan.start_burn_s)
        @test isfinite(signed_prograde_plan.stop_burn_s)
        @test model_signed_prograde.start_burn_time[1] == signed_prograde_plan.start_burn_s
        @test model_signed_prograde.stop_burn_time[1] == signed_prograde_plan.stop_burn_s

        model_signed_retrograde = make_base_thruster_model(thrust=2.0, direction=0.0, Δv=-20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        _clear_burn_plan!(p, 1)
        calcControlEffect!(model_signed_retrograde, state, p, 100.0, 1)
        @test model_signed_retrograde.Δv[1] == 20.0
        @test isapprox(model_signed_retrograde.direction[1], π; atol=1e-12, rtol=0.0)
        signed_retrograde_plan = _burn_plan(p, 1)
        @test signed_retrograde_plan.valid
        @test isfinite(signed_retrograde_plan.start_burn_s)
        @test isfinite(signed_retrograde_plan.stop_burn_s)
        @test model_signed_retrograde.start_burn_time[1] == signed_retrograde_plan.start_burn_s
        @test model_signed_retrograde.stop_burn_time[1] == signed_retrograde_plan.stop_burn_s

        # Pre-ignition tracking: a future scheduled window should be retimed.
        model_track = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=10_000.0, stop_burn_time=11_000.0)
        _clear_burn_plan!(p, 1)
        calcControlEffect!(model_track, state, p, 100.0, 1)
        plan_track = _burn_plan(p, 1)
        @test plan_track.valid
        @test plan_track.start_burn_s != 10_000.0
        @test plan_track.stop_burn_s != 11_000.0
        @test model_track.start_burn_time[1] == plan_track.start_burn_s
        @test model_track.stop_burn_time[1] == plan_track.stop_burn_s

        # Post-ignition lock: once burn start time has been reached, keep fixed.
        s_sched = plan.start_burn_s
        e_sched = plan.stop_burn_s
        calcControlEffect!(model, state, p, s_sched + 1e-6, 1)
        plan_locked = _burn_plan(p, 1)
        @test plan_locked.start_burn_s == s_sched
        @test plan_locked.stop_burn_s == e_sched

        sc_ineligible = make_spacecraft(ra_alt_m=600e3, rp_alt_m=400e3, ν_deg=210.0)
        args_ineligible = build_config(
            spacecraft=sc_ineligible,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_ineligible = ODEParams{1}(args=args_ineligible)
        state_ineligible = build_initial_conditions(args_ineligible)
        model_ineligible = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=11.0, stop_burn_time=22.0)
        _clear_burn_plan!(p_ineligible, 1)
        calcControlEffect!(model_ineligible, state_ineligible, p_ineligible, 100.0, 1)
        @test !_burn_plan(p_ineligible, 1).valid

        model_zero_thrust = make_base_thruster_model(thrust=0.0, Δv=20.0, start_burn_time=33.0, stop_burn_time=44.0)
        _clear_burn_plan!(p, 1)
        calcControlEffect!(model_zero_thrust, state, p, 100.0, 1)
        @test !_burn_plan(p, 1).valid

        sc_edge_block = make_spacecraft(ra_alt_m=600e3, rp_alt_m=400e3, ν_deg=180.0)
        r_edge_block, _ = orbitalelemtorv(sc_edge_block.initial_condition, EARTH)
        ei_edge_block_km = (norm(SVector{3, Float64}(r_edge_block)) - EARTH.Rp_e) / 1e3
        args_edge_block = build_config(
            spacecraft=sc_edge_block,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=ei_edge_block_km,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_edge_block = ODEParams{1}(args=args_edge_block)
        state_edge_block = build_initial_conditions(args_edge_block)
        model_edge_block = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        _clear_burn_plan!(p_edge_block, 1)
        calcControlEffect!(model_edge_block, state_edge_block, p_edge_block, 100.0, 1)
        @test !_burn_plan(p_edge_block, 1).valid

        sc_edge_allow = make_spacecraft(ra_alt_m=600e3, rp_alt_m=400e3, ν_deg=170.0)
        r_edge_allow, _ = orbitalelemtorv(sc_edge_allow.initial_condition, EARTH)
        ei_edge_allow_km = (norm(SVector{3, Float64}(r_edge_allow)) - EARTH.Rp_e) / 1e3
        args_edge_allow = build_config(
            spacecraft=sc_edge_allow,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=ei_edge_allow_km,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_edge_allow = ODEParams{1}(args=args_edge_allow)
        state_edge_allow = build_initial_conditions(args_edge_allow)
        model_edge_allow = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        _clear_burn_plan!(p_edge_allow, 1)
        calcControlEffect!(model_edge_allow, state_edge_allow, p_edge_allow, 100.0, 1)
        @test _burn_plan(p_edge_allow, 1).valid

        sc_circular = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=330.0)
        args_circular = build_config(
            spacecraft=sc_circular,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_circular = ODEParams{1}(args=args_circular)
        state_circular = build_initial_conditions(args_circular)
        model_circular = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        _clear_burn_plan!(p_circular, 1)
        calcControlEffect!(model_circular, state_circular, p_circular, 100.0, 1)
        @test _burn_plan(p_circular, 1).valid

        state_hyperbolic = build_initial_conditions(args)
        rmag = norm(SVector{3, Float64}(state_hyperbolic.sc[1].pos))
        escape_speed = sqrt(2.0 * EARTH.μ / rmag)
        vhat = normalize(SVector{3, Float64}(state_hyperbolic.sc[1].vel))
        state_hyperbolic.sc[1].vel .= (1.2 * escape_speed) .* vhat
        model_hyperbolic = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=91.0, stop_burn_time=92.0)
        _clear_burn_plan!(p, 1)
        @test_nowarn calcControlEffect!(model_hyperbolic, state_hyperbolic, p, 100.0, 1)
        @test !_burn_plan(p, 1).valid
        @test model_hyperbolic.start_burn_time[1] == -1.0
        @test model_hyperbolic.stop_burn_time[1] == -1.0

        state_near_parabolic = build_initial_conditions(args)
        rmag_parabolic = norm(SVector{3, Float64}(state_near_parabolic.sc[1].pos))
        vhat_parabolic = normalize(SVector{3, Float64}(state_near_parabolic.sc[1].vel))
        v_near_parabolic = 0.999999 * sqrt(2.0 * EARTH.μ / rmag_parabolic)
        state_near_parabolic.sc[1].vel .= v_near_parabolic .* vhat_parabolic
        model_near_parabolic = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=95.0, stop_burn_time=96.0)
        _clear_burn_plan!(p, 1)
        @test_nowarn calcControlEffect!(model_near_parabolic, state_near_parabolic, p, 100.0, 1)
        @test _burn_plan(p, 1).valid
        @test isfinite(_burn_plan(p, 1).start_burn_s)
        @test isfinite(_burn_plan(p, 1).stop_burn_s)
        @test model_near_parabolic.start_burn_time[1] == _burn_plan(p, 1).start_burn_s
        @test model_near_parabolic.stop_burn_time[1] == _burn_plan(p, 1).stop_burn_s

        state_singular = build_initial_conditions(args)
        state_singular.sc[1].pos .= 0.0
        state_singular.sc[1].vel .= 0.0
        model_singular = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=93.0, stop_burn_time=94.0)
        _clear_burn_plan!(p, 1)
        @test_nowarn calcControlEffect!(model_singular, state_singular, p, 100.0, 1)
        @test !_burn_plan(p, 1).valid
        @test model_singular.start_burn_time[1] == -1.0
        @test model_singular.stop_burn_time[1] == -1.0

        state_tiny_a = build_initial_conditions(args)
        state_tiny_a.sc[1].pos .= SVector{3, Float64}(1.0, 0.0, 0.0)
        state_tiny_a.sc[1].vel .= SVector{3, Float64}(0.0, 0.1, 0.0)
        model_tiny_a = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=97.0, stop_burn_time=98.0)
        _clear_burn_plan!(p, 1)
        @test_nowarn calcControlEffect!(model_tiny_a, state_tiny_a, p, 100.0, 1)
        @test !_burn_plan(p, 1).valid
        @test model_tiny_a.start_burn_time[1] == -1.0
        @test model_tiny_a.stop_burn_time[1] == -1.0

        sc_budget = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
        sc_budget.prop_mass = 5.0
        args_budget = build_config(
            spacecraft=sc_budget,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_budget = ODEParams{1}(args=args_budget)
        state_budget = build_initial_conditions(args_budget)
        model_budget = make_base_thruster_model(thrust=2.0, Δv=400.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        _clear_burn_plan!(p_budget, 1)
        calcControlEffect!(model_budget, state_budget, p_budget, 100.0, 1)
        @test !_burn_plan(p_budget, 1).valid

        sc_multi_1 = make_spacecraft(ra_alt_m=650e3, rp_alt_m=450e3, ν_deg=120.0)
        sc_multi_2 = make_spacecraft(ra_alt_m=700e3, rp_alt_m=500e3, ν_deg=150.0)
        args_multi = build_config_multi(
            spacecraft=[sc_multi_1, sc_multi_2],
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_multi = ODEParams{2}(args=args_multi)
        state_multi = build_initial_conditions(args_multi)
        model_multi = BaseThrusterModel(
            thrust=[2.0, 3.0],
            direction=[0.0, π],
            Δv=[20.0, 10.0],
            start_burn_time=[-1.0, -2.0],
            stop_burn_time=[-3.0, -4.0],
            Isp=[300.0, 300.0]
        )
        _clear_burn_plan!(p_multi, 1)
        _clear_burn_plan!(p_multi, 2)
        calcControlEffect!(model_multi, state_multi, p_multi, 100.0, 1)
        @test _burn_plan(p_multi, 1).valid
        @test !_burn_plan(p_multi, 2).valid
        @test model_multi.start_burn_time[1] == _burn_plan(p_multi, 1).start_burn_s
        @test model_multi.stop_burn_time[1] == _burn_plan(p_multi, 1).stop_burn_s
        @test model_multi.start_burn_time[2] == -2.0
        @test model_multi.stop_burn_time[2] == -4.0

        s1 = _burn_plan(p_multi, 1).start_burn_s
        e1 = _burn_plan(p_multi, 1).stop_burn_s
        calcControlEffect!(model_multi, state_multi, p_multi, 100.0, 2)
        @test _burn_plan(p_multi, 1).valid
        @test _burn_plan(p_multi, 1).start_burn_s == s1
        @test _burn_plan(p_multi, 1).stop_burn_s == e1
        @test _burn_plan(p_multi, 2).valid
        @test model_multi.start_burn_time[1] == s1
        @test model_multi.stop_burn_time[1] == e1
        @test model_multi.start_burn_time[2] == _burn_plan(p_multi, 2).start_burn_s
        @test model_multi.stop_burn_time[2] == _burn_plan(p_multi, 2).stop_burn_s

        model_oob = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        @test_nowarn calcControlEffect!(model_oob, state, p, 100.0, 2)
        @test model_oob.start_burn_time[1] == -1.0
        @test model_oob.stop_burn_time[1] == -1.0

        throw_planet = ThrowingOrbitPlanet(EARTH.Rp_e)
        args_throw = build_config(
            spacecraft=sc,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true,
            planet=throw_planet
        )
        p_throw = ODEParams{1}(args=args_throw)
        model_throw = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=123.0, stop_burn_time=124.0)
        withenv("SPACEAGORA_DEBUG_CONTROL" => "1", "SPACEAGORA_STRICT_CONTROL_EXCEPTIONS" => "0") do
            @test_logs (:warn, r"orbital-element conversion failed") calcControlEffect!(model_throw, state, p_throw, 100.0, 1)
        end
        @test model_throw.start_burn_time[1] == 123.0
        @test model_throw.stop_burn_time[1] == 124.0
        withenv("SPACEAGORA_STRICT_CONTROL_EXCEPTIONS" => "1") do
            @test_throws ErrorException calcControlEffect!(model_throw, state, p_throw, 100.0, 1)
        end

        trace_key_helper = SimulationModel.ControlHooks._maneuver_trace_key
        trace_bool_helper = SimulationModel.ControlHooks._trace_bool_enabled
        trace_path_helper = SimulationModel.ControlHooks._maneuver_trace_path
        safe_orbit_counter_helper = SimulationModel.ControlHooks._safe_orbit_counter
        trace_event_helper = SimulationModel.ControlHooks._trace_maneuver_event!

        @test trace_bool_helper("yes")
        @test trace_bool_helper("ON")
        @test !trace_bool_helper("0")
        @test trace_key_helper(model_throw, 1) isa Tuple{UInt64, Int64}

        p_no_orbit_counter = (shared_buffers=(debug_control=Ref(false),),)
        @test safe_orbit_counter_helper(p_no_orbit_counter, 1) == -1

        mktempdir() do tmp
            trace_csv = joinpath(tmp, "maneuver_trace.csv")
            withenv(
                "SPACEAGORA_TRACE_MANEUVERS" => "1",
                "SPACEAGORA_MANEUVER_TRACE_CSV" => trace_csv
            ) do
                @test trace_path_helper() == trace_csv

                model_trace = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
                state_trace = build_initial_conditions(args)
                p_trace = ODEParams{1}(args=args)

                _clear_burn_plan!(p_trace, 1)
                calcControlEffect!(model_trace, state_trace, p_trace, 100.0, 1)
                s_trace = _burn_plan(p_trace, 1).start_burn_s
                e_trace = _burn_plan(p_trace, 1).stop_burn_s
                @test s_trace < e_trace

                calcControlEffect!(model_trace, state_trace, p_trace, s_trace + 1e-6, 1)
                calcControlEffect!(model_trace, state_trace, p_trace, e_trace + 1.0, 1)
                trace_plan_after = _burn_plan(p_trace, 1)
                @test !trace_plan_after.valid ||
                      trace_plan_after.start_burn_s != s_trace ||
                      trace_plan_after.stop_burn_s != e_trace

                @test isfile(trace_csv)
                trace_text = read(trace_csv, String)
                @test occursin("event,t_s,spacecraft_idx", trace_text)
                @test occursin("schedule_set", trace_text)
                @test occursin("burn_start", trace_text)
                @test occursin("burn_end", trace_text)
                @test occursin("schedule_clear", trace_text)

                trace_event_helper("manual_event", model_trace, p_no_orbit_counter, 1, 0.0)
                trace_text_after_manual = read(trace_csv, String)
                @test occursin("manual_event", trace_text_after_manual)
            end
        end

        helper = SimulationModel.ControlHooks._control_effector_exception_fallback
        p_dummy = (shared_buffers=(debug_control=Ref(false),),)
        err_eff = ErrorException("control-effector-fallback")
        @test helper(p_dummy, 1, err_eff, stacktrace()) === nothing
        withenv("SPACEAGORA_DEBUG_CONTROL" => "1", "SPACEAGORA_STRICT_CONTROL_EXCEPTIONS" => "0") do
            @test_logs (:warn, r"orbital-element conversion failed") helper(p_dummy, 1, err_eff, stacktrace())
        end
        withenv("SPACEAGORA_STRICT_CONTROL_EXCEPTIONS" => "1") do
            @test_throws ErrorException helper(p_dummy, 1, err_eff, stacktrace())
        end

        burn_helper = SimulationModel.ControlHooks._constant_thrust_burn_duration_s
        @test burn_helper(100.0, 0.0, 2.0, 300.0) == 0.0
        @test isapprox(
            burn_helper(100.0, 20.0, 2.0, 300.0),
            100.0 * 300.0 * 9.80665 / 2.0 * (1.0 - exp(-20.0 / (300.0 * 9.80665)));
            atol=1e-12,
            rtol=0.0
        )
        @test isnan(burn_helper(0.0, 20.0, 2.0, 300.0))
        @test isnan(burn_helper(100.0, -1.0, 2.0, 300.0))
        @test isnan(burn_helper(100.0, 20.0, 0.0, 300.0))
        @test isnan(burn_helper(100.0, 20.0, 2.0, 0.0))

        @testset "typed maneuver and helper boundary branches" begin
            hooks = SimulationModel.ControlHooks
            guidance_command = hooks._guidance_maneuver_command
            burn_buffer = hooks._burn_plan_buffer
            active_plan = hooks._active_burn_plan
            set_plan! = hooks._set_burn_plan!
            clear_plan! = hooks._clear_burn_plan!
            commanded = hooks._commanded_maneuver
            model_window = hooks._model_burn_window
            effective_window = hooks._effective_burn_window
            effective_direction = hooks._effective_direction_rad
            effective_thrust_isp = hooks._effective_thrust_isp
            available_propellant = hooks._available_propellant_kg
            validated_plan = hooks._validated_burn_plan

            p_no_shared = (;)
            @test guidance_command(p_no_shared, 1) === nothing
            @test burn_buffer(p_no_shared) === nothing
            @test active_plan(p_no_shared, 1) === nothing
            @test set_plan!(p_no_shared, 1, PropulsiveBurnPlan(valid=true)) === nothing
            @test clear_plan!(p_no_shared, 1) === nothing

            p_helper = ODEParams{1}(args=args)
            model_helper = make_base_thruster_model(
                thrust=2.0,
                direction=0.25,
                Δv=3.0,
                start_burn_time=-1.0,
                stop_burn_time=-1.0,
                Isp=300.0
            )
            @test guidance_command(p_helper, 0) === nothing
            @test active_plan(p_helper, 0) === nothing
            @test set_plan!(p_helper, 0, PropulsiveBurnPlan(valid=true)) === nothing
            @test clear_plan!(p_helper, 0) === nothing
            @test model_window(model_helper, 0) === (NaN, NaN)
            @test isnan(effective_direction(model_helper, p_helper, 0))
            bad_thrust, bad_isp = effective_thrust_isp(model_helper, p_helper, 0)
            @test isnan(bad_thrust) && isnan(bad_isp)

            fallback_window = effective_window(model_helper, p_helper, 1)
            @test fallback_window == (-1.0, -1.0)
            plan_fallback = PropulsiveBurnPlan(
                valid=true,
                delta_v_mps=3.0,
                direction_rad=π,
                start_burn_s=10.0,
                stop_burn_s=20.0,
                thrust_n=5.0,
                isp_s=310.0
            )
            set_plan!(p_helper, 1, plan_fallback)
            @test active_plan(p_helper, 1).valid
            @test effective_window(model_helper, p_helper, 1) == (10.0, 20.0)
            @test isapprox(effective_direction(model_helper, p_helper, 1), π; atol=1e-12, rtol=0.0)
            @test effective_thrust_isp(model_helper, p_helper, 1) == (5.0, 310.0)
            clear_plan!(p_helper, 1)
            @test active_plan(p_helper, 1) === nothing

            p_helper.shared_buffers.maneuver_commands[1] = PropulsiveManeuverCommand(
                valid=true,
                delta_v_mps=-4.0,
                direction_rad=0.0,
                source_orbit=42
            )
            cmd = commanded(model_helper, p_helper, 1)
            @test cmd.delta_v_mps == 4.0
            @test cmd.direction_rad == π
            @test cmd.source_orbit == 42
            @test model_helper.Δv[1] == 4.0
            @test isapprox(model_helper.direction[1], π; atol=1e-12, rtol=0.0)

            model_nan = make_base_thruster_model(
                thrust=2.0,
                direction=0.5,
                Δv=NaN,
                start_burn_time=-1.0,
                stop_burn_time=-1.0,
                Isp=300.0
            )
            p_helper.shared_buffers.maneuver_commands[1] = PropulsiveManeuverCommand()
            nan_cmd = commanded(model_nan, p_helper, 1)
            @test isnan(nan_cmd.delta_v_mps)
            @test model_nan.direction[1] == 0.5
            @test commanded(model_nan, p_helper, 2) === nothing

            @test available_propellant(p_no_shared, 1, 100.0) === nothing
            @test available_propellant(p_helper, 2, 100.0) === nothing
            @test available_propellant(p_helper, 1, 100.0) === nothing

            maneuver_valid = (delta_v_mps=1.0, direction_rad=0.0, source_orbit=1)
            @test validated_plan(model_helper, p_helper, 1, -1.0, maneuver_valid) === nothing
            @test validated_plan(model_helper, p_helper, 1, 100.0, (delta_v_mps=0.0, direction_rad=0.0, source_orbit=1)) === nothing
            @test validated_plan(model_helper, p_helper, 1, 100.0, (delta_v_mps=1.0, direction_rad=NaN, source_orbit=1)) === nothing
            model_zero_isp = make_base_thruster_model(thrust=2.0, Δv=1.0, Isp=0.0)
            @test validated_plan(model_zero_isp, p_helper, 1, 100.0, maneuver_valid) === nothing
            model_zero_thrust_helper = make_base_thruster_model(thrust=0.0, Δv=1.0, Isp=300.0)
            @test validated_plan(model_zero_thrust_helper, p_helper, 1, 100.0, maneuver_valid) === nothing
            @test validated_plan(model_helper, p_helper, 1, 100.0, maneuver_valid).valid
        end
    end

    @testset "schmitt_trigger" begin
        @test schmitt_trigger(0.80, 0.75, 0.25) == 1.0
        @test schmitt_trigger(0.20, 0.75, 0.25) == 0.0
        @test schmitt_trigger(0.75, 0.75, 0.25) == 0.0
        @test schmitt_trigger(0.50, 0.75, 0.25) == 0.0

        # Sequential behavior check: current implementation is memoryless.
        state_high = schmitt_trigger(0.80, 0.75, 0.25)
        state_mid_after_high = schmitt_trigger(0.50, 0.75, 0.25)
        _ = schmitt_trigger(0.20, 0.75, 0.25)
        state_mid_after_low = schmitt_trigger(0.50, 0.75, 0.25)
        @test state_high == 1.0
        @test state_mid_after_high == 0.0
        @test state_mid_after_low == 0.0
        @test state_mid_after_high == state_mid_after_low
    end

    @testset "integrate_impulse!" begin
        link = Link{0}(root=true, attitude_control_rate=0.1)
        ω = 20.0

        thr_zero = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.0)
        impulse_zero = integrate_impulse!(link, thr_zero, 0.0, 0.0)
        @test impulse_zero == 0.0
        @test thr_zero.κ == 0.0

        thr_full = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.0)
        dt = link.attitude_control_rate
        expected_full = thr_full.max_thrust * (dt + (-1.0) / ω * (1 - exp(-ω * dt)))
        impulse_full = integrate_impulse!(link, thr_full, dt, 0.0)
        @test isapprox(impulse_full, expected_full; atol=1e-12, rtol=0.0)
        @test isapprox(thr_full.κ, 1 - exp(-ω * dt); atol=1e-12, rtol=0.0)

        thr_decay = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=1.0)
        expected_decay = thr_decay.max_thrust / ω * (1 - exp(-ω * dt))
        impulse_decay = integrate_impulse!(link, thr_decay, 0.0, 0.0)
        @test isapprox(impulse_decay, expected_decay; atol=1e-12, rtol=0.0)
        @test isapprox(thr_decay.κ, exp(-ω * dt); atol=1e-12, rtol=0.0)

        thr_small_ω = Thruster(max_thrust=10.0, cutoff_frequency=1e-12, κ=0.3)
        impulse_small_ω = integrate_impulse!(link, thr_small_ω, 0.05, 0.0)
        @test isfinite(impulse_small_ω)
        @test impulse_small_ω >= 0.0
        @test isfinite(thr_small_ω.κ)

        thr_ref_neg = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.4)
        thr_neg = deepcopy(thr_ref_neg)
        thr_zero_ref = deepcopy(thr_ref_neg)
        impulse_neg = integrate_impulse!(link, thr_neg, -0.01, 0.0)
        impulse_zero_ref = integrate_impulse!(link, thr_zero_ref, 0.0, 0.0)
        @test isapprox(impulse_neg, impulse_zero_ref; atol=1e-12, rtol=0.0)
        @test isapprox(thr_neg.κ, thr_zero_ref.κ; atol=1e-12, rtol=0.0)

        thr_ref_hi = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.4)
        thr_hi = deepcopy(thr_ref_hi)
        thr_dt_ref = deepcopy(thr_ref_hi)
        impulse_hi = integrate_impulse!(link, thr_hi, 10.0 * link.attitude_control_rate, 0.0)
        impulse_dt_ref = integrate_impulse!(link, thr_dt_ref, link.attitude_control_rate, 0.0)
        @test isapprox(impulse_hi, impulse_dt_ref; atol=1e-12, rtol=0.0)
        @test isapprox(thr_hi.κ, thr_dt_ref.κ; atol=1e-12, rtol=0.0)
    end

    @testset "thrust_calculation_schmitt_trigger!" begin
        mktempdir() do tmp
            cd(tmp) do
                link = Link{0}(root=true, attitude_control_rate=0.1)
                thr = Thruster(
                    max_thrust=10.0,
                    min_firing_time=0.05,
                    level_on=0.75,
                    level_off=0.25,
                    cutoff_frequency=100.0,
                    κ=0.0
                )
                debug_key = "SPACEAGORA_DEBUG_THRUSTER"
                old_debug = get(ENV, debug_key, nothing)
                try
                    if haskey(ENV, debug_key)
                        delete!(ENV, debug_key)
                    end
                    thrust_calculation_schmitt_trigger!(link, thr, 1.0, 0.0)
                    @test thr.thrust == 0.0
                    @test !isfile("thruster_debug.csv")

                    thrust_calculation_schmitt_trigger!(link, thr, 9.0, 0.1)
                    @test thr.thrust > 0.0
                    @test thr.thrust <= thr.max_thrust + 1e-9
                    @test !isfile("thruster_debug.csv")

                    ENV[debug_key] = "1"
                    thrust_calculation_schmitt_trigger!(link, thr, 9.0, 0.2)
                    @test isfile("thruster_debug.csv")
                finally
                    if old_debug === nothing
                        if haskey(ENV, debug_key)
                            delete!(ENV, debug_key)
                        end
                    else
                        ENV[debug_key] = old_debug
                    end
                end

                thr_zero_max = Thruster(max_thrust=0.0, cutoff_frequency=100.0, κ=0.0)
                @test_nowarn thrust_calculation_schmitt_trigger!(link, thr_zero_max, 1.0, 0.2)
                @test thr_zero_max.thrust == 0.0
            end
        end
    end

    @testset "update_thrusters!" begin
        mktempdir() do tmp
            cd(tmp) do
                link = Link{0}(root=true, attitude_control_rate=0.1)
                update_thrusters!(link, SVector{3, Float64}(0.0, 0.0, 0.0), 0.0)
                @test isempty(link.thrusters)
                @test size(link.J_thruster) == (3, 0)

                thr = Thruster(
                    max_thrust=1.0,
                    cutoff_frequency=100.0,
                    min_firing_time=0.0,
                    location=MVector{3, Float64}(0.0, 1.0, 0.0),
                    direction=MVector{3, Float64}(0.0, 0.0, 1.0)
                )
                push!(link.thrusters, thr)

                update_thrusters!(link, SVector{3, Float64}(-1.0, 0.0, 0.0), 0.1)
                @test link.thrusters[1].thrust >= 0.0

                link_full = Link{0}(root=true, attitude_control_rate=0.1)
                push!(link_full.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 1.0, 0.0), direction=MVector{3, Float64}(0.0, 0.0, 1.0)))
                push!(link_full.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 0.0, 1.0), direction=MVector{3, Float64}(1.0, 0.0, 0.0)))
                push!(link_full.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(1.0, 0.0, 0.0), direction=MVector{3, Float64}(0.0, 1.0, 0.0)))
                τ_req_full = SVector{3, Float64}(1.0, 2.0, 3.0)
                update_thrusters!(link_full, τ_req_full, 0.2)
                @test rank(link_full.J_thruster) == 3
                @test all(isfinite, link_full.J_thruster)
                @test all(thr_i -> isfinite(thr_i.thrust) && thr_i.thrust >= 0.0, link_full.thrusters)
                @test any(thr_i -> thr_i.thrust > 0.0, link_full.thrusters)
                thrust_full = [thr_i.thrust for thr_i in link_full.thrusters]
                τ_ach_full = link_full.J_thruster * thrust_full
                @test norm(τ_ach_full - τ_req_full) / norm(τ_req_full) < 1e-6

                link_singular = Link{0}(root=true, attitude_control_rate=0.1)
                push!(link_singular.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 1.0, 0.0), direction=MVector{3, Float64}(0.0, 0.0, 1.0)))
                push!(link_singular.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 2.0, 0.0), direction=MVector{3, Float64}(0.0, 0.0, 1.0)))
                τ_req_singular = SVector{3, Float64}(0.0, 1.0, 0.0)
                update_thrusters!(link_singular, τ_req_singular, 0.3)
                @test rank(link_singular.J_thruster) == 1
                @test all(isfinite, link_singular.J_thruster)
                @test all(thr_i -> isfinite(thr_i.thrust) && thr_i.thrust >= 0.0, link_singular.thrusters)
                thrust_singular = [thr_i.thrust for thr_i in link_singular.thrusters]
                τ_ach_singular = link_singular.J_thruster * thrust_singular
                τ_lsq_singular = link_singular.J_thruster * (pinv(link_singular.J_thruster) * τ_req_singular)
                @test norm(τ_ach_singular - τ_lsq_singular) < 1e-6

                link_degenerate = Link{0}(root=true, attitude_control_rate=0.1)
                push!(
                    link_degenerate.thrusters,
                    Thruster(
                        max_thrust=50.0,
                        cutoff_frequency=100.0,
                        min_firing_time=0.0,
                        location=MVector{3, Float64}(0.5, 0.0, 0.0),
                        direction=MVector{3, Float64}(0.0, 0.0, 0.0)
                    )
                )
                @test_nowarn update_thrusters!(link_degenerate, SVector{3, Float64}(0.5, 0.0, 0.0), 0.4)
                @test all(isfinite, link_degenerate.J_thruster)
                @test link_degenerate.J_thruster[:, 1] == zeros(3)
                @test isfinite(link_degenerate.thrusters[1].thrust)
                @test link_degenerate.thrusters[1].thrust >= 0.0
            end
        end
    end
end

@testset "Control Burn Energy Sign (End-to-End)" begin
    sc0 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    args0 = build_config(
        spacecraft=sc0,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(),
        control_rates=Float64[],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df0 = run_case_silent(args0)
    eps0 = specific_energy(df0, EARTH.μ)
    Δeps0 = last(eps0) - first(eps0)
    @test abs(Δeps0) < 2e3

    scp = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    thr_pro = TimedTangentialThrusterModel(800.0, +1.0, 100.0, 101.0)
    argsp = build_config(
        spacecraft=scp,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thr_pro,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    dfp = run_case_silent(argsp)
    epsp = specific_energy(dfp, EARTH.μ)
    Δepsp = last(epsp) - first(epsp)
    @test Δepsp > 5e3

    scr = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    thr_retro = TimedTangentialThrusterModel(800.0, -1.0, 100.0, 101.0)
    argsr = build_config(
        spacecraft=scr,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thr_retro,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    dfr = run_case_silent(argsr)
    epsr = specific_energy(dfr, EARTH.μ)
    Δepsr = last(epsr) - first(epsr)
    @test Δepsr < -5e3
end

@testset "Control Mass Flow Coupling (End-to-End)" begin
    sc_no_control = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    args_no_control = build_config(
        spacecraft=sc_no_control,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=180.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(),
        control_rates=Float64[],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df_no_control = run_case_silent(args_no_control)
    @test "sc1_mass" in names(df_no_control)
    mass_no_control = Vector{Float64}(df_no_control.sc1_mass)
    @test maximum(abs.(mass_no_control .- first(mass_no_control))) < 1e-9

    function run_burn_case(isp::Float64)
        sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
        thruster = make_base_thruster_model(
            thrust=600.0,
            direction=0.0,
            Δv=0.0,
            start_burn_time=0.0,
            stop_burn_time=80.0,
            Isp=isp
        )
        args = build_config(
            spacecraft=sc,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=180.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            control_effectors=(thruster,),
            control_rates=[1.0],
            keplerian=true,
            tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
        )
        return run_case_silent(args)
    end

    df_low_isp = run_burn_case(200.0)
    df_high_isp = run_burn_case(400.0)

    mass_low = Vector{Float64}(df_low_isp.sc1_mass)
    mass_high = Vector{Float64}(df_high_isp.sc1_mass)
    Δm_low = first(mass_low) - last(mass_low)
    Δm_high = first(mass_high) - last(mass_high)

    @test Δm_low > 5.0
    @test Δm_high > 2.0
    @test all(diff(mass_low) .<= 1e-7)
    @test all(diff(mass_high) .<= 1e-7)
    @test isapprox(Δm_low / Δm_high, 2.0; atol=0.08, rtol=0.0)
end

@testset "Control Callback Multi-Spacecraft Mapping" begin
    sc1 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
    sc2 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
    shared_thruster = BaseThrusterModel(
        thrust=[800.0, 800.0],
        direction=[0.0, π],
        Δv=[20.0, -20.0],
        start_burn_time=[-1.0, -1.0],
        stop_burn_time=[-1.0, -1.0],
        Isp=[300.0, 300.0]
    )
    args = build_config_multi(
        spacecraft=[sc1, sc2],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=1000.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(shared_thruster,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df = run_case_silent(args; isolate_state=false)

    for sat_idx in 1:2
        @test isfinite(shared_thruster.start_burn_time[sat_idx])
        @test isfinite(shared_thruster.stop_burn_time[sat_idx])
        @test shared_thruster.stop_burn_time[sat_idx] > shared_thruster.start_burn_time[sat_idx]
    end

    mass1 = Vector{Float64}(df.sc1_mass)
    mass2 = Vector{Float64}(df.sc2_mass)
    @test first(mass1) - last(mass1) > 0.1
    @test first(mass2) - last(mass2) > 0.1

    eps1 = 0.5 .* (df.sc1_vel_1.^2 .+ df.sc1_vel_2.^2 .+ df.sc1_vel_3.^2) .-
           EARTH.μ ./ sqrt.(df.sc1_pos_1.^2 .+ df.sc1_pos_2.^2 .+ df.sc1_pos_3.^2)
    eps2 = 0.5 .* (df.sc2_vel_1.^2 .+ df.sc2_vel_2.^2 .+ df.sc2_vel_3.^2) .-
           EARTH.μ ./ sqrt.(df.sc2_pos_1.^2 .+ df.sc2_pos_2.^2 .+ df.sc2_pos_3.^2)
    Δeps1 = last(eps1) - first(eps1)
    Δeps2 = last(eps2) - first(eps2)
    @test Δeps1 > 2e3
    @test Δeps2 < -2e3
end

include(joinpath(REPO_ROOT, "test", "probes", "coverage_parallel_telemetry_probes.jl"))
include(joinpath(REPO_ROOT, "test", "probes", "coverage_runtime_boundary_probes.jl"))
include(joinpath(REPO_ROOT, "test", "probes", "coverage_targeted_90_probes.jl"))
include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "e_edg_strategy_parity_tests.jl"))
include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "t_edg_strategy_parity_tests.jl"))
include(joinpath(REPO_ROOT, "test", "mission", "aerobraking_policy_selector_stub_tests.jl"))

@testset "Aqua Package Quality" begin
    if HAS_AQUA
        Aqua.test_all(
            SimulationModel;
            ambiguities=false,
            stale_deps=false,
            deps_compat=false,
            project_extras=false,
            piracies=false,
            persistent_tasks=false,
            undocumented_names=false
        )
    else
        @test_skip "Aqua is not available in this test environment"
    end
end
