@testset "Guidance Thruster Campaign Helpers" begin
    ensure_guidance_sandbox_loaded!()
    sandbox = GUIDANCE_SANDBOX

    guidance_model = sandbox.AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=[2, 3],
        maneuver_Δv=[-5.0, 20.0]
    )

    args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p = ODEParams(n_sats=1, args=args)
    u = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])

    p.orbit_counter[1] = 3
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test p.shared_buffers.maneuver_commands[1].valid == true
    @test p.shared_buffers.maneuver_commands[1].delta_v_mps == 20.0
    @test isapprox(p.shared_buffers.maneuver_commands[1].direction_rad, 0.0; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 2
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test p.shared_buffers.maneuver_commands[1].valid == true
    @test p.shared_buffers.maneuver_commands[1].delta_v_mps == 5.0
    @test isapprox(p.shared_buffers.maneuver_commands[1].direction_rad, π; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 1
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test p.shared_buffers.maneuver_commands[1].valid == true
    @test p.shared_buffers.maneuver_commands[1].delta_v_mps == 0.0
end

@testset "Odyssey Maneuver Schedule Bridge" begin
    maneuvers = odyssey_campaign_maneuvers(1:20)

    @test !(1 in maneuvers.maneuver_orbit_number)
    @test 7 in maneuvers.maneuver_orbit_number
    @test 14 in maneuvers.maneuver_orbit_number

    idx7 = findfirst(==(7), maneuvers.maneuver_orbit_number)
    idx14 = findfirst(==(14), maneuvers.maneuver_orbit_number)
    @test idx7 !== nothing
    @test idx14 !== nothing
    @test maneuvers.maneuver_Δv[idx7] > 0.0
    @test maneuvers.maneuver_Δv[idx14] < 0.0

    guidance_model = AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=maneuvers.maneuver_orbit_number,
        maneuver_Δv=maneuvers.maneuver_Δv
    )
    thruster = BaseThrusterModel(
        thrust=[4.0],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[-1.0],
        stop_burn_time=[-1.0],
        Isp=[300.0]
    )
    @test length(guidance_model.maneuver_orbit_number) == length(guidance_model.maneuver_Δv)
    @test length(thruster.thrust) == 1

    ensure_guidance_sandbox_loaded!()
    sandbox = GUIDANCE_SANDBOX
    sandbox_guidance = sandbox.AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=maneuvers.maneuver_orbit_number,
        maneuver_Δv=maneuvers.maneuver_Δv
    )

    args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p = ODEParams(n_sats=1, args=args)
    u = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])

    p.orbit_counter[1] = 7
    sandbox.calcGuidanceEffect!(sandbox_guidance, u, p, 0.0, 1)
    @test p.shared_buffers.maneuver_commands[1].delta_v_mps > 0.0
    @test isapprox(p.shared_buffers.maneuver_commands[1].direction_rad, 0.0; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 14
    sandbox.calcGuidanceEffect!(sandbox_guidance, u, p, 0.0, 1)
    @test p.shared_buffers.maneuver_commands[1].delta_v_mps > 0.0
    @test isapprox(p.shared_buffers.maneuver_commands[1].direction_rad, π; atol=1e-12, rtol=0.0)
end

@testset "Mission Maneuver Plan Tables" begin
    function probe_plan!(plan, passage)
        args = Dict{Symbol, Float64}(:delta_v => 0.0, :phi => 0.0)
        plan(nothing, 0.0, 0.0, passage, args)
        return args
    end

    odyssey_passages = [7, 14, 26, 30, 35, 47, 54, 69, 72, 80, 87, 110,
        128, 161, 179, 195, 211, 223, 239, 251, 263, 274, 287, 299, 311]
    odyssey_full = odyssey_campaign_maneuvers(odyssey_passages)
    @test odyssey_full.maneuver_orbit_number == odyssey_passages
    @test all(abs.(odyssey_full.maneuver_Δv) .> 0.0)
    @test any(>(0.0), odyssey_full.maneuver_Δv)
    @test any(<(0.0), odyssey_full.maneuver_Δv)

    true_beginning_passages = [7, 10, 12, 14, 16, 25, 32, 45, 49, 53, 66, 74,
        88, 91, 99, 105, 129, 146, 179, 197, 213, 229, 241, 257, 270, 282,
        295, 306, 318, 330]
    for passage in true_beginning_passages
        args = probe_plan!(Odyssey_firing_plan_true_beginning, passage)
        @test args[:delta_v] > 0.0
        @test args[:phi] in (0.0, deg2rad(180))
    end
    @test probe_plan!(Odyssey_firing_plan_true_beginning, -1)[:delta_v] == 0.0

    magellan_passages = [50, 147, 185, 212, 238, 297, 444, 508, 599]
    for passage in magellan_passages
        args = probe_plan!(Magellan_firing_plan, passage)
        @test args[:delta_v] > 0.0
        @test args[:phi] in (0.0, deg2rad(180))
    end
    @test probe_plan!(Magellan_firing_plan, 999)[:delta_v] == 0.0

    venus_express_passages = [6, 37, 42, 46]
    for passage in venus_express_passages
        args = probe_plan!(Venus_Express_firing_plan, passage)
        @test args[:delta_v] > 0.0
        @test args[:phi] in (0.0, deg2rad(180))
    end
    @test probe_plan!(Venus_Express_firing_plan, 999)[:delta_v] == 0.0

    earth_low = Dict{Symbol, Float64}(:delta_v => 0.0, :phi => 0.0)
    Earth_firing_plan(EARTH, EARTH.Rp_e + 500e3, EARTH.Rp_e + 100e3, 1, earth_low)
    @test earth_low[:delta_v] > 0.0
    @test isapprox(earth_low[:phi], deg2rad(180.0); atol=1e-12, rtol=0.0)
    earth_ok = Dict{Symbol, Float64}(:delta_v => -1.0, :phi => NaN)
    Earth_firing_plan(EARTH, EARTH.Rp_e + 500e3, EARTH.Rp_e + 130e3, 1, earth_ok)
    @test earth_ok[:delta_v] == 0.0
    @test earth_ok[:phi] == 0.0

    titan_raise = Dict{Symbol, Float64}(:delta_v => 0.0, :phi => 0.0)
    titan_firing_plan(nothing, 2575.5e3 + 2.0e6, 2575.5e3 + 545e3, 1, titan_raise)
    @test titan_raise[:delta_v] > 0.0
    @test isapprox(titan_raise[:phi], deg2rad(180); atol=1e-12, rtol=0.0)
    titan_ok = Dict{Symbol, Float64}(:delta_v => -1.0, :phi => NaN)
    titan_firing_plan(nothing, 2575.5e3 + 2.0e6, 2575.5e3 + 700e3, 1, titan_ok)
    @test titan_ok[:delta_v] == 0.0
    @test titan_ok[:phi] == 0.0

    @test odyssey_campaign_maneuvers(Int[]) == (maneuver_orbit_number=Int64[], maneuver_Δv=Float64[])
    @test_throws ArgumentError odyssey_campaign_maneuvers([1]; firing_plan=(planet, ra, rp, passage, args) -> (args[:delta_v] = Inf))
    @test_throws ArgumentError odyssey_campaign_maneuvers([1]; firing_plan=(planet, ra, rp, passage, args) -> (args[:delta_v] = 1.0; args[:phi] = π / 2))
end
