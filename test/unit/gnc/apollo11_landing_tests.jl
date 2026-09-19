module Apollo11LandingTests
# The Apollo 11 landing demo's vehicle, initial state and descent models,
# checked without propagating the descent: the seventeen-thruster four-quad
# geometry the control effector allocates over, the PDI state on the site's
# explicit datum, and the guidance, control and plume models sharing that datum
# on a synthetic flat site. SPICE kernels are needed only for the PDI state.
using Test, LinearAlgebra, StaticArrays, JSON, SpaceAGORA
const SM = SpaceAGORA.SimulationModel
const CH = SM.ControlHooks
module LandingDemo
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "viewer_demos", "apollo11_landing.jl"))
end
using .LandingDemo: pdi_state, lunar_module, descent_models, phase_code, PDI_ALTITUDE_M, PDI_UPRANGE_M, PDI_SPEED_MPS,
    DPS_THRUST_N, RCS_THRUST_N, RCS_QUAD_ARM_M, RCS_QUAD_Z_M, LM_MASS_KG, LM_PROPELLANT_KG, SPICE_PATH, PDI_UTC

const SITE_LAT = 0.67416
const SITE_LON = 23.47314
const SITE_HEIGHT_M = -1925.5
const RADIUS_M = 1_737_400.0

"A synthetic flat site bundle in the documented layout: a fine grid at the site inside a corridor-wide grid."
function synthetic_flat_site(dir)
    function grid(name, half_deg, n)
        header = Dict("rows" => n, "cols" => n, "lat_min" => SITE_LAT - half_deg, "lat_max" => SITE_LAT + half_deg,
            "lon_min" => SITE_LON - half_deg, "lon_max" => SITE_LON + half_deg, "reference_radius_m" => RADIUS_M,
            "source" => "synthetic flat test fixture", "units" => "m")
        open(joinpath(dir, name * ".json"), "w") do io; JSON.print(io, header); end
        bits = reinterpret(UInt32, Float32(SITE_HEIGHT_M))
        cell = UInt8[(bits >> shift) & 0xff for shift in (0, 8, 16, 24)]   # little-endian, as the loader documents
        write(joinpath(dir, name * ".f32"), repeat(cell, n * n))
        return Dict("name" => name, "reference_radius_m" => RADIUS_M)
    end
    entries = [grid("dem_site", 0.05, 9), grid("dem_corridor", 12.0, 25)]
    path = joinpath(dir, "site.json")
    open(path, "w") do io
        JSON.print(io, Dict("site" => Dict("lat_deg" => SITE_LAT, "lon_deg" => SITE_LON, "name" => "synthetic flat site"), "dem" => entries))
    end
    return path
end

@testset "Apollo 11 landing demo: requested calendar epoch" begin
    @test PDI_UTC == "1969-07-20T20:05:05"
    requested = SM.InitialTime(year=1969, month=7, day=20, hour=20, minute=5, second=5.0)
    @test SM.EphemeridesModels._initial_time_utc_string(requested) == "1969-07-20T20:05:05.000"
end

@testset "Apollo 11 landing demo: vehicle geometry" begin
    q = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
    sc = lunar_module(SM.CartesianInitialCondition(SVector(1.8e6, 0.0, 0.0), SVector(0.0, 1.7e3, 0.0); q=q))
    thrusters = sc.root.thrusters
    @test length(thrusters) == 17
    @test count(t -> t.max_thrust == DPS_THRUST_N, thrusters) == 1
    @test count(t -> t.max_thrust == RCS_THRUST_N, thrusters) == 16
    dps = only(filter(t -> t.max_thrust == DPS_THRUST_N, thrusters))
    @test SVector{3, Float64}(dps.direction) == SVector(0.0, 0.0, 1.0)        # display exhaust direction; the control effector applies thrust along body -z
    quads = unique(SVector{3, Float64}(t.location) for t in thrusters if t.max_thrust == RCS_THRUST_N)
    @test length(quads) == 4
    @test all(q -> abs(q[1]) == RCS_QUAD_ARM_M && abs(q[2]) == RCS_QUAD_ARM_M && q[3] == RCS_QUAD_Z_M, quads)
    @test Set((sign(q[1]), sign(q[2])) for q in quads) == Set([(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)])
    for q in quads
        jets = [SVector{3, Float64}(t.direction) for t in thrusters if t.max_thrust == RCS_THRUST_N && SVector{3, Float64}(t.location) == q]
        @test length(jets) == 4
        @test SVector(0.0, 0.0, 1.0) in jets && SVector(0.0, 0.0, -1.0) in jets           # one up, one down
        @test SVector(sign(q[1]), 0.0, 0.0) in jets && SVector(0.0, sign(q[2]), 0.0) in jets   # two lateral, outward from the quad
    end
    @test all(t -> t.thrust == 0.0, thrusters)   # the control effector applies the forces; the engine's thruster state stays idle
    @test sc.prop_mass == LM_PROPELLANT_KG && sc.root.m == LM_MASS_KG
    # the control effector's allocation over this layout: the DPS is the engine and the sixteen jets reach all three axes
    layout = CH.descent_thruster_layout(sc)
    @test layout !== nothing
    @test layout.engine_max_thrust_n == DPS_THRUST_N
    @test length(layout.jets) == 16
    @test rank(layout.torque_arms_nm) == 3
    @test size(layout.torque_to_levels) == (16, 3)
    levels = zeros(17)
    CH.descent_thruster_levels!(levels, layout, 0.5 * DPS_THRUST_N, SVector(0.0, 0.0, 0.0))
    @test levels[layout.engine] ≈ 0.5
    @test all(iszero, levels[layout.jets])
    CH.descent_thruster_levels!(levels, layout, 0.5 * DPS_THRUST_N, SVector(300.0, -200.0, 100.0))
    @test all(0.0 .<= levels .<= 1.0)
    @test any(>(0.0), levels[layout.jets])
end

@testset "Apollo 11 landing demo: descent models on the site's explicit datum" begin
    mktempdir() do dir
        site_json = synthetic_flat_site(dir)
        terrain, site = load_site_terrain(site_json)
        @test site.reference_radius_m == RADIUS_M
        @test site.height_m ≈ SITE_HEIGHT_M
        @test terrain_height(terrain, SITE_LAT + 5.0, SITE_LON - 5.0) ≈ SITE_HEIGHT_M   # the corridor grid, not the fallback
        guidance, control, plume, state = descent_models(site, terrain)
        @test guidance.config.reference_radius_m == RADIUS_M
        @test guidance.config.site_height_m ≈ SITE_HEIGHT_M
        @test control.state === state && guidance.state === state
        @test plume.control === control
        @test SM.touchdown_spec(control, 1) !== nothing
        @test SM.touchdown_spec(control, 2) === nothing
        @test state.phase[1] === :braking && phase_code(state.phase[1]) == 1.0
        @test phase_code(:landed) == 4.0
        frame = SM.descent_site_frame(guidance.config, terrain, RADIUS_M)
        @test norm(frame.origin_p) ≈ RADIUS_M + SITE_HEIGHT_M
        @test frame.height_m ≈ SITE_HEIGHT_M
        onset = SM.plume_erosion_onset_height(plume.config, 11_500.0)
        @test isfinite(onset) && onset > 0.0
    end
end

const SPICE_READY = isdir(SPICE_PATH) && isfile(joinpath(SPICE_PATH, "lsk", "naif0012.tls")) &&
    isfile(joinpath(SPICE_PATH, "spk", "satellites", "SPICELunaCurrentKernel.bpc")) &&
    isfile(joinpath(SPICE_PATH, "tf", "SPICELunaFrameKernel.tf"))
if SPICE_READY
    @testset "Apollo 11 landing demo: PDI state" begin
        mktempdir() do dir
            terrain, site = load_site_terrain(synthetic_flat_site(dir))
            planet = SM.Moon("", SPICE_PATH)
            et0 = LandingDemo.et_of(LandingDemo.initial_time_of(LandingDemo.et_of(PDI_UTC)))
            # Pin the calendar request, not a self-consistent but shifted round trip.
            @test LandingDemo.utc_of(et0) == "1969-07-20T20:05:05.000"
            requested_et = lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
                LandingDemo.str2et("1969-07-20T20:05:05.000")
            end
            @test isapprox(et0, requested_et; atol=1e-6, rtol=0.0)
            r_i, v_i, q_pdi = pdi_state(site, planet, et0)
            @test norm(r_i) ≈ RADIUS_M + site.height_m + PDI_ALTITUDE_M rtol=1e-9
            @test norm(q_pdi) ≈ 1.0
            # back in the rotating body frame: horizontal at the local speed, and the uprange arc to the site
            l_pi = SM.planet_frame_lpi(planet, et0, SM.SpiceEphemeridesModel())
            r_p = SVector{3, Float64}(l_pi * r_i)
            v_p = SVector{3, Float64}(l_pi * v_i) - cross(SVector{3, Float64}(planet.ω), r_p)
            @test abs(dot(v_p, normalize(r_p))) < 1e-6 * PDI_SPEED_MPS
            @test norm(v_p) ≈ PDI_SPEED_MPS rtol=1e-9
            φ = deg2rad(site.lat_deg); λ = deg2rad(site.lon_deg)
            up_site = SVector(cos(φ) * cos(λ), cos(φ) * sin(λ), sin(φ))
            arc = acos(clamp(dot(normalize(r_p), up_site), -1.0, 1.0)) * RADIUS_M
            @test arc ≈ PDI_UPRANGE_M rtol=1e-9
            @test dot(v_p, up_site - dot(up_site, normalize(r_p)) * normalize(r_p)) > 0   # moving toward the site
            # engine retrograde at PDI: body -z along the thrust, opposite the inertial velocity.
            # rot(q) is the passive inertial-to-body matrix, so its transpose maps body axes out.
            R_bi = SM.rot(q_pdi)'
            @test dot(R_bi * SVector(0.0, 0.0, -1.0), -normalize(v_i)) ≈ 1.0 atol=1e-9
            @test dot(R_bi * SVector(1.0, 0.0, 0.0), normalize(r_i)) > 0.9       # windows up
        end
    end
else
    @info "SPICE kernels absent; set SPACEAGORA_SPICE_PATH to run the Apollo 11 PDI-state test" SPICE_PATH
    @test_skip SPICE_READY
end
end # module
