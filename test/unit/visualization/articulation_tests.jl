using Test
using SpaceAGORA
using StaticArrays
using Arrow
using DataFrames
using JSON
using Base64

import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const SV = SM.SceneVisualization

# Aerobraking energy-depletion control on a Mars Odyssey-like orbit
# (examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl), started 930 s
# before periapsis so one drag pass fits in a 2400 s run. The control model
# articulates panel links 2 and 3 through `rotate_link`, which is exactly the
# in-place `Link` mutation the `link_pose` save field exists to capture.
function _aerobraking_control_config(; results_directory::String, mission_time::Float64=2400.0, data_rate::Float64=2.0)
    planet = make_no_gram_planet(:mars)
    ra = 28_559.615e3
    rp = planet.Rp_e + 77e3
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.2, 2.6, 1.7),
        panel_dims=(0.01, 3.76 / 2.0, 1.93),
        bus_mass=391.0,
        panel_mass_each=10.0,
        panel_offset_y=2.6 / 2.0 + 3.76 / 4.0,
        ic=SM.InitialCondition(ra=ra, rp=rp, i=93.6, ω=109.7454, Ω=28.1517, ν=300.0),
        reflection_coefficient=0.9,
        prop_mass=50.0,
        id=101
    )
    config = SM.AerobrakingEnergyDepletionConfig(
        guidance_modes=(:targeting, :max_energy_depletion),
        max_energy_submodes=(:heat_rate, :heat_load),
        heat_load_switch_solver=:tpbvp_integration,
        controlled_panel_links=(2, 3),
        target_apoapsis_radius_m=26_750e3,
        max_alpha_rad=pi / 2,
        min_alpha_rad=1e-4,
        heat_rate_limit_w_cm2=0.15,
        heat_load_limit_j_cm2=30.0,
        structural_load_limit_pa=0.5
    )
    state = SM.AerobrakingEnergyDepletionState(num_sats=1)
    for link in spacecraft.links
        link.α = config.min_alpha_rad
    end
    base = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=mission_time,
        initial_time=SM.InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredGravityModel(), SM.AerodynamicCoefficientfM()),
        density_model=SM.ExponentialAtmosphereModel(planet),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=false,
        keplerian=true,
        EI_km=160.0,
        verbose=false,
        results=true,
        results_directory=results_directory
    )
    return SM.SimulationConfiguration(
        file_paths=base.file_paths,
        simulation_settings=base.simulation_settings,
        mission_configuration=SM.MissionConfiguration(
            mission_type=base.mission_configuration.mission_type,
            keplerian=true,
            number_of_orbits=1,
            mission_time=mission_time,
            orientation_sim=false,
            num_steps_to_save=1000,
            data_rate=data_rate
        ),
        environment_model=base.environment_model,
        dynamics_model=base.dynamics_model,
        guidance_model=SM.GuidanceModel(guidance_effectors=(SM.AerobrakingEnergyDepletionGuidanceModel(config, state),), guidance_rates=[3.0]),
        navigation_model=base.navigation_model,
        control_model=SM.ControlModel(control_effectors=(SM.AerobrakingEnergyDepletionControlModel(config, state),), control_rates=[0.1]),
        initial_time=base.initial_time,
        integration_tolerances=base.integration_tolerances
    )
end

_decode_f32(b64) = collect(reinterpret(Float32, base64decode(b64)))

@testset "Articulation: aerobraking panel control recorded as link poses" begin
    dir = mktempdir()
    args = _aerobraking_control_config(results_directory=dir)
    run_simulation(args; visualization=true)

    feather = joinpath(dir, "simulation_results.feather")
    @test isfile(feather)
    @test isfile(joinpath(dir, "simulation_results_scene.json"))
    @test isfile(joinpath(dir, "simulation_results_viewer.html"))
    df = DataFrame(Arrow.Table(feather))
    @test nrow(df) > 100
    # The pass happens: altitude dips below the 160 km entry interface.
    @test minimum(df.sc1_altitude) < 160e3
    for k in 1:14
        @test "sc1_link_pose_$(k)" in names(df)
    end

    # Layout: [rx ry rz qx qy qz qw] for link 2 then link 3, in meters and scalar-last.
    scene = read_visualization_scene(joinpath(dir, "simulation_results_scene.json"))
    @test length(scene.spacecraft[1].links) == 3
    @test df[1, "sc1_link_pose_2"] ≈ scene.spacecraft[1].links[2].r_m[2]
    @test df[1, "sc1_link_pose_9"] ≈ scene.spacecraft[1].links[3].r_m[2]

    # Panel positions stay fixed on the bus ...
    for k in (1, 2, 3, 8, 9, 10)
        @test all(df[!, "sc1_link_pose_$(k)"] .== df[1, "sc1_link_pose_$(k)"])
    end
    # ... while both panels rotate about their hinge (body y) through the pass:
    # qy and qw change, qx and qz stay zero, and every quaternion stays unit.
    for (qx, qy, qz, qw) in ((4, 5, 6, 7), (11, 12, 13, 14))
        @test all(df[!, "sc1_link_pose_$(qx)"] .== 0.0)
        @test all(df[!, "sc1_link_pose_$(qz)"] .== 0.0)
        @test count(!=(0.0), diff(df[!, "sc1_link_pose_$(qy)"])) > 10
        @test maximum(df[!, "sc1_link_pose_$(qy)"]) > 0.5
        @test minimum(df[!, "sc1_link_pose_$(qw)"]) < 0.8
        norms = sqrt.(df[!, "sc1_link_pose_$(qy)"] .^ 2 .+ df[!, "sc1_link_pose_$(qw)"] .^ 2)
        @test all(abs.(norms .- 1.0) .< 1e-9)
        # The rotation is largest near periapsis and relaxed at the start.
        i_peri = argmin(df.sc1_altitude)
        @test abs(df[i_peri, "sc1_link_pose_$(qy)"]) > abs(df[1, "sc1_link_pose_$(qy)"])
    end

    # The viewer payload carries the same articulation.
    html = read(joinpath(dir, "simulation_results_viewer.html"), String)
    start = findfirst("window.SPACEAGORA_VIEWER = ", html)
    stop = findnext(";\n</script>", html, last(start))
    payload = JSON.parse(html[last(start)+1:first(stop)-1])
    lp = payload["frames"]["link_pose"]
    @test lp["counts"] == [2] && lp["total"] == 14
    data = _decode_f32(lp["data"])
    qy_link2 = data[5:14:end]
    @test length(qy_link2) == payload["frames"]["count"]
    @test maximum(qy_link2) > 0.5f0 && minimum(qy_link2) < 0.1f0
    @test payload["scene"]["spacecraft"][1]["links"][2]["name"] == "link1"

    # Phase 6: density is saved through the pass and the page carries the
    # diagnostics the trail is coloured by, plus the atmosphere description.
    @test "sc1_density" in names(df)
    @test maximum(df.sc1_density) > 1e-9
    @test df.sc1_density[argmin(df.sc1_altitude)] > df.sc1_density[1]
    for key in ("density_kg_m3", "heat_rate_w_m2", "drag_n", "wind_ms")
        @test payload["frames"][key] !== nothing
    end
    dens = _decode_f32(payload["frames"]["density_kg_m3"])
    @test length(dens) == payload["frames"]["count"] && maximum(dens) > 0.0f0
    @test length(_decode_f32(payload["frames"]["wind_ms"])) == 3 * payload["frames"]["count"]
    atmo = payload["scene"]["atmosphere"]
    @test atmo["model"] == "ExponentialAtmosphereModel"
    @test atmo["ei_altitude_m"] == 160e3
    @test length(atmo["profile"]["altitude_m"]) == SV.ATMOSPHERE_PROFILE_POINTS
    @test atmo["map"] === nothing

    # Explicit save_fields (the pattern the example scripts use) still get the
    # link poses appended when the flag is on.
    dir2 = mktempdir()
    args2 = _aerobraking_control_config(results_directory=dir2, mission_time=600.0)
    custom = vcat(SM.default_save_fields(args2), [SM.SaveField(:probe, (u, t, integrator) -> 1.0)])
    @test !any(f -> f.name === :link_pose, custom)
    run_simulation(args2; save_fields=custom, visualization=true)
    df2 = DataFrame(Arrow.Table(joinpath(dir2, "simulation_results.feather")))
    @test "probe" in names(df2)
    @test "sc1_link_pose_14" in names(df2)
    # A user who already asked for link_pose does not get it twice.
    args3 = with_visualization_scene(_aerobraking_control_config(results_directory=mktempdir(), mission_time=300.0), true)
    doubled = SM.default_save_fields(args3)
    @test count(f -> f.name === :link_pose, doubled) == 1
    run_simulation(args3; save_fields=doubled)
    df3 = DataFrame(Arrow.Table(joinpath(args3.simulation_settings.results_directory, "simulation_results.feather")))
    @test count(startswith("sc1_link_pose_"), names(df3)) == 14
end
