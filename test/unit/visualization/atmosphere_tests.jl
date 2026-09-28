using Test
using SpaceAGORA
using StaticArrays
using JSON

import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const SV = SM.SceneVisualization

function _mars_config(density_model; results_directory::String=mktempdir(), EI_km::Float64=160.0)
    planet = make_no_gram_planet(:mars)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.5),
        panel_dims=(0.01, 2.5, 1.0),
        bus_mass=500.0,
        panel_mass_each=10.0,
        panel_offset_y=2.3,
        ic=SM.InitialCondition(ra=4_500.0e3, rp=3_800.0e3, i=30.0, ω=0.0, Ω=45.0, ν=0.0),
        prop_mass=0.0,
        id=1
    )
    return make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=300.0,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredJ2GravityModel(),),
        density_model=density_model,
        ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=false,
        keplerian=true,
        EI_km=EI_km,
        verbose=false,
        results=true,
        results_directory=results_directory
    )
end

# A stand-in model that varies with longitude and needs the 7-argument form, like GRAM.
struct _LonVaryingModel <: SM.AbstractDensityModel end
SM.EnvironmentModels.getDensity(::_LonVaryingModel, h::Float64, lat::Float64, lon::Float64, t::Float64, wind::Bool, p) =
    (1e-9 * exp(-h / 20e3) * (1.5 + cos(lon)), 200.0, SVector{3, Float64}(0.0, 0.0, 0.0))

# A custom model can mutate on either density API. Merely supplying its
# parameters to an exporter must not opt into those calls.
struct _MutatingVisualizationDensityModel <: SM.AbstractDensityModel
    calls::Base.RefValue{Int}
end
function SM.EnvironmentModels.getDensity(model::_MutatingVisualizationDensityModel,
    h::Float64, lat::Float64, lon::Float64, t::Float64, wind::Bool)
    model.calls[] += 1
    return (1e-9, 200.0, SVector{3, Float64}(0.0, 0.0, 0.0))
end
SM.EnvironmentModels.getDensity(model::_MutatingVisualizationDensityModel,
    h::Float64, lat::Float64, lon::Float64, t::Float64, wind::Bool, p) =
    SM.EnvironmentModels.getDensity(model, h, lat, lon, t, wind)

@testset "Atmosphere spec" begin
    planet = make_no_gram_planet(:mars)

    @testset "no atmosphere means no spec and no density field" begin
        args = with_visualization_scene(_mars_config(SM.NoAtmosphereModel()), true)
        @test atmosphere_spec(args) === nothing
        scene = build_visualization_scene(args; rotation_max_samples=4)
        @test scene.atmosphere === nothing
        @test SV.scene_dict(scene)["atmosphere"] === nothing
        @test read_visualization_scene(write_visualization_scene(joinpath(mktempdir(), "s.json"), scene)) == scene
        @test !any(f -> f.name === :density, SM.default_save_fields(args))
        @test !any(f -> f.name === :link_pose, SM.SimulationCallbacks.visualization_save_fields(_mars_config(SM.NoAtmosphereModel())))
    end

    @testset "exponential model: profile, no map, density field" begin
        model = SM.ExponentialAtmosphereModel(planet)
        args = with_visualization_scene(_mars_config(model; EI_km=160.0), true)
        spec = atmosphere_spec(args; profile_points=16)
        @test spec !== nothing
        @test spec.model == "ExponentialAtmosphereModel"
        @test spec.ei_altitude_m == 160e3
        @test length(spec.profile_altitude_m) == 16
        @test spec.profile_altitude_m[1] == 0.0 && spec.profile_altitude_m[end] == 1.5 * 160e3
        @test all(diff(spec.profile_density_kg_m3) .<= 0.0)      # monotone decreasing with altitude
        @test spec.profile_density_kg_m3[1] > spec.profile_density_kg_m3[end] > 0.0
        # Matches the model directly.
        ρ_direct = SM.getDensity(model, 100e3, 0.0, 0.0, 0.0, false)[1]
        k = findfirst(==(100e3), spec.profile_altitude_m)
        @test k === nothing || spec.profile_density_kg_m3[k] ≈ ρ_direct
        @test isempty(spec.map_density_kg_m3) && isempty(spec.map_lat_deg)   # altitude-only model: no map
        @test spec.map_altitude_m == 0.6 * 160e3

        fields = SM.default_save_fields(args)
        names_on = Symbol[f.name for f in fields]
        @test :density in names_on && :link_pose in names_on
        @test length(SM.SimulationCallbacks.visualization_save_fields(args)) == 2
        @test SM.SimulationCallbacks.visualization_save_fields(args)[end].column_prefix == "density"

        scene = build_visualization_scene(args; rotation_max_samples=4)
        @test scene.atmosphere == spec || length(scene.atmosphere.profile_altitude_m) == SV.ATMOSPHERE_PROFILE_POINTS
        d = SV.scene_dict(scene)
        @test d["atmosphere"]["model"] == "ExponentialAtmosphereModel"
        @test d["atmosphere"]["map"] === nothing
        @test length(d["atmosphere"]["profile"]["altitude_m"]) == SV.ATMOSPHERE_PROFILE_POINTS
        path = write_visualization_scene(joinpath(mktempdir(), "scene.json"), scene)
        back = read_visualization_scene(path)
        @test back == scene
        @test back.atmosphere.profile_density_kg_m3 == scene.atmosphere.profile_density_kg_m3
    end

    @testset "standalone live-model sampling requires explicit opt-in" begin
        args = with_visualization_scene(_mars_config(_LonVaryingModel(); EI_km=100.0), true)
        # Without integrator parameters the 7-argument-only model cannot be sampled: empty profile, no map.
        bare = atmosphere_spec(args; profile_points=8)
        @test bare !== nothing && isempty(bare.profile_altitude_m) && isempty(bare.map_density_kg_m3)
        # Integrator parameters alone never permit automatic live sampling.
        guarded = atmosphere_spec(args; density_params=(dummy=true,), profile_points=8)
        @test isempty(guarded.profile_altitude_m) && isempty(guarded.map_density_kg_m3)
        # Standalone opt-in enables the requested profile and map.
        spec = atmosphere_spec(args; density_params=(dummy=true,), sample_model=true, profile_points=8, map_step_deg=30.0)
        @test length(spec.profile_density_kg_m3) == 8
        @test length(spec.map_lat_deg) == 6 && length(spec.map_lon_deg) == 12
        @test length(spec.map_density_kg_m3) == 72
        @test spec.map_altitude_m == 60e3
        @test spec.map_lat_deg[1] == -75.0 && spec.map_lon_deg[1] == -165.0
        # Row-major by latitude; density peaks near longitude 0 on every row.
        row = spec.map_density_kg_m3[1:12]
        @test argmax(row) in (6, 7)
        # grid centers sit at ±15° and ±165°, so the extremes are 1.5 ± cos(15°)
        @test maximum(row) / minimum(row) ≈ (1.5 + cosd(15)) / (1.5 - cosd(15)) rtol=1e-6
        d = SV._scene_dict(spec)
        @test d["map"]["altitude_m"] == 60e3
        @test length(d["map"]["density_kg_m3"]) == length(d["map"]["lat_deg"]) * length(d["map"]["lon_deg"])
    end
    @testset "scene export never samples a stateful custom atmosphere" begin
        mktempdir() do dir
            calls = Ref(0)
            args = with_visualization_scene(
                _mars_config(_MutatingVisualizationDensityModel(calls); results_directory=dir), true,
            )
            params = (args=args,)
            spec = atmosphere_spec(args; density_params=params)
            @test calls[] == 0
            @test spec.ei_altitude_m == 160e3
            @test isempty(spec.profile_altitude_m) && isempty(spec.map_density_kg_m3)
            scene = build_visualization_scene(args; density_params=params, rotation_max_samples=4)
            @test calls[] == 0
            @test scene.atmosphere == spec
            path = SV.write_visualization_scene!(args; density_params=params)
            @test isfile(path)
            @test read_visualization_scene(path).atmosphere == spec
            @test calls[] == 0
            # Trajectory density remains available through the existing buffer
            # snapshot; an empty display profile does not remove the save field.
            @test any(field -> field.name === :density, SM.default_save_fields(args))
            @test calls[] == 0
            sampled = atmosphere_spec(args; density_params=params, sample_model=true,
                profile_points=4, map_step_deg=90.0)
            @test calls[] > 0
            @test length(sampled.profile_density_kg_m3) == 4
            @test length(sampled.map_density_kg_m3) == 8
        end
    end

end
