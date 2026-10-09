using SpaceAGORA, GRAMSuite
using Test, Serialization, SHA, StaticArrays, LinearAlgebra, Libdl, TOML

const EM = SpaceAGORA.SimulationModel.EnvironmentModels
const CB = SpaceAGORA.SimulationModel.SimulationCallbacks
const EN = SpaceAGORA.SimulationEngine
const ES = SpaceAGORA.SimulationModel.EffectorSampling

# Tiny analytic fields on a spherical reference body. These are interface
# fixtures, not Mars data or a scientific validation of the dependency's rule.
function near_surface_payload(; stored_winds::Bool)
    lower_levels = vcat(collect(-5.0:1.0:10.0), collect(15.0:5.0:75.0))
    levels = stored_winds ? vcat(lower_levels, [80.0323, 85.0323]) : lower_levels
    nl, ni, nj = length(levels), 3, 4
    payload = Dict{String,Any}(
        "format" => stored_winds ? "spaceagora_mars_near_surface_v2" : "spaceagora_mars_near_surface_scalars_v1",
        "planet" => "mars", "status" => "synthetic interface fixture", "radii_km" => (3390.0, 3390.0),
        "levels_km" => levels,
        "terrain" => Dict("lat0_deg" => -90.0, "lon0_deg" => 0.0, "step_deg" => 90.0,
            "surface_height_km" => fill(0.5, ni, nj), "areoid_radius_km" => fill(3390.0, ni, nj)),
        "lattice" => Dict("lat0_deg" => -90.0, "step_deg" => 90.0, "nlat" => ni, "nlon" => nj),
        "level_T_K" => fill(200.0, nl, ni, nj), "level_R" => fill(192.0, nl, ni, nj),
        "level_lnp" => [log(700.0) - z / 11.0 for z in levels, i in 1:ni, j in 1:nj],
        "level_source" => ones(UInt8, nl, ni, nj),
        "surface_T30_K" => fill(200.0, ni, nj), "surface_T5_K" => fill(200.0, ni, nj),
        # Queries are above the first level, so no surface-fit rows are needed.
        "q_models" => Dict("band" => Int[], "cell" => Int[], "L" => Int[], "order" => Int[],
            "phic_center" => Float64[], "lam_center" => Float64[], "coef" => zeros(0, 6)),
        "support" => Dict("zs_refuse_km" => 9.0, "phic_max_deg" => 85.0,
            "min_clearance_km" => 0.005, "top_areoid_km" => stored_winds ? 81.0 : 75.0))
    stored_winds || return payload
    for (name, dims, value) in (
        ("wind_level_U", (29, 2, nj), 10.0), ("wind_level_V", (29, 2, nj), -5.0),
        ("wind_mtgcm_U", (2, 2, nj), 10.0), ("wind_mtgcm_V", (2, 2, nj), -5.0),
        ("wind_surface_U", (2, 2, nj + 1), 10.0), ("wind_surface_V", (2, 2, nj), -5.0),
        ("sound_offset_level", (nl, ni, nj), 0.0), ("sound_offset_surface", (2, ni, nj), 0.0))
        payload[name] = fill(value, dims)
    end
    for name in ("wind_level_U", "wind_level_V", "wind_mtgcm_U", "wind_mtgcm_V", "sound_offset_level")
        payload[name * "_source"] = ones(UInt8, size(payload[name]))
    end
    temperatures = collect(50.0:50.0:350.0)
    payload["wind_meta"] = Dict("levels_km" => lower_levels, "longitudes_deg" => [0.0,90.0,180.0,270.0],
        "s_knots_deg" => [-90.0,0.0,90.0], "u_knots_deg" => [-90.0,90.0], "u_sides" => [0,0],
        "v_knots_deg" => [-90.0,90.0], "mtgcm_rows_deg" => [-90.0,90.0], "ho_km" => 0.0323, "solar_offset_h" => 12.0)
    payload["sound_reference"] = Dict("K" => fill(3.5, 56), "undetermined" => Int[],
        "temperature_levels_K" => temperatures, "pressure_levels_Pa" => [1e-2,1e-1,1.0,10.0,100.0,1000.0,1e4,1e5])
    payload["sound_composition"] = Dict("knee_kgkmol" => 43.5, "slope_per_kgkmol" => zeros(7),
        "slope_levels_K" => temperatures, "undetermined" => Int[], "R_u" => 8314.46261815324)
    payload["wind_rules"] = Dict("clip_fraction" => 0.7, "composition_switch_km" => 80.0, "switch_band_km" => 1e-9)
    return payload
end

function engine_params(model, wind)
    cfg = withenv("SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "0", "SPACEAGORA_VACUUM_GRAM_CACHE" => "0",
            "SPACEAGORA_GRAM_TRACK_CACHE" => "0", "SPACEAGORA_GRAM_RUNTIME_STATS" => "0") do
        CB._snapshot_callback_env_config()
    end
    return (
        args=(environment_model=(density_model=model, wind=wind, planet=(T_ref=200.0,), EI=100.0),
            dynamics_model=(dynamic_effectors=(),),),
        shared_buffers=(density_models=SpaceAGORA.AbstractDensityModel[], callback_env_config=Ref(cfg),
            gram_density_cache=Union{Nothing,CB.GramTrackCache}[nothing],
            densities=[-1.0], temperatures=[-1.0], winds=[SVector(999.0,999.0,999.0)],
            density_sample_t=[NaN], current_time=Ref(0.0)))
end

@testset "Pinned GRAM near-surface dependency boundary" begin
    @test Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) !== nothing
    @test GRAMSuite._GRAM_WRAPPER[] === nothing
    mktempdir() do dir
        models = map((false, true)) do stored
            file = joinpath(dir, stored ? "format2.jls" : "format1.jls")
            serialize(file, near_surface_payload(; stored_winds=stored))
            EM.GRAMNearSurfaceAtmosphereModel(; surrogate_file=file, expected_sha256=bytes2hex(sha256(read(file))))
        end
        legacy, lower = models
        @test !GRAMSuite.near_surface_winds_available(legacy.core)
        @test GRAMSuite.near_surface_winds_available(lower.core)
        @test legacy.core.winds === nothing
        for h in (3_000.0, 12_000.0), wind in (false, true)
            state = EM.getDensity(lower, h, 0.0, 0.0, 123.0, wind)
            @test state == GRAMSuite.density_state(lower.core, h, 0.0, 0.0, 123.0, wind)
            @test state[1] ≈ exp(log(700.0) - h / 11_000.0) / (192.0 * 200.0)
            @test state[2] == 200.0
            @test state[3] ≈ SVector(10.0, -5.0, 0.0)
            @test EM.getDensity(legacy, h, 0.0, 0.0, 123.0, wind) == (state[1], state[2], SVector(0.0,0.0,0.0))
            @test EM.getDensity(lower, h, 0.0, 0.0, -10.0, wind, nothing) == state
        end
        @test_throws DomainError EM.getDensity(lower, 82_000.0, 0.0, 0.0, 0.0, true)
        @test_throws DomainError EM.getDensity(legacy, 76_000.0, 0.0, 0.0, 0.0, true)

        grid_file = joinpath(dir, "upper.jls")
        serialize(grid_file, Dict("status" => "ok", "type" => "surrogate_trilinear", "planet" => "mars",
            "grid" => Dict("alt_km" => [10.0,20.0], "lat_deg" => [-10.0,10.0], "lon_deg" => [0.0,180.0]),
            "fields" => Dict(k => fill(v,2,2,2) for (k,v) in zip(
                ("density_kgm3","temperature_K","wind_ew_ms","wind_ns_ms","wind_up_ms"), (1e-8,180.0,2.0,3.0,4.0)))))
        upper = EM.GRAMGridAtmosphereModel(; planet="Mars", surrogate_file=grid_file)
        combined = EM.CombinedAtmosphereModel(lower, upper; handover_height_m=10_000.0)
        for (h, component) in ((prevfloat(10_000.0),lower), (10_000.0,upper), (12_000.0,upper))
            @test EM.getDensity(combined,h,0.0,0.0,0.0,true) == EM.getDensity(component,h,0.0,0.0,0.0,true)
        end
        @test_throws DomainError EM.getDensity(combined,20_001.0,0.0,0.0,0.0,true)

        # Exercise the actual engine's snapshot route and written wind buffer,
        # not just the mask helper. Direct frozen queries return stored winds
        # even for wind=false; the environment flag is applied downstream.
        x = SVector(1.0,2.0,3.0,4.0,5.0,6.0,100.0)
        for model in (lower, combined), h in (3_000.0, 12_000.0), enabled in (false, true)
            p = engine_params(model, enabled)
            frame = ES.PlanetFrameSample(SMatrix{3,3,Float64}(I),SVector(1.0,2.0,3.0),SVector(4.0,5.0,6.0),h,0.0,0.0)
            raw = EM.getDensity(model,h,0.0,0.0,0.0,enabled)
            sample = EN._sample_atmosphere_from_planet_frame(x,frame,p,1,0.0)
            @test (sample.rho_kg_m3,sample.temperature_k) == raw[1:2]
            @test sample.wind_pp == (enabled ? raw[3] : SVector(0.0,0.0,0.0))
            @test p.shared_buffers.winds[1] == sample.wind_pp
        end

        # The dependency API upgrade does not add a new named preset or widen
        # the existing format-1 catalog contract. No artifact is downloaded.
        catalog = TOML.parsefile(EM._SURROGATE_CATALOG)
        entry = deepcopy(first(filter(p -> get(p,"kind","gram_grid") == "gram_near_surface_scalars", catalog["presets"])))
        @test entry["required_metadata"]["format"] == "spaceagora_mars_near_surface_scalars_v1"
        catalog["presets"] = [entry]
        entry["payload"]["format"] = entry["required_metadata"]["format"] = "spaceagora_mars_near_surface_v2"
        entry["atmosphere"]["winds_available"] = true
        file = joinpath(dir,"catalog.toml")
        open(io -> TOML.print(io,catalog),file,"w")
        error = try SpaceAGORA.available_surrogate_presets(; catalog_file=file); nothing catch ex; ex end
        @test error isa ArgumentError
        @test occursin("Near-surface presets require payload format", sprint(showerror,error))
    end
    @test GRAMSuite._GRAM_WRAPPER[] === nothing
    @test !any(path -> occursin("libgram",lowercase(path)), Libdl.dllist())
end
