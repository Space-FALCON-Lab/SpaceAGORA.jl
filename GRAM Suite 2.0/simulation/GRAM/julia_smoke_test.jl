function resolve_gram_root()
    gram_root = get(ENV, "GRAM_ROOT", "")
    if isempty(gram_root)
        candidates = [
            normpath(joinpath(@__DIR__, "..", "..")),
            normpath(joinpath(@__DIR__, "..", "..", "GRAM Suite 2.0")),
            normpath(joinpath(@__DIR__, "..", "..", "GRAM")),
        ]
        for c in candidates
            if isdir(joinpath(c, "Build")) && isdir(joinpath(c, "Julia"))
                gram_root = c
                break
            end
        end
    end

    isempty(gram_root) && error("Set GRAM_ROOT environment variable to your GRAM Suite root path.")
    return gram_root
end

gram_root = resolve_gram_root()
include(joinpath(gram_root, "Julia", "GRAM.jl"))
using .GRAM

function main(gram_root::AbstractString)
    libext = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")
    libpath = get(ENV, "GRAM_LIB", joinpath(gram_root, "Build", "lib", "libGRAM.$libext"))
    set_library!(libpath)
    initialize!(joinpath(gram_root, "SPICE"))
    enable_console_output!(false)
    spice_guess = try_get_spice_path()
    spice_guess2, data_guess2 = try_get_data_paths()
    mars_spice_guess = mars_try_get_spice_path()
    println("spice_path_guess_len=$(length(spice_guess))")
    println("data_path_guess_len=$(length(data_guess2))")
    println("mars_spice_path_guess_len=$(length(mars_spice_guess))")
    @assert length(spice_guess2) >= 0

    atmos = create_atmosphere(BODY_MARS; data_path = joinpath(gram_root, "Mars", "data"))
    atmos_direct = create_mars_atmosphere(data_path = joinpath(gram_root, "Mars", "data"))
    close!(atmos_direct)
    atmos_copy = copy_atmosphere(atmos)
    close!(atmos_copy)
    atmos_copy_mars = copy_mars_atmosphere(atmos)
    close!(atmos_copy_mars)
    set_start_time!(atmos; year = 2020, month = 3, day = 15, hour = 0, minute = 0, seconds = 0.0, scale = 1, frame = 1)
    set_delta!(atmos; elapsed_time = 5.0)
    set_seed!(atmos, 1001)
    set_min_relative_step_size!(atmos, 0.5)
    set_perturbation_scales!(atmos; density_scale = 1.2, ew_wind_scale = 1.0, ns_wind_scale = 1.0, vertical_wind_scale = 1.0)
    set_perturbation_action!(atmos, true)
    set_ephemeris_fast_mode!(atmos, true)
    set_subsolar_update_time!(atmos, 0.0)
    set_map_year!(atmos, 0)
    set_planetary_radii!(atmos; equatorial = 3396.19, polar = 3376.2)
    set_height_offset_model!(atmos; model = 0, height_offset = 0.0)
    set_height_above_surface!(atmos, 0.0)
    set_mgcm_dust_levels!(atmos; constant_level = 2.0, min_level = 0.0, max_level = 0.0)
    set_dust_storm!(atmos; longitude_sun = 250.0, duration = 30.0, intensity = 1.0, max_radius = 20.0, latitude = 10.0, longitude = 120.0)
    set_perturbation_wavelength_scale!(atmos, 1.0)
    set_exospheric_temperature!(atmos; offset = 0.0, factor = 1.0)
    set_wave_defaults!(
        atmos;
        date = 0.0, scale = 1.0, mean = 0.0,
        a1 = 0.0, p1 = 0.0, r1 = 0.0,
        a2 = 0.0, p2 = 0.0, r2 = 0.0,
        a3 = 0.0, p3 = 0.0, r3 = 0.0
    )
    set_wind_scales!(atmos; mean_winds = 1.0, boundary_layer_winds = 1.0)
    set_mola_heights!(atmos, false)
    set_min_max!(atmos, 0)
    println("mars_version=$(get_version_string(atmos))")
    set_position!(atmos; height = 50.0, latitude = 22.0, longitude = 48.0, elapsed_time = 100.0)

    err = mars_update!(atmos)
    if err != 0
        error("GRAM update failed: $(mars_get_error_message())")
    end

    st = get_start_time(atmos)
    pos = get_position(atmos)
    dyn = get_dynamics_state(atmos)
    den = get_density_state(atmos)
    wnd = get_winds_state(atmos)
    gas = get_gases_state(atmos)
    ar = get_constituent_gas(atmos, GAS_ARGON)
    eph = get_ephemeris_state(atmos)
    set_ephemeris_state!(atmos, eph)
    pert = get_perturbation_state(atmos)
    mars_daily = get_daily_dynamics_state(atmos)
    mars_state = get_mars_state(atmos)
    mars_gases = get_mars_gases_state(atmos)
    println("GRAM smoke test passed")
    println("start_year=$(st.year)")
    println("height_km=$(pos.height)")
    println("temperature_K=$(dyn.temperature)")
    println("density_kg_m3=$(dyn.density)")
    println("density_state_kg_m3=$(den.density)")
    println("winds_ew_m_s=$(wnd.ewWind)")
    println("avg_mw=$(gas.averageMolecularWeight)")
    println("argon_mole_fraction=$(ar.moleFraction)")
    println("solar_time=$(eph.solarTime)")
    println("density_rand=$(pert.densityRandomNumber)")
    println("daily_temperature=$(mars_daily.temperatureDaily)")
    println("dust_optical_depth=$(mars_state.dustOpticalDepth)")
    println("mars_gases_avg_mw=$(mars_gases.state.averageMolecularWeight)")

    track = generate_trajectory(
        atmos;
        initial_height = 50.0,
        initial_latitude = 22.0,
        initial_longitude = 48.0,
        initial_elapsed_time = 0.0,
        delta_height = 0.1,
        delta_latitude = 0.02,
        delta_longitude = 0.01,
        delta_elapsed_time = 5.0,
        n_points = 8
    )
    println("trajectory_points=$(length(track))")
    println("trajectory_first_density_kg_m3=$(track[1].dynamics.density)")

    tracks_mc = generate_monte_carlo_trajectories(
        atmos;
        initial_height = 50.0,
        initial_latitude = 22.0,
        initial_longitude = 48.0,
        initial_elapsed_time = 0.0,
        delta_height = 0.1,
        delta_latitude = 0.02,
        delta_longitude = 0.01,
        delta_elapsed_time = 5.0,
        n_points = 8,
        n_runs = 3,
        initial_seed = 1001
    )
    println("trajectory_montecarlo_runs=$(length(tracks_mc))")
    println("trajectory_montecarlo_first_density_kg_m3=$(tracks_mc[1][1].dynamics.density)")

    close!(atmos)
end

main(gram_root)
