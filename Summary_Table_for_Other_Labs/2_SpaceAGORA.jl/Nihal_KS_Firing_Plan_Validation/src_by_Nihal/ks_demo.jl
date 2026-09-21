function ks_state_error_timeseries(
    truth_states::AbstractVector{<:AbstractVector{<:Real}},
    approx_states::AbstractVector{<:AbstractVector{<:Real}},
)
    length(truth_states) == length(approx_states) ||
        throw(ArgumentError("truth_states and approx_states must have the same length."))

    nstate = length(truth_states[1])
    nsamples = length(truth_states)
    errors = zeros(Float64, nstate, nsamples)

    for k in 1:nsamples
        length(truth_states[k]) == nstate || throw(ArgumentError("truth_states must share a common state dimension."))
        length(approx_states[k]) == nstate || throw(ArgumentError("approx_states must share a common state dimension."))
        errors[:, k] .= Float64.(truth_states[k]) .- Float64.(approx_states[k])
    end

    return errors
end

function ks_reduced_state_error_timeseries(
    state_errors::AbstractMatrix{<:Real},
    nsat::Int,
)
    size(state_errors, 1) == 10 * nsat || throw(ArgumentError("state_errors must have 10 rows per satellite."))
    nsamples = size(state_errors, 2)
    reduced = zeros(Float64, 4 * nsat, nsamples)

    for sat in 1:nsat
        src = 10 * (sat - 1)
        dst = 4 * (sat - 1)
        reduced[dst + 1, :] .= sqrt.(vec(sum(abs2, state_errors[(src + 1):(src + 4), :]; dims = 1)))
        reduced[dst + 2, :] .= sqrt.(vec(sum(abs2, state_errors[(src + 5):(src + 8), :]; dims = 1)))
        reduced[dst + 3, :] .= state_errors[src + 9, :]
        reduced[dst + 4, :] .= state_errors[src + 10, :]
    end

    return reduced
end

function ks_reduced_state_component_labels(nsat::Int)
    base_labels = ("||Δp||", "||Δq||", "Δh", "Δt")
    labels = String[]
    for sat in 1:nsat
        for label in base_labels
            push!(labels, "Sat $sat $label")
        end
    end
    return labels
end

function plot_ks_state_error_grid(
    physical_times::AbstractVector{<:Real},
    state_errors::AbstractMatrix{<:Real};
    component_labels::AbstractVector{<:AbstractString},
    case_label::AbstractString,
    filename::AbstractString,
)
    nstate, nsamples = size(state_errors)
    length(physical_times) == nsamples || throw(ArgumentError("physical_times length must match the number of state samples."))
    length(component_labels) == nstate || throw(ArgumentError("component_labels length must match the number of state rows."))

    ncols = 2
    nrows = cld(nstate, ncols)
    plt = plot(
        layout = (nrows, ncols),
        size = (1800, 340 * nrows),
        legend = false,
    )

    for idx in 1:nstate
        plot!(
            plt,
            physical_times,
            state_errors[idx, :];
            subplot = idx,
            linewidth = 2.5,
            title = component_labels[idx],
            xlabel = "Time [s]",
            ylabel = "Truth - linearized",
            formatter = :scientific,
        )
    end

    plot!(plt; plot_title = "KS State Error Comparison: $case_label")
    savefig(plt, filename)
    return plt
end

function ks_cartesian_trajectory(sol)
    nsamples = length(sol.u)
    nsat = ks_satellite_count(sol.u[1])
    positions = zeros(Float64, 3, nsamples, nsat)
    physical_times = zeros(Float64, nsamples)

    for k in 1:nsamples
        physical_times[k] = ks_mean_physical_time(sol.u[k])
        for sat in 1:nsat
            positions[:, k, sat], _, _, _ = ks_cartesian_components(sol.u[k], sat)
        end
    end

    return positions, physical_times
end

function ks_plot_earth_sphere!(plt; radius_km::Real = R_EARTH_KM)
    theta = range(0.0, 2pi; length = 40)
    phi = range(0.0, pi; length = 20)
    x = [float(radius_km) * cos(th) * sin(ph) for ph in phi, th in theta]
    y = [float(radius_km) * sin(th) * sin(ph) for ph in phi, th in theta]
    z = [float(radius_km) * cos(ph) for ph in phi, th in theta]

    surface!(
        plt,
        x,
        y,
        z;
        color = cgrad([RGB(0.2, 0.45, 0.85), RGB(0.65, 0.85, 1.0)]),
        fillalpha = 0.55,
        linealpha = 0.08,
        colorbar = false,
        label = "Earth",
    )
    return plt
end

function run_ks_two_satellite_orbit_animation(;
    output_dir = "outputs/ks",
    dt_s = 45.0,
    num_orbits = 0.35,
    anomaly_spacing_deg = 0.5,
    altitude_km = 550.0,
    inclination_deg = 53.0,
    eccentricity = 0.0,
    argp_deg = 0.0,
    raan_deg = 0.0,
    nu0_deg = 0.0,
    link_max_range_km = 200.0,
    atmosphere_top_km = 100.0,
    cr = 2.0,
    laser_power_w = 1.0e4,
    satellite_mass_kg = 300.0,
    frame_stride = 2,
    fps = 15,
)
    init = initialize_coplanar_ks_constellation(
        2;
        altitude_km = altitude_km,
        inclination_deg = inclination_deg,
        eccentricity = eccentricity,
        argp_deg = argp_deg,
        raan_deg = raan_deg,
        nu0_deg = nu0_deg,
        anomaly_spacing_deg = anomaly_spacing_deg,
    )

    orbital_period_s = 2pi * sqrt(init.a_km^3 / MU_EARTH_KM3_S2)
    physical_tspan = (0.0, float(num_orbits) * orbital_period_s)
    nominal_s_step = float(dt_s) / init.a_km
    controls = ones(Float64, init.N)

    sol = propagate_ks_constellation(
        init.u0,
        physical_tspan;
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = controls,
        saveat_sundman = nominal_s_step,
    )

    positions, physical_times = ks_cartesian_trajectory(sol)
    nsamples = size(positions, 2)
    frame_indices = collect(1:max(1, frame_stride):nsamples)
    last(frame_indices) == nsamples || push!(frame_indices, nsamples)

    max_extent = maximum(abs, positions)
    axis_limit = 1.15 * max(max_extent, R_EARTH_KM)

    project_root = normpath(joinpath(@__DIR__, ".."))
    output_dir_resolved = isabspath(output_dir) ? output_dir : normpath(joinpath(project_root, output_dir))
    mkpath(output_dir_resolved)
    gif_path = joinpath(output_dir_resolved, "ks_two_satellite_orbits.gif")

    anim = @animate for k in frame_indices
        r1 = positions[:, k, 1]
        r2 = positions[:, k, 2]
        pair_data = ks_pair_geometry_data(
            sol.u[k];
            link_max_range_km = link_max_range_km,
            atmosphere_top_km = atmosphere_top_km,
        )
        link_active = !isempty(pair_data) && pair_data[1].zeta == 1.0

        plt = plot(
            size = (900, 800),
            legend = :topright,
            background_color = :white,
            plot_background_color = RGB(0.98, 0.99, 1.0),
            foreground_color = :black,
            xlabel = "x [km]",
            ylabel = "y [km]",
            zlabel = "z [km]",
            xlims = (-axis_limit, axis_limit),
            ylims = (-axis_limit, axis_limit),
            zlims = (-axis_limit, axis_limit),
            aspect_ratio = :equal,
            camera = (35, 25),
            title = "Two-Satellite Cartesian Orbit Animation\n t = $(round(physical_times[k]; digits = 1)) s",
        )

        ks_plot_earth_sphere!(plt)

        plot!(
            plt,
            positions[1, 1:k, 1],
            positions[2, 1:k, 1],
            positions[3, 1:k, 1];
            color = :crimson,
            linewidth = 2.5,
            label = "Satellite 1 path",
        )
        plot!(
            plt,
            positions[1, 1:k, 2],
            positions[2, 1:k, 2],
            positions[3, 1:k, 2];
            color = :darkorange,
            linewidth = 2.5,
            label = "Satellite 2 path",
        )

        scatter!(
            plt,
            [r1[1]],
            [r1[2]],
            [r1[3]];
            color = :crimson,
            markerstrokecolor = :white,
            markerstrokewidth = 0.8,
            markersize = 7,
            label = "Satellite 1",
        )
        scatter!(
            plt,
            [r2[1]],
            [r2[2]],
            [r2[3]];
            color = :darkorange,
            markerstrokecolor = :white,
            markerstrokewidth = 0.8,
            markersize = 7,
            label = "Satellite 2",
        )

        if link_active
            plot!(
                plt,
                [r1[1], r2[1]],
                [r1[2], r2[2]],
                [r1[3], r2[3]];
                color = :deepskyblue3,
                linestyle = :dash,
                linewidth = 2.5,
                label = "Laser link",
            )
        end

        plt
    end

    gif(anim, gif_path; fps = fps)
    println("Saved KS orbit animation: $gif_path")

    return (
        sol = sol,
        positions = positions,
        physical_times = physical_times,
        gif_path = gif_path,
    )
end

function ks_max_matrix_error(error_matrix::AbstractMatrix{<:Real})
    linear_idx = argmax(error_matrix)
    cart_idx = Tuple(CartesianIndices(error_matrix)[linear_idx])
    return (
        value = error_matrix[linear_idx],
        row = cart_idx[1],
        col = cart_idx[2],
    )
end

function plot_ks_jacobian_error_heatmaps(
    a_error::AbstractMatrix{<:Real},
    b_error::AbstractMatrix{<:Real};
    filename::AbstractString,
    case_label::AbstractString,
)
    a_plot = Float64.(a_error)
    b_plot = Float64.(b_error)
    clim_lo = min(minimum(a_plot), minimum(b_plot))
    clim_hi = max(maximum(a_plot), maximum(b_plot))
    if isapprox(clim_lo, clim_hi; atol = 1e-12, rtol = 0.0)
        clim_lo -= 1.0
        clim_hi += 1.0
    end
    clims = (
        clim_lo,
        clim_hi,
    )

    a_max = ks_max_matrix_error(a_error)
    b_max = ks_max_matrix_error(b_error)

    plt = plot(
        layout = (1, 2),
        size = (1800, 700),
        background_color = :white,
        plot_background_color = :white,
        foreground_color = :black,
        foreground_color_text = :black,
        foreground_color_axis = :black,
        foreground_color_grid = RGBA(0, 0, 0, 0.15),
    )
    heatmap!(
        plt,
        1:size(a_plot, 2),
        1:size(a_plot, 1),
        a_plot;
        subplot = 1,
        xlabel = "Column",
        ylabel = "Row",
        title = "Ac abs error\nmax=$(round(a_max.value; sigdigits = 4)) @ ($(a_max.row), $(a_max.col))",
        clims = clims,
        color = cgrad(:viridis),
        colorbar_title = "abs error",
        aspect_ratio = :equal,
        yflip = true,
    )
    heatmap!(
        plt,
        1:size(b_plot, 2),
        1:size(b_plot, 1),
        b_plot;
        subplot = 2,
        xlabel = "Column",
        ylabel = "Row",
        title = "Bc abs error\nmax=$(round(b_max.value; sigdigits = 4)) @ ($(b_max.row), $(b_max.col))",
        clims = clims,
        color = cgrad(:viridis),
        colorbar_title = "abs error",
        aspect_ratio = :auto,
        xlims = (0.5, size(b_plot, 2) + 0.5),
        xticks = 1:size(b_plot, 2),
        yflip = true,
    )
    plot!(plt; plot_title = "KS Continuous Jacobian Verification: $case_label")
    savefig(plt, filename)
    return plt
end

function run_ks_jacobian_verification_demo(;
    output_dir = "outputs/ks",
    anomaly_spacing_deg = 0.5,
    altitude_km = 550.0,
    inclination_deg = 53.0,
    eccentricity = 0.0,
    argp_deg = 0.0,
    raan_deg = 0.0,
    nu0_deg = 0.0,
    link_max_range_km = 200.0,
    atmosphere_top_km = 100.0,
    cr = 2.0,
    laser_power_w = 1.0e4,
    satellite_mass_kg = 300.0,
)
    init = initialize_coplanar_ks_constellation(
        2;
        altitude_km = altitude_km,
        inclination_deg = inclination_deg,
        eccentricity = eccentricity,
        argp_deg = argp_deg,
        raan_deg = raan_deg,
        nu0_deg = nu0_deg,
        anomaly_spacing_deg = anomaly_spacing_deg,
    )

    u = Float64.(init.u0)
    controls = ones(Float64, ks_satellite_count(u))
    analytic = ks_continuous_jacobians(
        u;
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = controls,
    )
    fd = ks_forwarddiff_continuous_jacobians(
        u;
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = controls,
    )

    first_pair = analytic.pair_data[1]
    first_pair.zeta == 1.0 ||
        throw(ArgumentError("Chosen verification geometry does not activate the laser link."))

    a_error = abs.(analytic.A .- fd.A)
    b_error = abs.(analytic.B .- fd.B)
    a_max = ks_max_matrix_error(a_error)
    b_max = ks_max_matrix_error(b_error)

    project_root = normpath(joinpath(@__DIR__, ".."))
    output_dir_resolved = isabspath(output_dir) ? output_dir : normpath(joinpath(project_root, output_dir))
    mkpath(output_dir_resolved)
    plot_path = joinpath(output_dir_resolved, "ks_continuous_jacobian_forwarddiff_error.png")
    plot_ks_jacobian_error_heatmaps(
        a_error,
        b_error;
        filename = plot_path,
        case_label = "J2 + laser on",
    )

    println("Saved KS Jacobian verification plot: $plot_path")
    println("Laser-on pair distance [km]: $(round(first_pair.distance_km; digits = 6))")
    println("Max |Ac_analytic - Ac_forwarddiff|: $(a_max.value) at ($(a_max.row), $(a_max.col))")
    println("Max |Bc_analytic - Bc_forwarddiff|: $(b_max.value) at ($(b_max.row), $(b_max.col))")

    return (
        analytic = analytic,
        forwarddiff = fd,
        a_error = a_error,
        b_error = b_error,
        plot_path = plot_path,
        max_a_error = a_max,
        max_b_error = b_max,
        init = init,
    )
end

function run_ks_radial_control_linearization_verification_demo(;
    output_dir = "outputs/ks",
    anomaly_spacing_deg = 17.5,
    altitude_km = 550.0,
    inclination_deg = 53.0,
    eccentricity = 0.01,
    argp_deg = 0.0,
    raan_deg = 0.0,
    nu0_deg = 35.0,
    radial_control_accelerations = [2.0e-8, -1.5e-8],
    state_step = 1e-6,
    control_step = 1e-8,
    kwargs...,
)
    init = initialize_coplanar_ks_constellation(
        2;
        altitude_km = altitude_km,
        inclination_deg = inclination_deg,
        eccentricity = eccentricity,
        argp_deg = argp_deg,
        raan_deg = raan_deg,
        nu0_deg = nu0_deg,
        anomaly_spacing_deg = anomaly_spacing_deg,
    )

    u = Float64.(init.u0)
    controls = Float64.(radial_control_accelerations)
    length(controls) == ks_satellite_count(u) ||
        throw(ArgumentError("radial_control_accelerations must have one entry per satellite."))

    analytic = ks_radial_control_continuous_jacobians(
        u;
        radial_control_accelerations = controls,
    )
    ad = ks_forwarddiff_radial_control_continuous_jacobians(
        u;
        radial_control_accelerations = controls,
    )
    fd = ks_finite_difference_radial_control_continuous_jacobians(
        u;
        radial_control_accelerations = controls,
        state_step = state_step,
        control_step = control_step,
    )

    ad_a_error = abs.(analytic.A .- ad.A)
    ad_b_error = abs.(analytic.B .- ad.B)
    fd_a_error = abs.(analytic.A .- fd.A)
    fd_b_error = abs.(analytic.B .- fd.B)
    ad_a_max = ks_max_matrix_error(ad_a_error)
    ad_b_max = ks_max_matrix_error(ad_b_error)
    fd_a_max = ks_max_matrix_error(fd_a_error)
    fd_b_max = ks_max_matrix_error(fd_b_error)

    project_root = normpath(joinpath(@__DIR__, ".."))
    output_dir_resolved = isabspath(output_dir) ? output_dir : normpath(joinpath(project_root, output_dir))
    mkpath(output_dir_resolved)
    forwarddiff_plot_path = joinpath(output_dir_resolved, "ks_radial_control_forwarddiff_error.png")
    finite_difference_plot_path = joinpath(output_dir_resolved, "ks_radial_control_finite_difference_error.png")

    plot_ks_jacobian_error_heatmaps(
        ad_a_error,
        ad_b_error;
        filename = forwarddiff_plot_path,
        case_label = "J2 + radial control vs ForwardDiff",
    )
    plot_ks_jacobian_error_heatmaps(
        fd_a_error,
        fd_b_error;
        filename = finite_difference_plot_path,
        case_label = "J2 + radial control vs finite difference",
    )

    println("Saved KS radial-control ForwardDiff error plot: $forwarddiff_plot_path")
    println("Saved KS radial-control finite-difference error plot: $finite_difference_plot_path")
    println("Radial control accelerations [km/s^2]: $controls")
    println("Max |Ac_radial_analytic - Ac_radial_forwarddiff|: $(ad_a_max.value) at ($(ad_a_max.row), $(ad_a_max.col))")
    println("Max |Bc_radial_analytic - Bc_radial_forwarddiff|: $(ad_b_max.value) at ($(ad_b_max.row), $(ad_b_max.col))")
    println("Max |Ac_radial_analytic - Ac_radial_finite_difference|: $(fd_a_max.value) at ($(fd_a_max.row), $(fd_a_max.col))")
    println("Max |Bc_radial_analytic - Bc_radial_finite_difference|: $(fd_b_max.value) at ($(fd_b_max.row), $(fd_b_max.col))")

    return (
        analytic = analytic,
        forwarddiff = ad,
        finite_difference = fd,
        forwarddiff_a_error = ad_a_error,
        forwarddiff_b_error = ad_b_error,
        finite_difference_a_error = fd_a_error,
        finite_difference_b_error = fd_b_error,
        forwarddiff_plot_path = forwarddiff_plot_path,
        finite_difference_plot_path = finite_difference_plot_path,
        max_forwarddiff_a_error = ad_a_max,
        max_forwarddiff_b_error = ad_b_max,
        max_finite_difference_a_error = fd_a_max,
        max_finite_difference_b_error = fd_b_max,
        init = init,
    )
end

function run_ks_cartesian_controllability_verification_demo(;
    no_laser_anomaly_spacing_deg = 0.0,
    laser_anomaly_spacing_deg = 0.0,
    altitude_km = 550.0,
    inclination_deg = 0.0,
    eccentricity = 0.0,
    argp_deg = 0.0,
    raan_deg = 0.0,
    nu0_deg = 0.0,
    altitude_offsets_km = [0.0, 50.0],
    j2 = J2_EARTH,
    atmosphere_top_km = 100.0,
    cr = 2.0,
    laser_power_w = 1.0e12,
    satellite_mass_kg = 300.0,
    rtol = 1e-9,
    residual_tolerance = 1e-8,
    require_match = false,
    kwargs...,
)
    no_laser_init = initialize_coplanar_ks_constellation(
        2;
        altitude_km = altitude_km,
        inclination_deg = inclination_deg,
        eccentricity = eccentricity,
        argp_deg = argp_deg,
        raan_deg = raan_deg,
        nu0_deg = nu0_deg,
        anomaly_spacing_deg = no_laser_anomaly_spacing_deg,
        altitude_offsets_km = altitude_offsets_km,
    )
    laser_init = initialize_coplanar_ks_constellation(
        2;
        altitude_km = altitude_km,
        inclination_deg = inclination_deg,
        eccentricity = eccentricity,
        argp_deg = argp_deg,
        raan_deg = raan_deg,
        nu0_deg = nu0_deg,
        anomaly_spacing_deg = laser_anomaly_spacing_deg,
        altitude_offsets_km = altitude_offsets_km,
    )

    no_laser = ks_cartesian_controllability_comparison(
        no_laser_init.u0;
        j2 = j2,
        link_max_range_km = 1.0,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = ones(Float64, 2),
        rtol = rtol,
    )
    laser = ks_cartesian_controllability_comparison(
        laser_init.u0;
        j2 = j2,
        link_max_range_km = 200.0,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = ones(Float64, 2),
        rtol = rtol,
    )

    println("KS/Cartesian controllability comparison: no laser interaction")
    println("  active links: $(sum(data.zeta for data in no_laser.pair_data))")
    println("  KS rank: $(no_laser.ks_rank)")
    println("  projected KS Cartesian rank: $(no_laser.projected_ks_cartesian_rank)")
    println("  Cartesian rank: $(no_laser.cartesian_rank)")
    println("  input projection error: $(no_laser.input_projection_error)")
    println("  KS -> Cartesian subspace residual: $(no_laser.ks_to_cartesian_residual)")
    println("  Cartesian -> KS subspace residual: $(no_laser.cartesian_to_ks_residual)")
    println("KS/Cartesian controllability comparison: laser interaction")
    println("  active links: $(sum(data.zeta for data in laser.pair_data))")
    println("  KS rank: $(laser.ks_rank)")
    println("  projected KS Cartesian rank: $(laser.projected_ks_cartesian_rank)")
    println("  Cartesian rank: $(laser.cartesian_rank)")
    println("  input projection error: $(laser.input_projection_error)")
    println("  KS -> Cartesian subspace residual: $(laser.ks_to_cartesian_residual)")
    println("  Cartesian -> KS subspace residual: $(laser.cartesian_to_ks_residual)")

    if require_match
        no_laser.rank_match || throw(ArgumentError("No-laser projected KS rank does not match Cartesian rank."))
        laser.rank_match || throw(ArgumentError("Laser projected KS rank does not match Cartesian rank."))
        no_laser.ks_to_cartesian_residual <= residual_tolerance ||
            throw(ArgumentError("No-laser KS-to-Cartesian controllability subspace residual is too large."))
        no_laser.cartesian_to_ks_residual <= residual_tolerance ||
            throw(ArgumentError("No-laser Cartesian-to-KS controllability subspace residual is too large."))
        laser.ks_to_cartesian_residual <= residual_tolerance ||
            throw(ArgumentError("Laser KS-to-Cartesian controllability subspace residual is too large."))
        laser.cartesian_to_ks_residual <= residual_tolerance ||
            throw(ArgumentError("Laser Cartesian-to-KS controllability subspace residual is too large."))
    end

    return (
        no_laser = no_laser,
        laser = laser,
        no_laser_init = no_laser_init,
        laser_init = laser_init,
    )
end

function ks_controllability_rank_history(
    a::AbstractMatrix{<:Real},
    b::AbstractMatrix{<:Real};
    projection::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
    horizon::Union{Nothing, Int} = nothing,
    rtol::Real = 1e-9,
    atol::Real = 0.0,
)
    n = size(a, 1)
    steps = isnothing(horizon) ? n : horizon
    a_mat = Matrix{Float64}(a)
    q = ks_orthonormal_column_basis(b; rtol = rtol, atol = atol)
    ranks = zeros(Int, steps)

    for k in 1:steps
        basis_for_rank = isnothing(projection) ? q : ks_orthonormal_column_basis(projection * q; rtol = rtol, atol = atol)
        ranks[k] = size(basis_for_rank, 2)
        k == steps && break
        size(q, 2) == 0 && continue
        q = ks_orthonormal_column_basis([q a_mat * q]; rtol = rtol, atol = atol)
    end

    return ranks
end

function ks_principal_angles_deg(
    q_a::AbstractMatrix{<:Real},
    q_b::AbstractMatrix{<:Real},
)
    (size(q_a, 2) == 0 || size(q_b, 2) == 0) && return Float64[]
    singular_values = svdvals(transpose(Matrix{Float64}(q_a)) * Matrix{Float64}(q_b))
    return acosd.(clamp.(singular_values, -1.0, 1.0))
end

function ks_rtn_input_direction_metric(
    b_rtn::AbstractMatrix{<:Real},
    sat::Int,
    control::Int,
)
    rows = (6 * (sat - 1) + 4):(6 * (sat - 1) + 6)
    components = Float64.(b_rtn[rows, control])
    component_norm = norm(components)
    fractions = component_norm == 0.0 ? zeros(Float64, 3) : abs.(components) ./ component_norm
    dominant_index = argmax(fractions)

    return (
        components = components,
        fractions = fractions,
        dominant_axis = ("R", "T", "N")[dominant_index],
        dominant_fraction = fractions[dominant_index],
    )
end

function ks_rtn_direction_test_case(
    name::AbstractString;
    expected_axis::AbstractString,
    altitude_offsets_km::AbstractVector{<:Real},
    anomaly_spacing_deg::Real,
    inclination_deg::Real = 45.0,
    altitude_km::Real = 550.0,
    j2::Real = J2_EARTH,
    link_max_range_km::Real = 200.0,
    rtol::Real = 1e-9,
    minimum_fraction::Real = 0.99,
)
    init = initialize_coplanar_ks_constellation(
        2;
        altitude_km = altitude_km,
        inclination_deg = inclination_deg,
        eccentricity = 0.0,
        anomaly_spacing_deg = anomaly_spacing_deg,
        altitude_offsets_km = altitude_offsets_km,
    )

    rtn = ks_laser_rtn_linearization(
        init.u0;
        j2 = j2,
        link_max_range_km = link_max_range_km,
        control_amplitudes = ones(Float64, 2),
        rtol = rtol,
    )

    pair = rtn.comparison.pair_data[1]
    pair.zeta == 1.0 || throw(ArgumentError("RTN direction test case $name does not have an active laser link."))

    sat1_metric = ks_rtn_input_direction_metric(rtn.B_RTN, 1, 1)
    sat2_metric = ks_rtn_input_direction_metric(rtn.B_RTN, 2, 2)
    classical_rank = ks_numerical_rank(rtn.controllability; rtol = rtol)
    passes = sat1_metric.dominant_axis == expected_axis &&
        sat2_metric.dominant_axis == expected_axis &&
        sat1_metric.dominant_fraction >= minimum_fraction &&
        sat2_metric.dominant_fraction >= minimum_fraction

    return (
        name = String(name),
        expected_axis = String(expected_axis),
        passes = passes,
        distance_km = pair.distance_km,
        inclination_deg = float(inclination_deg),
        anomaly_spacing_deg = float(anomaly_spacing_deg),
        altitude_offsets_km = Float64.(altitude_offsets_km),
        classical_rank = classical_rank,
        iterative_rank = rtn.controllability_rank,
        sat1_control_components = sat1_metric.components,
        sat1_control_fractions = sat1_metric.fractions,
        sat1_dominant_axis = sat1_metric.dominant_axis,
        sat1_dominant_fraction = sat1_metric.dominant_fraction,
        sat2_control_components = sat2_metric.components,
        sat2_control_fractions = sat2_metric.fractions,
        sat2_dominant_axis = sat2_metric.dominant_axis,
        sat2_dominant_fraction = sat2_metric.dominant_fraction,
        rtn = rtn,
    )
end

function run_ks_rtn_direction_controllability_tests(;
    inclination_deg::Real = 45.0,
    altitude_km::Real = 550.0,
    link_distance_km::Real = 50.0,
    j2::Real = J2_EARTH,
    rtol::Real = 1e-9,
    minimum_fraction::Real = 0.99,
)
    a_km = R_EARTH_KM + float(altitude_km)
    tangential_spacing_deg = rad2deg(2.0 * asin(float(link_distance_km) / (2.0 * a_km)))
    cases = (
        ks_rtn_direction_test_case(
            "radial_$(round(link_distance_km; digits = 3))km";
            expected_axis = "R",
            altitude_offsets_km = [0.0, float(link_distance_km)],
            anomaly_spacing_deg = 0.0,
            inclination_deg = inclination_deg,
            altitude_km = altitude_km,
            j2 = j2,
            rtol = rtol,
            minimum_fraction = minimum_fraction,
        ),
        ks_rtn_direction_test_case(
            "tangential_$(round(link_distance_km; digits = 3))km";
            expected_axis = "T",
            altitude_offsets_km = [0.0, 0.0],
            anomaly_spacing_deg = tangential_spacing_deg,
            inclination_deg = inclination_deg,
            altitude_km = altitude_km,
            j2 = j2,
            rtol = rtol,
            minimum_fraction = minimum_fraction,
        ),
    )

    for case in cases
        println("RTN direction test: $(case.name)")
        println("  pass: $(case.passes)")
        println("  link distance [km]: $(case.distance_km)")
        println("  expected dominant axis: $(case.expected_axis)")
        println("  sat 1 control components [v_R, v_T, v_N]: $(case.sat1_control_components)")
        println("  sat 1 fractions [R, T, N]: $(case.sat1_control_fractions)")
        println("  sat 1 dominant axis: $(case.sat1_dominant_axis), fraction: $(case.sat1_dominant_fraction)")
        println("  sat 2 control components [v_R, v_T, v_N]: $(case.sat2_control_components)")
        println("  sat 2 fractions [R, T, N]: $(case.sat2_control_fractions)")
        println("  sat 2 dominant axis: $(case.sat2_dominant_axis), fraction: $(case.sat2_dominant_fraction)")
        println("  classical RTN controllability rank: $(case.classical_rank)")
        println("  iterative RTN controllability rank: $(case.iterative_rank)")
    end

    all(case.passes for case in cases) ||
        throw(ArgumentError("At least one RTN direction controllability test failed."))

    return cases
end

function run_ks_laser_controllability_visualization(;
    output_dir = "outputs/ks",
    anomaly_spacing_deg = 0.5,
    altitude_km = 550.0,
    inclination_deg = 53.0,
    eccentricity = 0.0,
    argp_deg = 0.0,
    raan_deg = 0.0,
    nu0_deg = 0.0,
    atmosphere_top_km = 100.0,
    link_max_range_km = 200.0,
    cr = 2.0,
    laser_power_w = 1.0e12,
    satellite_mass_kg = 300.0,
    rtol = 1e-9,
)
    init = initialize_coplanar_ks_constellation(
        2;
        altitude_km = altitude_km,
        inclination_deg = inclination_deg,
        eccentricity = eccentricity,
        argp_deg = argp_deg,
        raan_deg = raan_deg,
        nu0_deg = nu0_deg,
        anomaly_spacing_deg = anomaly_spacing_deg,
    )

    comparison = ks_cartesian_controllability_comparison(
        init.u0;
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = ones(Float64, 2),
        rtol = rtol,
    )

    pair = comparison.pair_data[1]
    pair.zeta == 1.0 || throw(ArgumentError("Visualization case must have one active laser link."))

    x = comparison.cartesian_state
    r1 = x[1:3]
    r2 = x[7:9]
    midpoint = 0.5 .* (r1 .+ r2)
    span_km = max(1.5 * norm(r2 .- r1), 5.0)
    a_km = R_EARTH_KM + float(altitude_km)
    i_rad = deg2rad(float(inclination_deg))
    raan_rad = deg2rad(float(raan_deg))
    argp_rad = deg2rad(float(argp_deg))
    nu0_rad = deg2rad(float(nu0_deg))
    orbit_nu = range(0.0, 2pi; length = 240)
    orbit_points = reduce(hcat, first(coe_to_rv(a_km, eccentricity, i_rad, raan_rad, argp_rad, nu; mu = MU_EARTH_KM3_S2)) for nu in orbit_nu)
    local_nu = range(-deg2rad(3.0), deg2rad(3.0); length = 80)
    sat1_arc = reduce(hcat, first(coe_to_rv(a_km, eccentricity, i_rad, raan_rad, argp_rad, nu0_rad + dnu; mu = MU_EARTH_KM3_S2)) for dnu in local_nu)
    sat2_arc = reduce(hcat, first(coe_to_rv(a_km, eccentricity, i_rad, raan_rad, argp_rad, nu0_rad + deg2rad(float(anomaly_spacing_deg)) + dnu; mu = MU_EARTH_KM3_S2)) for dnu in local_nu)

    orbit_limit = 1.15 * a_km
    orbit_plot = plot(
        orbit_points[1, :],
        orbit_points[2, :],
        orbit_points[3, :];
        linewidth = 2.5,
        color = :gray45,
        label = "orbital path",
        xlabel = "x [km]",
        ylabel = "y [km]",
        zlabel = "z [km]",
        xlims = (-orbit_limit, orbit_limit),
        ylims = (-orbit_limit, orbit_limit),
        zlims = (-orbit_limit, orbit_limit),
        camera = (35, 25),
        legend = :topright,
        aspect_ratio = :equal,
        title = "Orbit context",
    )
    scatter!(orbit_plot, [0.0], [0.0], [0.0]; markersize = 7, color = :deepskyblue4, label = "Earth center")
    plot!(orbit_plot, [r1[1], r2[1]], [r1[2], r2[2]], [r1[3], r2[3]]; linewidth = 3, color = :deepskyblue3, label = "active laser link")
    scatter!(orbit_plot, [r1[1]], [r1[2]], [r1[3]]; markersize = 6, color = :crimson, label = "sat 1")
    scatter!(orbit_plot, [r2[1]], [r2[2]], [r2[3]]; markersize = 6, color = :darkorange, label = "sat 2")

    geom_plot = plot(
        sat1_arc[1, :],
        sat1_arc[2, :],
        sat1_arc[3, :];
        linewidth = 3,
        color = :crimson,
        label = "sat 1 local orbit",
        xlabel = "x [km]",
        ylabel = "y [km]",
        zlabel = "z [km]",
        xlims = (midpoint[1] - span_km, midpoint[1] + span_km),
        ylims = (midpoint[2] - span_km, midpoint[2] + span_km),
        zlims = (midpoint[3] - span_km, midpoint[3] + span_km),
        camera = (35, 25),
        title = "Active laser geometry\nrange = $(round(pair.distance_km; digits = 3)) km",
    )
    plot!(geom_plot, sat2_arc[1, :], sat2_arc[2, :], sat2_arc[3, :]; linewidth = 3, color = :darkorange, label = "sat 2 local orbit")
    plot!(
        geom_plot,
        [r1[1], r2[1]],
        [r1[2], r2[2]],
        [r1[3], r2[3]];
        linewidth = 3,
        color = :deepskyblue3,
        marker = :circle,
        markersize = 7,
        label = "active laser link",
    )

    horizon = max(size(comparison.ks.A, 1), size(comparison.cartesian.A, 1))
    ks_rank_history = ks_controllability_rank_history(
        comparison.ks.A,
        comparison.ks.B;
        horizon = horizon,
        rtol = rtol,
    )
    projected_rank_history = ks_controllability_rank_history(
        comparison.ks.A,
        comparison.ks.B;
        projection = comparison.cartesian_projection,
        horizon = horizon,
        rtol = rtol,
    )
    cartesian_rank_history = ks_controllability_rank_history(
        comparison.cartesian.A,
        comparison.cartesian.B;
        horizon = horizon,
        rtol = rtol,
    )

    rank_plot = plot(
        1:horizon,
        projected_rank_history;
        linewidth = 3,
        marker = :circle,
        label = "projected KS rank",
        xlabel = "Krylov step",
        ylabel = "reachable rank",
        ylims = (0, max(maximum(ks_rank_history), maximum(cartesian_rank_history)) + 1),
        title = "Reachable rank growth",
    )
    plot!(rank_plot, 1:horizon, cartesian_rank_history; linewidth = 3, linestyle = :dash, marker = :diamond, label = "Cartesian rank")
    plot!(rank_plot, 1:horizon, ks_rank_history; linewidth = 2, linestyle = :dot, label = "raw KS rank")

    angles = ks_principal_angles_deg(
        comparison.projected_ks_cartesian_basis,
        comparison.cartesian_controllability_basis,
    )
    angle_plot = bar(
        1:length(angles),
        angles;
        label = false,
        xlabel = "canonical angle index",
        ylabel = "angle [deg]",
        title = "Projected KS vs Cartesian reachable subspaces\nmax angle = $(round(maximum(angles); sigdigits = 4)) deg",
    )

    input_error = abs.(comparison.cartesian_projection * comparison.ks.B .- comparison.cartesian.B)
    error_plot = heatmap(
        input_error;
        xlabel = "laser input",
        ylabel = "Cartesian state row",
        colorbar_title = "abs error",
        title = "Direct input projection error\nmax = $(round(maximum(input_error); sigdigits = 4))",
    )

    plt = plot(
        orbit_plot,
        geom_plot,
        rank_plot,
        angle_plot,
        error_plot,
        plot(
            framestyle = :none,
            legend = false,
            xlims = (0, 1),
            ylims = (0, 1),
            annotations = [
                (0.05, 0.78, text("Projected KS rank = $(comparison.projected_ks_cartesian_rank)", 12, :left)),
                (0.05, 0.62, text("Cartesian rank = $(comparison.cartesian_rank)", 12, :left)),
                (0.05, 0.46, text("Max principal angle = $(round(maximum(angles); sigdigits = 4)) deg", 12, :left)),
                (0.05, 0.30, text("Max input projection error = $(round(maximum(input_error); sigdigits = 4))", 12, :left)),
            ],
            title = "Numerical agreement",
        );
        layout = (2, 3),
        size = (1900, 1100),
        plot_title = "Laser Controllability: KS Projection Matches Cartesian Test",
    )

    project_root = normpath(joinpath(@__DIR__, ".."))
    output_dir_resolved = isabspath(output_dir) ? output_dir : normpath(joinpath(project_root, output_dir))
    mkpath(output_dir_resolved)
    plot_path = joinpath(output_dir_resolved, "ks_laser_controllability_visualization.png")
    savefig(plt, plot_path)

    println("Saved KS laser controllability visualization: $plot_path")
    println("Projected KS Cartesian rank: $(comparison.projected_ks_cartesian_rank)")
    println("Cartesian rank: $(comparison.cartesian_rank)")
    println("Max principal angle [deg]: $(maximum(angles))")
    println("Input projection max error: $(maximum(input_error))")

    return (
        comparison = comparison,
        plot_path = plot_path,
        rank_history = (
            ks = ks_rank_history,
            projected_ks = projected_rank_history,
            cartesian = cartesian_rank_history,
        ),
        principal_angles_deg = angles,
        input_projection_error = input_error,
    )
end

function run_ks_laser_rtn_analysis_visualization(;
    output_dir = "outputs/ks",
    output_filename = "ks_laser_rtn_analysis.png",
    anomaly_spacing_deg = 0.0,
    altitude_km = 550.0,
    inclination_deg = 0.0,
    eccentricity = 0.0,
    argp_deg = 0.0,
    raan_deg = 0.0,
    nu0_deg = 0.0,
    altitude_offsets_km = [0.0, 50.0],
    j2 = J2_EARTH,
    atmosphere_top_km = 100.0,
    link_max_range_km = 200.0,
    cr = 2.0,
    laser_power_w = 1.0e12,
    satellite_mass_kg = 300.0,
    rtol = 1e-9,
)
    init = initialize_coplanar_ks_constellation(
        2;
        altitude_km = altitude_km,
        inclination_deg = inclination_deg,
        eccentricity = eccentricity,
        argp_deg = argp_deg,
        raan_deg = raan_deg,
        nu0_deg = nu0_deg,
        anomaly_spacing_deg = anomaly_spacing_deg,
        altitude_offsets_km = altitude_offsets_km,
    )

    rtn = ks_laser_rtn_linearization(
        init.u0;
        j2 = j2,
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = ones(Float64, 2),
        rtol = rtol,
    )
    pair = rtn.comparison.pair_data[1]
    pair.zeta == 1.0 || throw(ArgumentError("RTN visualization requires an active laser link."))

    x = rtn.comparison.cartesian_state
    r1 = x[1:3]
    v1 = x[4:6]
    r2 = x[7:9]
    v2 = x[10:12]
    midpoint = 0.5 .* (r1 .+ r2)
    span_km = max(1.7 * norm(r2 .- r1), 8.0)
    axis_scale = 0.22 * span_km

    geometry_plot = plot(
        [r1[1], r2[1]],
        [r1[2], r2[2]],
        [r1[3], r2[3]];
        linewidth = 3,
        color = :deepskyblue3,
        marker = :circle,
        markersize = 7,
        label = "active laser link",
        xlabel = "x [km]",
        ylabel = "y [km]",
        zlabel = "z [km]",
        xlims = (midpoint[1] - span_km, midpoint[1] + span_km),
        ylims = (midpoint[2] - span_km, midpoint[2] + span_km),
        zlims = (midpoint[3] - span_km, midpoint[3] + span_km),
        camera = (35, 25),
        legend = :topright,
        title = "RTN frames at active laser link\nrange = $(round(pair.distance_km; digits = 3)) km",
    )

    for (sat, r_sat, v_sat) in ((1, r1, v1), (2, r2, v2))
        m = ks_rtn_basis(r_sat, v_sat)
        colors = (:crimson, :darkgreen, :purple)
        labels = ("sat $sat R", "sat $sat T", "sat $sat N")
        for axis in 1:3
            endpoint = r_sat .+ axis_scale .* m[:, axis]
            plot!(
                geometry_plot,
                [r_sat[1], endpoint[1]],
                [r_sat[2], endpoint[2]],
                [r_sat[3], endpoint[3]];
                linewidth = 2.5,
                color = colors[axis],
                label = labels[axis],
            )
        end
    end

    eig_plot = scatter(
        real.(rtn.eigenvalues),
        imag.(rtn.eigenvalues);
        marker = :circle,
        markersize = 5,
        label = false,
        xlabel = "real(lambda)",
        ylabel = "imag(lambda)",
        title = "Eigenvalues of A_RTN",
    )
    hline!(eig_plot, [0.0]; color = :gray70, linestyle = :dash, label = false)
    vline!(eig_plot, [0.0]; color = :gray70, linestyle = :dash, label = false)

    eigvec_plot = heatmap(
        abs.(rtn.eigenvectors);
        xlabel = "eigenvector",
        ylabel = "RTN state row",
        colorbar_title = "abs",
        title = "|eigenvectors(A_RTN)|",
    )

    a_plot = heatmap(
        abs.(rtn.A_RTN);
        xlabel = "state column",
        ylabel = "state row",
        colorbar_title = "abs",
        title = "|A_RTN|",
    )

    b_plot = heatmap(
        abs.(rtn.B_RTN);
        xlabel = "laser input",
        ylabel = "state row",
        colorbar_title = "abs",
        title = "|B_RTN|",
    )

    c_plot = heatmap(
        log10.(abs.(rtn.controllability) .+ eps());
        xlabel = "controllability column",
        ylabel = "state row",
        colorbar_title = "log10(abs)",
        title = "RTN controllability matrix\nrank = $(rtn.controllability_rank)",
    )

    j2_label = iszero(float(j2)) ? "J2 off" : "J2 on"
    plt = plot(
        geometry_plot,
        eig_plot,
        eigvec_plot,
        a_plot,
        b_plot,
        c_plot;
        layout = (2, 3),
        size = (1900, 1100),
        plot_title = "Laser Link RTN Linearization And Controllability ($j2_label)",
    )

    project_root = normpath(joinpath(@__DIR__, ".."))
    output_dir_resolved = isabspath(output_dir) ? output_dir : normpath(joinpath(project_root, output_dir))
    mkpath(output_dir_resolved)
    plot_path = joinpath(output_dir_resolved, output_filename)
    savefig(plt, plot_path)

    println("Saved KS laser RTN analysis visualization: $plot_path")
    println("A_RTN size: $(size(rtn.A_RTN))")
    println("B_RTN size: $(size(rtn.B_RTN))")
    println("RTN controllability matrix size: $(size(rtn.controllability))")
    println("RTN controllability rank: $(rtn.controllability_rank)")
    println("Eigenvalue real range: $(extrema(real.(rtn.eigenvalues)))")
    println("Eigenvalue imag range: $(extrema(imag.(rtn.eigenvalues)))")

    return (
        rtn = rtn,
        plot_path = plot_path,
    )
end

function run_ks_two_satellite_linearization_demo(;
    output_dir = "outputs/ks",
    dt_s = 30.0,
    num_orbits = 0.5,
    anomaly_spacing_deg = 0.5,
    perturb_satellite = 2,
    perturb_true_anomaly_deg = 0.01,
    altitude_km = 550.0,
    inclination_deg = 53.0,
    eccentricity = 0.0,
    argp_deg = 0.0,
    raan_deg = 0.0,
    nu0_deg = 0.0,
    link_max_range_km = 200.0,
    atmosphere_top_km = 100.0,
    cr = 2.0,
    laser_power_w = 1.0e4,
    satellite_mass_kg = 300.0,
)
    N = 2
    1 <= perturb_satellite <= N || throw(ArgumentError("perturb_satellite must be 1 or 2 for the demo."))

    reference_init = initialize_coplanar_ks_constellation(
        N;
        altitude_km = altitude_km,
        inclination_deg = inclination_deg,
        eccentricity = eccentricity,
        argp_deg = argp_deg,
        raan_deg = raan_deg,
        nu0_deg = nu0_deg,
        anomaly_spacing_deg = anomaly_spacing_deg,
    )

    perturb_offsets = zeros(Float64, N)
    perturb_offsets[perturb_satellite] = float(perturb_true_anomaly_deg)
    perturbed_init = initialize_coplanar_ks_constellation(
        N;
        altitude_km = altitude_km,
        inclination_deg = inclination_deg,
        eccentricity = eccentricity,
        argp_deg = argp_deg,
        raan_deg = raan_deg,
        nu0_deg = nu0_deg,
        anomaly_spacing_deg = anomaly_spacing_deg,
        anomaly_offsets_deg = perturb_offsets,
    )

    orbital_period_s = 2pi * sqrt(reference_init.a_km^3 / MU_EARTH_KM3_S2)
    physical_tspan = (0.0, float(num_orbits) * orbital_period_s)
    nominal_s_step = float(dt_s) / reference_init.a_km

    reference_sol = propagate_ks_constellation(
        reference_init.u0,
        physical_tspan;
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        saveat_sundman = nominal_s_step,
    )
    perturbed_sol = propagate_ks_constellation(
        perturbed_init.u0,
        physical_tspan;
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        saveat_sundman = nominal_s_step,
    )

    linearization = compute_ks_linearization_sequences(
        reference_sol;
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
    )

    frozen_prediction = rollout_ks_discrete_linearization(
        reference_sol,
        perturbed_sol;
        relinearize = false,
        linearization = linearization,
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
    )
    relinearized_prediction = rollout_ks_discrete_linearization(
        reference_sol,
        perturbed_sol;
        relinearize = true,
        linearization = linearization,
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
    )

    truth_states = [Float64.(x) for x in perturbed_sol.u]
    frozen_errors = ks_position_error_timeseries(truth_states, frozen_prediction)
    relinearized_errors = ks_position_error_timeseries(truth_states, relinearized_prediction)
    frozen_state_errors = ks_state_error_timeseries(truth_states, frozen_prediction)
    relinearized_state_errors = ks_state_error_timeseries(truth_states, relinearized_prediction)
    frozen_reduced_state_errors = ks_reduced_state_error_timeseries(frozen_state_errors, N)
    relinearized_reduced_state_errors = ks_reduced_state_error_timeseries(relinearized_state_errors, N)

    frozen_rms = ks_rms_position_error(frozen_errors)
    relinearized_rms = ks_rms_position_error(relinearized_errors)
    mean_physical_times = [ks_mean_physical_time(x) for x in reference_sol.u]
    sundman_times = reference_sol.t
    state_labels = ks_reduced_state_component_labels(N)

    project_root = normpath(joinpath(@__DIR__, ".."))
    output_dir_resolved = isabspath(output_dir) ? output_dir : normpath(joinpath(project_root, output_dir))
    mkpath(output_dir_resolved)
    error_plot_path = joinpath(output_dir_resolved, "ks_linearization_error_comparison.png")
    frozen_state_error_plot_path = joinpath(output_dir_resolved, "ks_frozen_state_error_grid.png")
    relinearized_state_error_plot_path = joinpath(output_dir_resolved, "ks_relinearized_state_error_grid.png")

    plt = plot(
        mean_physical_times,
        frozen_rms;
        linewidth = 2,
        label = "Frozen first-step linearization",
        xlabel = "Time [s]",
        ylabel = "RMS position error [km]",
        title = "KS Linearization Error vs Nonlinear Truth (Sundman rollout)",
    )
    plot!(
        plt,
        mean_physical_times,
        relinearized_rms;
        linewidth = 2,
        label = "Relinearized each step",
    )
    savefig(plt, error_plot_path)
    plot_ks_state_error_grid(
        mean_physical_times,
        frozen_reduced_state_errors;
        component_labels = state_labels,
        case_label = "Frozen first-step linearization",
        filename = frozen_state_error_plot_path,
    )
    plot_ks_state_error_grid(
        mean_physical_times,
        relinearized_reduced_state_errors;
        component_labels = state_labels,
        case_label = "Relinearized each step",
        filename = relinearized_state_error_plot_path,
    )

    first_pair = isempty(linearization.pair_data[1]) ? nothing : linearization.pair_data[1][1]
    println("Saved KS error comparison plot: $error_plot_path")
    println("Saved frozen state error grid: $frozen_state_error_plot_path")
    println("Saved relinearized state error grid: $relinearized_state_error_plot_path")
    println("Reference samples: $(length(reference_sol.t))")
    println("Final mean time [s]: $(mean_physical_times[end])")
    println("Final Sundman time: $(sundman_times[end])")
    println("First Ac size: $(size(linearization.A_c[1]))")
    println("First Bc size: $(size(linearization.B_c[1]))")
    println("First Ad size: $(size(linearization.A_d[1]))")
    println("First Bd size: $(size(linearization.B_d[1]))")
    if !isnothing(first_pair)
        println("Initial pair zeta: $(first_pair.zeta), distance [km]: $(round(first_pair.distance_km; digits = 3))")
    end
    println("Max frozen RMS position error [km]: $(maximum(frozen_rms))")
    println("Max relinearized RMS position error [km]: $(maximum(relinearized_rms))")

    return (
        reference_sol = reference_sol,
        perturbed_sol = perturbed_sol,
        linearization = linearization,
        frozen_prediction = frozen_prediction,
        relinearized_prediction = relinearized_prediction,
        frozen_errors = frozen_errors,
        relinearized_errors = relinearized_errors,
        frozen_state_errors = frozen_state_errors,
        relinearized_state_errors = relinearized_state_errors,
        frozen_reduced_state_errors = frozen_reduced_state_errors,
        relinearized_reduced_state_errors = relinearized_reduced_state_errors,
        frozen_rms = frozen_rms,
        relinearized_rms = relinearized_rms,
        error_plot_path = error_plot_path,
        frozen_state_error_plot_path = frozen_state_error_plot_path,
        relinearized_state_error_plot_path = relinearized_state_error_plot_path,
        reference_init = reference_init,
        perturbed_init = perturbed_init,
    )
end
