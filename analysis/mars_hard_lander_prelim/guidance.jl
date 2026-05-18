struct TargetGuidanceCase
    beta_high_kg_m2::Float64
    beta_ratio::Float64
    h_jettison_km::Float64
    min_range_km::Float64
    max_range_km::Float64
    target_range_km::Float64
end

struct TargetGuidanceResult
    success::Bool
    status::Symbol
    target_range_km::Float64
    achieved_range_km::Float64
    range_error_km::Float64
    h_switch_m::Float64
    h_jettison_m::Float64
    bracket_low_m::Float64
    bracket_high_m::Float64
    iterations::Int
    trajectory_result::TrajectoryResult
end

function select_target_guidance_case(authority_df::DataFrame, config::PrelimConfig)
    subset = filter(
        row -> _approx_eq(row.target_beta_high_kg_m2, config.representative_beta_high) &&
            _approx_eq(row.target_beta_ratio, config.representative_beta_ratio) &&
            _approx_eq(row.h_jettison_km, config.representative_h_jettison_m / 1e3),
        authority_df,
    )
    target_row = nrow(subset) > 0 ? subset[1, :] : authority_df[argmax(authority_df.downrange_authority_km), :]
    target_range_km = target_row.min_impact_downrange_km +
        config.target_range_fraction * (target_row.max_impact_downrange_km - target_row.min_impact_downrange_km)
    return TargetGuidanceCase(
        target_row.target_beta_high_kg_m2,
        target_row.target_beta_ratio,
        target_row.h_jettison_km,
        target_row.min_impact_downrange_km,
        target_row.max_impact_downrange_km,
        target_range_km,
    )
end

function _trajectory_for_switch(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
    h_switch_m::Float64,
    h_jettison_m::Float64;
    save_trajectory::Bool=false,
    h0_m::Real=config.h0_m,
    V0_mps::Real=config.V0_mps,
    gamma0_deg::Real=config.gamma0_deg,
    theta0_rad::Real=config.theta0_rad,
    density_scale::Real=1.0,
)
    policy = ControlPolicy("target_guided", :bang_bang, h_switch_m, h_jettison_m)
    return solve_entry_trajectory(
        config,
        adapter,
        aero_case,
        policy;
        save_trajectory=save_trajectory,
        h0_m=h0_m,
        V0_mps=V0_mps,
        gamma0_deg=gamma0_deg,
        theta0_rad=theta0_rad,
        density_scale=density_scale,
    )
end

function _find_root_bracket(h_values_m::Vector{Float64}, f_values_km::Vector{Float64})
    n = length(h_values_m)
    n == length(f_values_km) || throw(ArgumentError("Bracket vectors must have equal length."))
    for i in 1:(n - 1)
        f1 = f_values_km[i]
        f2 = f_values_km[i + 1]
        if f1 == 0.0
            return (i, i)
        elseif f1 * f2 < 0.0
            return (i, i + 1)
        elseif f2 == 0.0
            return (i + 1, i + 1)
        end
    end
    return nothing
end

function find_switch_altitude_for_target_range(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
    h_jettison_m::Float64,
    target_range_km::Float64;
    h_switch_grid_m::AbstractVector{<:Real}=config.h_switch_grid_m,
    h0_m::Real=config.h0_m,
    V0_mps::Real=config.V0_mps,
    gamma0_deg::Real=config.gamma0_deg,
    theta0_rad::Real=config.theta0_rad,
    density_scale::Real=1.0,
    tolerance_km::Real=config.target_range_tolerance_km,
    max_iterations::Integer=config.target_range_max_iterations,
    save_trajectory::Bool=true,
)::TargetGuidanceResult
    h_grid = sort(Float64.(collect(h_switch_grid_m)))
    coarse_results = TargetGuidanceResult[]
    coarse_ranges = Float64[]
    for h_switch_m in h_grid
        traj = _trajectory_for_switch(
            config,
            adapter,
            aero_case,
            h_switch_m,
            h_jettison_m;
            save_trajectory=false,
            h0_m=h0_m,
            V0_mps=V0_mps,
            gamma0_deg=gamma0_deg,
            theta0_rad=theta0_rad,
            density_scale=density_scale,
        )
        push!(coarse_ranges, traj.summary.impact_downrange_km)
    end

    f_grid = coarse_ranges .- Float64(target_range_km)
    nearest_idx = argmin(abs.(f_grid))
    nearest_h_m = h_grid[nearest_idx]
    nearest_result = _trajectory_for_switch(
        config,
        adapter,
        aero_case,
        nearest_h_m,
        h_jettison_m;
        save_trajectory=save_trajectory,
        h0_m=h0_m,
        V0_mps=V0_mps,
        gamma0_deg=gamma0_deg,
        theta0_rad=theta0_rad,
        density_scale=density_scale,
    )
    nearest_error_km = nearest_result.summary.impact_downrange_km - target_range_km
    abs(nearest_error_km) <= tolerance_km && return TargetGuidanceResult(
        true,
        :grid_hit,
        target_range_km,
        nearest_result.summary.impact_downrange_km,
        nearest_error_km,
        nearest_h_m,
        h_jettison_m,
        nearest_h_m,
        nearest_h_m,
        0,
        nearest_result,
    )

    bracket = _find_root_bracket(h_grid, f_grid)
    if bracket === nothing
        status = target_range_km < minimum(coarse_ranges) ? :clamped_low : :clamped_high
        return TargetGuidanceResult(
            false,
            status,
            target_range_km,
            nearest_result.summary.impact_downrange_km,
            nearest_error_km,
            nearest_h_m,
            h_jettison_m,
            h_grid[1],
            h_grid[end],
            0,
            nearest_result,
        )
    end

    i_lo, i_hi = bracket
    if i_lo == i_hi
        return TargetGuidanceResult(
            true,
            :grid_exact,
            target_range_km,
            nearest_result.summary.impact_downrange_km,
            nearest_error_km,
            h_grid[i_lo],
            h_jettison_m,
            h_grid[i_lo],
            h_grid[i_lo],
            0,
            nearest_result,
        )
    end

    h_lo = h_grid[i_lo]
    h_hi = h_grid[i_hi]
    f_lo = f_grid[i_lo]
    f_hi = f_grid[i_hi]
    final_result = nearest_result
    final_h_m = nearest_h_m
    final_error_km = nearest_error_km
    for iter in 1:Int(max_iterations)
        h_mid = 0.5 * (h_lo + h_hi)
        traj_mid = _trajectory_for_switch(
            config,
            adapter,
            aero_case,
            h_mid,
            h_jettison_m;
            save_trajectory=save_trajectory,
            h0_m=h0_m,
            V0_mps=V0_mps,
            gamma0_deg=gamma0_deg,
            theta0_rad=theta0_rad,
            density_scale=density_scale,
        )
        f_mid = traj_mid.summary.impact_downrange_km - target_range_km
        final_result = traj_mid
        final_h_m = h_mid
        final_error_km = f_mid
        abs(f_mid) <= tolerance_km && return TargetGuidanceResult(
            true,
            :converged,
            target_range_km,
            traj_mid.summary.impact_downrange_km,
            f_mid,
            h_mid,
            h_jettison_m,
            h_lo,
            h_hi,
            iter,
            traj_mid,
        )
        if f_lo * f_mid <= 0.0
            h_hi = h_mid
            f_hi = f_mid
        else
            h_lo = h_mid
            f_lo = f_mid
        end
    end

    return TargetGuidanceResult(
        false,
        :max_iterations,
        target_range_km,
        final_result.summary.impact_downrange_km,
        final_error_km,
        final_h_m,
        h_jettison_m,
        h_lo,
        h_hi,
        Int(max_iterations),
        final_result,
    )
end
