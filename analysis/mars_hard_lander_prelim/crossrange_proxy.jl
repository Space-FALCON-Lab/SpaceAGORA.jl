function run_crossrange_proxy(
    trajectory_df::DataFrame,
    sigmas::AbstractVector{<:Real},
)
    rows = NamedTuple[]
    n = nrow(trajectory_df)
    if n < 2
        return DataFrame(rows)
    end

    time_s = trajectory_df.time_s
    drag_accel_mps2 = trajectory_df.drag_accel_mps2
    deployed = trajectory_df.deployed
    for sigma in sigmas
        y_m = 0.0
        vy_mps = 0.0
        for i in 1:(n - 1)
            dt = time_s[i + 1] - time_s[i]
            a1 = deployed[i] ? Float64(sigma) * drag_accel_mps2[i] : 0.0
            a2 = deployed[i + 1] ? Float64(sigma) * drag_accel_mps2[i + 1] : 0.0
            a_bar = 0.5 * (a1 + a2)
            y_m += vy_mps * dt + 0.5 * a_bar * dt^2
            vy_mps += a_bar * dt
        end
        push!(rows, (
            sigma=Float64(sigma),
            impact_crossrange_km=y_m / 1e3,
            impact_lateral_velocity_mps=vy_mps,
        ))
    end
    return DataFrame(rows)
end

