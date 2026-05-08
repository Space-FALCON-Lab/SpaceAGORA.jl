const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_TRAJ_DIR = joinpath(
    REPO_ROOT,
    "output", "performance", "isolated_multi_gravity", "r1_a_outer_only",
    "multi_4_gravity", "pass_01", "r1_a_outer_only", "trajectories",
    "multi_4_gravity", "repeat_001", "attempt_01",
)

using Arrow
using DataFrames
using Plots

function main(argv::Vector{String})
    traj_dir = isempty(argv) ? DEFAULT_TRAJ_DIR : abspath(argv[1])
    isdir(traj_dir) || throw(ArgumentError("Missing trajectory directory: $(traj_dir)"))

    files = sort(collect(joinpath(traj_dir, f) for f in readdir(traj_dir) if occursin(r"^sat_00\d+\.feather$", f)))
    isempty(files) && throw(ArgumentError("No sat_00x.feather files found in $(traj_dir)"))

    p = plot(
        layout=(3, 1),
        size=(1300, 900),
        dpi=150,
        xlabel="time [s]",
        legend=:outerright,
        gridalpha=0.25,
    )

    components = [(:pos_x_m, "x [km]"), (:pos_y_m, "y [km]"), (:pos_z_m, "z [km]")]
    for file in files
        df = DataFrame(Arrow.Table(file))
        sort!(df, :time_s)
        sat_label = "sat $(lpad(string(first(skipmissing(df.sat_id))), 3, '0'))"

        for (subplot, (col, ylabel)) in enumerate(components)
            plot!(
                p[subplot],
                df.time_s,
                df[!, col] ./ 1e3;
                label=sat_label,
                linewidth=1.4,
                ylabel=ylabel,
            )
        end
    end

    out = joinpath(traj_dir, "sat_00x_position_components.png")
    savefig(p, out)
    println("Wrote $(out)")
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
