using CSV
using DataFrames
using Plots

const CAMPAIGNS = (
    (name="Earth", directory="mpc_earth_campaign",
        drag_limit=500.0, heat_limit=20.0),
    (name="Mars", directory="mpc_mars_campaign",
        drag_limit=12.1, heat_limit=0.15),
    (name="Venus", directory="mpc_venus_campaign",
        drag_limit=7.55, heat_limit=0.29),
)

function campaign_plot(case)
    directory = normpath(joinpath(@__DIR__, "..", "output", case.directory))
    data = CSV.read(joinpath(directory, "simulation_results.csv"), DataFrame)
    time_hours = data.time ./ 3600.0
    phase = data.sc1_mpc_campaign_phase_code

    figure = plot(layout=(2, 2), size=(1250, 800), dpi=180)
    plot!(figure[1], time_hours, data.sc1_mpc_commanded_area_m2,
        xlabel="Elapsed time (h)", ylabel="Commanded area (m²)",
        title="$(case.name) commanded area", label="Area", linewidth=1.6)
    plot!(figure[2], time_hours, data.sc1_mpc_drag_n,
        xlabel="Elapsed time (h)", ylabel="Drag force (N)",
        title="Drag constraint", label="Drag", linewidth=1.6)
    hline!(figure[2], [case.drag_limit], label="Limit",
        color=:red, linestyle=:dash)
    plot!(figure[3], time_hours, data.sc1_mpc_heat_rate_w_cm2,
        xlabel="Elapsed time (h)", ylabel="Heat rate (W/cm²)",
        title="Heat-rate constraint", label="Heat rate", linewidth=1.6)
    hline!(figure[3], [case.heat_limit], label="Limit",
        color=:red, linestyle=:dash)
    plot!(figure[4], time_hours, phase,
        xlabel="Elapsed time (h)", ylabel="Campaign phase",
        yticks=([1.0, 2.0], ["MED", "Targeting"]),
        title="MED-to-targeting transition", label="Phase", linewidth=2)

    savefig(figure, joinpath(directory, "campaign_summary.png"))
    savefig(figure, joinpath(directory, "campaign_summary.pdf"))
    summary = DataFrame(
        planet=[case.name],
        samples=[nrow(data)],
        maximum_drag_n=[maximum(filter(isfinite, data.sc1_mpc_drag_n))],
        maximum_heat_rate_w_cm2=[maximum(filter(
            isfinite, data.sc1_mpc_heat_rate_w_cm2))],
        minimum_area_m2=[minimum(data.sc1_mpc_commanded_area_m2)],
        maximum_area_m2=[maximum(data.sc1_mpc_commanded_area_m2)],
        targeting_reached=[maximum(phase) >= 2.0],
    )
    CSV.write(joinpath(directory, "campaign_summary.csv"), summary)
    return summary
end

available_campaigns = filter(CAMPAIGNS) do case
    isfile(joinpath(@__DIR__, "..", "output", case.directory,
        "simulation_results.csv"))
end
isempty(available_campaigns) && error(
    "No campaign results found. Run AGORA_MPC_Aerobraking_Campaign.jl first.")
println(vcat(campaign_plot.(available_campaigns)...))
