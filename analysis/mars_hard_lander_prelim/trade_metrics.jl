function compute_authority_trade_metrics(
    config::PrelimConfig,
    summary_df::DataFrame,
    authority_df::DataFrame,
)
    bang_df = filter(row -> row.policy == "bang_bang", summary_df)
    geometry_rows = combine(
        groupby(bang_df, [:target_beta_high_kg_m2, :target_beta_ratio]),
        :mass_kg => first => :mass_kg,
        :panel_area_total_m2 => first => :panel_area_total_m2,
        :panel_area_each_m2 => first => :panel_area_each_m2,
    )

    trade_df = leftjoin(
        authority_df,
        geometry_rows;
        on=[:target_beta_high_kg_m2, :target_beta_ratio],
    )

    added_mass_kg = trade_df.panel_area_total_m2 .* config.panel_system_areal_density_kg_m2
    trade_df[!, :panel_system_areal_density_kg_m2] .= config.panel_system_areal_density_kg_m2
    trade_df[!, :added_panel_system_mass_kg] = added_mass_kg
    trade_df[!, :added_mass_fraction_pct] = 100.0 .* added_mass_kg ./ trade_df.mass_kg
    trade_df[!, :authority_per_added_area_km_per_m2] = [
        area > 0.0 ? auth / area : NaN
        for (auth, area) in zip(trade_df.downrange_authority_km, trade_df.panel_area_total_m2)
    ]
    trade_df[!, :authority_per_added_mass_km_per_kg] = [
        mass > 0.0 ? auth / mass : NaN
        for (auth, mass) in zip(trade_df.downrange_authority_km, trade_df.added_panel_system_mass_kg)
    ]

    baseline_rows = combine(
        groupby(sort(trade_df, :h_jettison_km), [:target_beta_high_kg_m2, :target_beta_ratio]),
        :h_jettison_km => first => :baseline_h_jettison_km,
        :downrange_authority_km => first => :baseline_authority_km,
    )
    trade_df = leftjoin(
        trade_df,
        baseline_rows;
        on=[:target_beta_high_kg_m2, :target_beta_ratio],
    )
    trade_df[!, :authority_loss_vs_latest_jettison_km] = trade_df.baseline_authority_km .- trade_df.downrange_authority_km
    trade_df[!, :authority_retained_vs_latest_jettison_fraction] = [
        base > 0.0 ? auth / base : NaN
        for (auth, base) in zip(trade_df.downrange_authority_km, trade_df.baseline_authority_km)
    ]
    trade_df[!, :authority_loss_vs_latest_jettison_pct] = 100.0 .* (1.0 .- trade_df.authority_retained_vs_latest_jettison_fraction)
    return trade_df
end

function _proposal_ratio_judgement(beta_ratio::Float64)
    if _approx_eq(beta_ratio, 2.0)
        return "efficient low-complexity option"
    elseif _approx_eq(beta_ratio, 4.0)
        return "recommended balance"
    elseif _approx_eq(beta_ratio, 8.0)
        return "aggressive, diminishing returns"
    end
    return "reference / sensitivity case"
end

function _write_markdown_table(io, headers::Vector{String}, rows::Vector{Vector{String}})
    println(io, "| " * join(headers, " | ") * " |")
    println(io, "| " * join(fill("---", length(headers)), " | ") * " |")
    for row in rows
        println(io, "| " * join(row, " | ") * " |")
    end
end

function write_proposal_trade_table(
    config::PrelimConfig,
    trade_df::DataFrame,
    target_guidance_df::DataFrame,
    monte_carlo_summary_df::DataFrame,
)
    representative_trade = sort(
        filter(
            row -> _approx_eq(row.target_beta_high_kg_m2, config.representative_beta_high) &&
                _approx_eq(row.h_jettison_km, config.representative_h_jettison_m / 1e3),
            trade_df,
        ),
        :target_beta_ratio,
    )
    nominal_trade = filter(row -> _approx_eq(row.target_beta_ratio, config.representative_beta_ratio), representative_trade)
    nominal_jettison_trade = sort(
        filter(
            row -> _approx_eq(row.target_beta_high_kg_m2, config.representative_beta_high) &&
                _approx_eq(row.target_beta_ratio, config.representative_beta_ratio),
            trade_df,
        ),
        :h_jettison_km,
    )
    guided = filter(row -> row.policy == "guided_targeted_optimistic", monte_carlo_summary_df)
    body = filter(row -> row.policy == "body_only_monte_carlo", monte_carlo_summary_df)

    open(joinpath(config.output_root, "proposal_trade_table.md"), "w") do io
        println(io, "# Proposal-Facing Trade Table")
        println(io)
        println(io, "## Why Promising")
        why_rows = Vector{Vector{String}}()
        if nrow(nominal_trade) > 0
            row = nominal_trade[1, :]
            push!(why_rows, ["Representative authority", @sprintf("%.2f km", row.downrange_authority_km), "Nominal SHIELD-like case at β_high = $(Int(round(row.target_beta_high_kg_m2))) kg/m², β_ratio = $(Int(round(row.target_beta_ratio))), h_jettison = $(Int(round(row.h_jettison_km))) km."])
        end
        if nrow(target_guidance_df) > 0
            push!(why_rows, ["Target-range solve", @sprintf("%.3f km error", abs(target_guidance_df.range_error_km[1])), "Simple switch-altitude solve can hit a nominal target range closely in the deterministic model."])
        end
        if nrow(guided) > 0 && nrow(body) > 0
            push!(why_rows, ["Downrange scatter reduction", @sprintf("%.2fx", guided.std_reduction_factor_vs_body_only[1]), "Optimistic guided Monte Carlo reduced downrange standard deviation relative to body-only."])
        end
        push!(why_rows, ["Mission niche", "Strong fit", "Targets the gap between passive SHIELD-like rough landers and high-cost precision EDL systems."])
        _write_markdown_table(io, ["Item", "Result", "Interpretation"], why_rows)
        println(io)

        println(io, "## Minimum Viable Architecture")
        mva_rows = [
            ["Entry body", "SHIELD-like passive sphere-cone", "Keep body simple and impact-tolerant."],
            ["Control surfaces", "2 deployable symmetric plates", "Primary function is drag-area modulation, not full lift steering."],
            ["Control logic", "Discrete stow / deploy / jettison", "Avoid continuous flap control to preserve the low-cost story."],
            ["Navigation", "IMU-first", "Coarse targeting and switch timing, not pinpoint precision landing."],
            ["Terminal phase", "Impactor with jettisoned panels", "Do not rely on parachutes or propulsive descent in this proposal scope."],
        ]
        _write_markdown_table(io, ["Element", "Recommended Minimum", "Reason"], mva_rows)
        println(io)

        println(io, "## Mass / Authority Trade")
        println(io, @sprintf("Assumed illustrative panel-system areal density: %.1f kg/m².", config.panel_system_areal_density_kg_m2))
        mass_rows = Vector{Vector{String}}()
        for row in eachrow(filter(r -> !(_approx_eq(r.target_beta_ratio, 1.0)), representative_trade))
            push!(mass_rows, [
                @sprintf("%.0f", row.target_beta_ratio),
                @sprintf("%.3f", row.panel_area_total_m2),
                @sprintf("%.2f", row.added_panel_system_mass_kg),
                @sprintf("%.1f%%", row.added_mass_fraction_pct),
                @sprintf("%.2f", row.downrange_authority_km),
                @sprintf("%.1f", row.authority_per_added_area_km_per_m2),
                @sprintf("%.1f", row.authority_per_added_mass_km_per_kg),
                _proposal_ratio_judgement(row.target_beta_ratio),
            ])
        end
        _write_markdown_table(
            io,
            ["β_ratio", "Added Area (m²)", "Added Mass (kg)", "Mass Increase", "Authority (km)", "Authority / m²", "Authority / kg", "Judgement"],
            mass_rows,
        )
        println(io)

        println(io, "## Early-Jettison Penalty")
        if nrow(nominal_jettison_trade) > 0
            jettison_rows = Vector{Vector{String}}()
            for row in eachrow(nominal_jettison_trade)
                push!(jettison_rows, [
                    @sprintf("%.0f", row.h_jettison_km),
                    @sprintf("%.2f", row.downrange_authority_km),
                    @sprintf("%.2f", row.authority_loss_vs_latest_jettison_km),
                    @sprintf("%.1f%%", row.authority_loss_vs_latest_jettison_pct),
                ])
            end
            _write_markdown_table(
                io,
                ["h_jettison (km)", "Authority (km)", "Loss vs Latest Jettison (km)", "Loss (%)"],
                jettison_rows,
            )
        end
        println(io)

        println(io, "## What Mass / Cost Increase Looks Acceptable")
        acceptable_rows = [
            ["Acceptable", "β_ratio ≈ 4", "Added mass stays modest in the nominal class while still delivering large authority gains."],
            ["Possibly acceptable", "β_ratio ≈ 2", "Very efficient per kg and per area, but less total authority."],
            ["Aggressive", "β_ratio ≈ 8", "Highest authority, but diminishing return per added kg and more panel / mechanism burden."],
        ]
        _write_markdown_table(io, ["Category", "Case", "Interpretation"], acceptable_rows)
        println(io)

        println(io, "## Claims We Should Make")
        should_rows = [
            ["Should make", "Meaningful downrange targeting authority with simple discrete drag-area control."],
            ["Should make", "Proposal-grade evidence for footprint reduction relative to passive body-only entry."],
            ["Should make", "A plausible low-cost architecture exists if actuation stays discrete and limited."],
            ["Should not make", "Strong crossrange steering or precision landing based on the current model."],
            ["Should not make", "Full 6DOF flight-readiness or validated transition-flow performance."],
            ["Should not make", "That the current geometry is exact SHIELD hardware rather than SHIELD-like."],
        ]
        _write_markdown_table(io, ["Claim", "Guidance"], should_rows)
    end
    return nothing
end
