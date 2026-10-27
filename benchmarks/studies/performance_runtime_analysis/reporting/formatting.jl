@inline function _safe_unique_join(values; delimiter::String=",")
    vec = collect(skipmissing(values))
    isempty(vec) && return missing
    return join(sort(unique(string.(vec))), delimiter)
end

@inline function _fmt(v; digits::Int=3)
    if v isa Missing
        return "n/a"
    elseif v isa AbstractFloat
        return isfinite(v) ? string(round(v; digits=digits)) : "n/a"
    else
        return string(v)
    end
end

function _scenario_metric(summary_df::DataFrame, scenario::String, metric::Symbol)
    idx = findfirst(==(scenario), summary_df.scenario)
    idx === nothing && return nothing
    return summary_df[idx, metric]
end
