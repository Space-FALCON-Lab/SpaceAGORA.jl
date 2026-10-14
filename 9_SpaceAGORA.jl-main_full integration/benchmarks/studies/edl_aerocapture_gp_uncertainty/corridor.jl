struct TrajectoryPoint
    dt::DateTime
    elapsed_s::Float64
    lat_deg::Float64
    lon_deg::Float64
    alt_m::Float64
end

struct TrajectoryPair
    aerocapture::Vector{TrajectoryPoint}
    edl::Vector{TrajectoryPoint}
    q_weights::Vector{Float64}
end

struct StudyCase
    case_id::String
    anchor_time::DateTime
    gap_s::Float64
    lat_offset_deg::Float64
    lon_offset_deg::Float64
end

const AEROCAPTURE_ENTRY_ALT_M = 125.0e3
const AEROCAPTURE_EXIT_ALT_M = 60.0e3
const EDL_ENTRY_ALT_M = 125.0e3
const EDL_EXIT_ALT_M = 10.0e3
const SAMPLE_PERIOD_S = 2.0
const AEROCAPTURE_DURATION_S = 260.0
const EDL_DURATION_S = 420.0
const NOMINAL_QMAX_ALT_KM = 35.0
const NOMINAL_QMAX_SIGMA_KM = 11.0

@inline function _lerp(a::Float64, b::Float64, x::Float64)::Float64
    return muladd(x, b - a, a)
end

@inline function _normalize_lon_deg(lon_deg::Float64)::Float64
    lon = mod(lon_deg + 180.0, 360.0) - 180.0
    return lon == -180.0 ? 180.0 : lon
end

function _build_corridor(
    start_dt::DateTime;
    duration_s::Float64,
    alt0_m::Float64,
    alt1_m::Float64,
    lat0_deg::Float64,
    lat1_deg::Float64,
    lon0_deg::Float64,
    lon1_deg::Float64,
    sample_period_s::Float64=SAMPLE_PERIOD_S
)::Vector{TrajectoryPoint}
    n_steps = max(2, floor(Int, duration_s / sample_period_s) + 1)
    pts = Vector{TrajectoryPoint}(undef, n_steps)
    for k in 1:n_steps
        x = (k - 1) / max(n_steps - 1, 1)
        elapsed_s = x * duration_s
        pts[k] = TrajectoryPoint(
            start_dt + Millisecond(round(Int, 1000 * elapsed_s)),
            elapsed_s,
            _lerp(lat0_deg, lat1_deg, x),
            _normalize_lon_deg(_lerp(lon0_deg, lon1_deg, x)),
            _lerp(alt0_m, alt1_m, x)
        )
    end
    return pts
end

function _q_weight_profile(points::Vector{TrajectoryPoint})::Vector{Float64}
    weights = Vector{Float64}(undef, length(points))
    for i in eachindex(points)
        alt_km = points[i].alt_m * 1e-3
        bell = exp(-0.5 * ((alt_km - NOMINAL_QMAX_ALT_KM) / NOMINAL_QMAX_SIGMA_KM)^2)
        weights[i] = 0.25 + bell
    end
    return weights
end

"""
    build_trajectory_pair(case; aerocapture_exit_alt_m, edl_exit_alt_m) -> TrajectoryPair

`aerocapture_exit_alt_m` is worth setting deliberately: it is the lowest
altitude the sensing pass reaches, and therefore the top of the band the
estimator can only extrapolate into. The `60 km` default sits above MERRA-2's
`0.1 mb` ceiling near `64 km`, so a MERRA-2 truth source has almost nothing for
the pass to sense unless this is lowered.
"""
function build_trajectory_pair(
    case::StudyCase;
    aerocapture_entry_alt_m::Float64=AEROCAPTURE_ENTRY_ALT_M,
    aerocapture_exit_alt_m::Float64=AEROCAPTURE_EXIT_ALT_M,
    edl_exit_alt_m::Float64=EDL_EXIT_ALT_M,
)::TrajectoryPair
    aero_start = case.anchor_time
    aero_end = aero_start + Millisecond(round(Int, 1000 * AEROCAPTURE_DURATION_S))
    edl_start = aero_end + Millisecond(round(Int, 1000 * case.gap_s))

    base_lat0 = 17.5
    base_lat1 = 23.0
    base_lon0 = -72.0
    base_lon1 = -57.0

    aerocapture = _build_corridor(
        aero_start;
        duration_s=AEROCAPTURE_DURATION_S,
        alt0_m=aerocapture_entry_alt_m,
        alt1_m=aerocapture_exit_alt_m,
        lat0_deg=base_lat0,
        lat1_deg=base_lat1,
        lon0_deg=base_lon0,
        lon1_deg=base_lon1
    )
    edl = _build_corridor(
        edl_start;
        duration_s=EDL_DURATION_S,
        alt0_m=EDL_ENTRY_ALT_M,
        alt1_m=edl_exit_alt_m,
        lat0_deg=base_lat0 + case.lat_offset_deg,
        lat1_deg=base_lat1 + case.lat_offset_deg,
        lon0_deg=base_lon0 + case.lon_offset_deg,
        lon1_deg=base_lon1 + case.lon_offset_deg
    )
    return TrajectoryPair(aerocapture, edl, _q_weight_profile(edl))
end

function default_study_cases()::Vector{StudyCase}
    anchors = (
        DateTime(2024, 3, 20, 18, 0, 0),
        DateTime(2024, 10, 5, 18, 0, 0),
        DateTime(2025, 3, 21, 18, 0, 0),
    )
    gaps_s = (3600.0, 3.0 * 3600.0, 6.0 * 3600.0)
    offsets = (
        (-0.75, -1.25),
        (0.0, 0.0),
        (0.9, 1.2),
        (1.6, -0.8),
    )

    cases = StudyCase[]
    idx = 1
    for anchor in anchors
        for gap_s in gaps_s
            for (lat_off, lon_off) in offsets
                push!(
                    cases,
                    StudyCase(
                        @sprintf("case_%03d", idx),
                        anchor,
                        gap_s,
                        lat_off,
                        lon_off
                    )
                )
                idx += 1
            end
        end
    end
    return cases
end
