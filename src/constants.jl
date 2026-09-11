const MU_EARTH_KM3_S2 = 398600.4418
const R_EARTH_KM = 6378.1363
const J2_EARTH = 1.08262668e-3
const OMEGA_EARTH_RAD_S = 7.2921159e-5

const PLANE_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#4e79a7", "#f28e2b",
]

const NATURAL_EARTH_COASTLINE_50M_URL =
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_coastline.geojson"

plane_color(plane_id::Int) = PLANE_COLORS[mod1(plane_id, length(PLANE_COLORS))]
