# Ground stations: fixed points on the rotating planet whose antennas are link
# terminals. A station is a graph node with a negative node number: antenna `k`
# of station `s` is the vertex `(-s, k)`. This file holds the station, its
# geometry and its link end's state sample; its registration on an
# `InterlinkModel` (`add_ground_station!`) is in interlink_graph.jl. Included by
# `InterLinkModels` (interlink_graph.jl), after the link kinds.

"""
    GroundStation(; name, lat_deg, lon_deg, alt_m=0.0, min_elevation_deg=10.0, antennas)

A ground station: a fixed point on the run's planet, at geodetic latitude
`lat_deg` and longitude `lon_deg` (degrees, east positive) and `alt_m` (m)
above the planet's ellipsoid, with radio `antennas` (a vector of
`RadioTerminal`s, numbered `1:length(antennas)`). Register it with
`add_ground_station!`, which returns its station number `s`; antenna `k` is
then the link endpoint `(-s, k)`.

- **Availability.** A link between a spacecraft and the station is available
  while the spacecraft is at least `min_elevation_deg` above the station's local
  horizon, the plane normal to its geodetic up direction, and the link kind's
  `link_available` allows it. The mask is the link's line of sight: a
  `DataLink` does not test a station end against the planet.
- **The antennas** steer freely above the mask, so their `boresight` is not
  read. A station has no power budget: its antennas draw nothing and are never
  shed.
- **Its data.** In a run with data links, the station counts the bits it
  receives and sends (`gs{s}_received_bits` and `gs{s}_sent_bits`). With data
  stores, it is a store of unbounded capacity and supply: it never limits a
  link.
- **Its motion.** The station turns with the planet's frame (`L_PI`, the
  rotation from J2000 to planet-fixed axes) and spin `ω`: its inertial
  position is `L_PIᵀ r` and its velocity `L_PIᵀ (ω × r)`, with `r` its
  planet-fixed position. Its state sample has no mass and no attitude.

`lat_deg` must be within [-90, 90], `lon_deg` and `alt_m` finite,
`min_elevation_deg` within [0, 90), and `antennas` nonempty.
"""
struct GroundStation
    name::String
    lat_deg::Float64
    lon_deg::Float64
    alt_m::Float64
    min_elevation_deg::Float64
    antennas::Vector{RadioTerminal}

    function GroundStation(name::AbstractString, lat_deg::Real, lon_deg::Real, alt_m::Real, min_elevation_deg::Real,
                           antennas::AbstractVector)
        isempty(name) && throw(ArgumentError("A GroundStation needs a name."))
        isfinite(lat_deg) && -90 <= lat_deg <= 90 || throw(ArgumentError(
            "GroundStation $(name): lat_deg must be within [-90, 90], got $(lat_deg)."))
        isfinite(lon_deg) || throw(ArgumentError("GroundStation $(name): lon_deg must be finite, got $(lon_deg)."))
        isfinite(alt_m) || throw(ArgumentError("GroundStation $(name): alt_m must be finite, got $(alt_m)."))
        isfinite(min_elevation_deg) && 0 <= min_elevation_deg < 90 || throw(ArgumentError(
            "GroundStation $(name): min_elevation_deg must be within [0, 90), got $(min_elevation_deg)."))
        isempty(antennas) && throw(ArgumentError("GroundStation $(name) needs at least one antenna (a RadioTerminal)."))
        for antenna in antennas
            antenna isa RadioTerminal || throw(ArgumentError(
                "GroundStation $(name): its antennas are RadioTerminals, got a $(nameof(typeof(antenna)))."))
        end
        new(String(name), Float64(lat_deg), Float64(lon_deg), Float64(alt_m), Float64(min_elevation_deg),
            Vector{RadioTerminal}(antennas))
    end
end
GroundStation(; name::AbstractString, lat_deg::Real, lon_deg::Real, alt_m::Real=0.0, min_elevation_deg::Real=10.0,
              antennas::AbstractVector) = GroundStation(name, lat_deg, lon_deg, alt_m, min_elevation_deg, antennas)

# The antenna number of an endpoint's terminal on `station`: a number as given,
# or one of its antennas.
_station_antenna_number(station::GroundStation, terminal::Integer) = Int(terminal)
function _station_antenna_number(station::GroundStation, terminal)
    matches = findall(antenna -> antenna === terminal, station.antennas)
    isempty(matches) && throw(ArgumentError(
        "The $(nameof(typeof(terminal))) given as a link endpoint is not an antenna of ground station $(station.name)."))
    length(matches) == 1 || throw(ArgumentError(
        "Ground station $(station.name) has $(length(matches)) antennas equal to the one given (numbers " *
        "$(join(matches, ", "))): give the antenna number."))
    return only(matches)
end

# The station's planet-fixed position (m) on `planet`'s ellipsoid and its
# geodetic up direction.
function _station_site(station::GroundStation, planet)
    lat, lon = deg2rad(station.lat_deg), deg2rad(station.lon_deg)
    position = geodetic_to_pp(lat, lon, station.alt_m, planet)
    up = SVector{3, Float64}(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
    return position, up
end

# The rotation from J2000 to planet-fixed axes that a station end reads: the
# one given (the planet's frame at the state's time, at a right-hand-side
# evaluation), or the planet's current `L_PI`, which the planet-frame callback
# refreshes at every accepted step.
function _station_frame(planet, planet_frame)
    planet === nothing && throw(ArgumentError(
        "A ground station's link end needs the run's planet: pass LinkContext(planet=...)."))
    return SMatrix{3, 3, Float64}(planet_frame === nothing ? planet.L_PI : planet_frame)
end

# The state sample of station number `s`'s link end: its inertial position and
# velocity, no mass and no attitude, the station as its model and `-s` as its
# index.
@noinline function _station_sample(station::GroundStation, s::Int, planet, planet_frame)
    frame = _station_frame(planet, planet_frame)
    position, _ = _station_site(station, planet)
    pos_ii = transpose(frame) * position
    vel_ii = transpose(frame) * cross(SVector{3, Float64}(planet.ω), position)
    return StateSample(pos_ii, vel_ii, NaN, nothing, nothing, station, -s)
end

# The elevation (rad) of the inertial position `position_ii` above the station's
# local horizon: `asin(n̂ · d̂)`, with `n̂` its geodetic up direction and `d̂` the
# planet-fixed direction to the position.
@noinline function _station_elevation_rad(station::GroundStation, planet, planet_frame, position_ii)
    frame = _station_frame(planet, planet_frame)
    position, up = _station_site(station, planet)
    direction = frame * SVector{3, Float64}(position_ii) - position
    distance = norm(direction)
    distance == 0 && return -pi / 2
    return asin(clamp(dot(up, direction) / distance, -1.0, 1.0))
end

# Whether a state sample is a ground station's link end.
_station_end(sample) = sample.spacecraft isa GroundStation

"""
    station_pairs(spacecraft, stations) -> Vector{Tuple{Int, Int}}

Every pair `(i, -s)` of a spacecraft vector index `i` and a ground station
number `s`, spacecraft by spacecraft: a link from every spacecraft to every
station, for `register_links!`. `stations` is the `InterlinkModel` or its
`ground_stations`. For a directed link the spacecraft transmits; for uplinks,
reverse each pair: `reverse.(station_pairs(spacecraft, model))`.
"""
station_pairs(spacecraft, stations::AbstractVector{GroundStation}) =
    [(i, -s) for i in eachindex(spacecraft) for s in eachindex(stations)]
