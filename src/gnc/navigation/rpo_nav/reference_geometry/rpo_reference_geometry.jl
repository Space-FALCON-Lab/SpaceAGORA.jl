"""Combined target-station and chaser geometry used by RPO planners."""
struct RPOReferenceGeometry
    station::RPOStationGeometry
    chaser::RPOCubeSatGeometry
end

"""Combined target-station and chaser geometry used by RPO planners."""
function RPOReferenceGeometry(station::RPOStationGeometry; chaser::RPOCubeSatGeometry=RPOCubeSatGeometry())
    return RPOReferenceGeometry(station, chaser)
end

