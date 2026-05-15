struct RPOReferenceGeometry
    station::RPOStationGeometry
    chaser::RPOCubeSatGeometry
end

function RPOReferenceGeometry(station::RPOStationGeometry; chaser::RPOCubeSatGeometry=RPOCubeSatGeometry())
    return RPOReferenceGeometry(station, chaser)
end

