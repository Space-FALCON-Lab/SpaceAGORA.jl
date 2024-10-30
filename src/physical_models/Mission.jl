mutable struct Mission
    e::Int64
    d::Int64
    l::Int64
    a::Int64
    planet::Int64
end

mutable struct InitalParameters
    M::Mission
    gm::Int64
    dm::Int64
    wm::Int64
    am::Int64
    tm::Int64
    cm::Int64
    tc::Int64
    mc::Int64
end

function mission_def(mission)
    M = Mission(mission[:e], mission[:d], mission[:l], mission[:a], mission[:planet])

    ip = InitalParameters(M, mission[:gm], mission[:dm], mission[:wm], mission[:am], mission[:tm], mission[:cm], mission[:tc], mission[:mc])

    return ip
end