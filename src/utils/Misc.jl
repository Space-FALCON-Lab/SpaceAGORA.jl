using LinearAlgebra

function find_neighbours(point, list)
    dist_list = [abs(point - p) for p in list]
    temp = sort(Set(dist_list))

    index_1 = findfirst(==(temp[1]), dist_list)
    index_2 = findfirst(==(temp[2]), dist_list)

    temp = (point-list[index_1]), (point-list[index_2])
    temp_index = [index_1, index_2]
    index_sorted = sort(1:length(temp), by = x -> temp[x])

    return temp_index[index_sorted[1]], temp_index[index_sorted[2]], temp[index_sorted[1]], temp[index_sorted[2]]
end

function new_periapsis(m, r, v, args)
    Energy = (norm(v)^2 / 2) - (m.planet.μ / norm(r))
    a = - m.planet.μ / (2 * Energy)
    h = cross(r, v)
    e = sqrt(1 + (2 * Energy * dot(h,h) / m.planet.μ^2))
    r_p = a * (1 - e)

    if Bool(args[:print_res])
        println("NEW PERIAPSIS VACUUM ALTITUDE: " * string((r_p - m.planet.Rp_e)*1e-3) * " km")
    end
end