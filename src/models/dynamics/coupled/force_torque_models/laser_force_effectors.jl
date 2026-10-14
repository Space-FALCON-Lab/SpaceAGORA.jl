module LaserForceEffectors

using LinearAlgebra
using StaticArrays
using ...InterLinkModels: ConstellationInterLinkModel

export force_on_endpoint, laser_force_on_spacecraft

const SPEED_OF_LIGHT_MPS = 299_792_458.0

function force_on_endpoint(current, other, parameters)
	separation = SVector{3, Float64}(current.pos) - SVector{3, Float64}(other.pos)
	distance = norm(separation)
	distance == 0 && return zero(separation)
	return (parameters.B * parameters.P / SPEED_OF_LIGHT_MPS) * (separation / distance)
end

function laser_force_on_spacecraft(model::ConstellationInterLinkModel, state, satellite::Int)
	force = SVector{3, Float64}(0.0, 0.0, 0.0)
	1 <= satellite <= length(model.graph.satellite_edges) || return force
	for index in model.graph.satellite_edges[satellite]
		key = map(v -> model.graph.vertices[v], model.graph.edges[index])
		model.graph.active[index] || continue
		partner = if key[1][1] == satellite
			key[2][1]
		elseif key[2][1] == satellite
			key[1][1]
		else
			continue
		end
		force += force_on_endpoint(state.sc[satellite], state.sc[partner], model.graph.edge_properties[index])
	end
	return force
end

end
