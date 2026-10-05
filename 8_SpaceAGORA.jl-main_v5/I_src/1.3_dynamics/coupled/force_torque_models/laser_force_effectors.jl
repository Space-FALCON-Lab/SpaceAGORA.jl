module LaserForceEffectors

using LinearAlgebra
using StaticArrays
using ...InterLinkModels: InterLinkModel, InterLinkParameters

export force_on_endpoint, laser_force_on_spacecraft

const SPEED_OF_LIGHT_MPS = 299_792_458.0

function force_on_endpoint(current, other, parameters::InterLinkParameters)
	separation = SVector{3, Float64}(current.pos) - SVector{3, Float64}(other.pos)
	distance = norm(separation)
	distance == 0 && return zero(separation)
	return (parameters.B * parameters.P / SPEED_OF_LIGHT_MPS) * (separation / distance)
end

function laser_force_on_spacecraft(model::InterLinkModel, state, satellite::Int)
	force = SVector{3, Float64}(0.0, 0.0, 0.0)
	for (key, connection) in model.linkgraph
		connection.state.active || continue
		partner = if key[1][1] == satellite
			key[2][1]
		elseif key[2][1] == satellite
			key[1][1]
		else
			continue
		end
		force += force_on_endpoint(state.sc[satellite], state.sc[partner], connection.parameters)
	end
	return force
end

end
