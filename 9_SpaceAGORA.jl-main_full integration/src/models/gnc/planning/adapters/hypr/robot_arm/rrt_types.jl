"""RRT-Connect tree storage for robot-arm joint-space warm starts."""
struct RobotArmRRTConnectTree
    nodes::Vector{Vector{Float64}}
    parents::Vector{Int}
    costs::Vector{Float64}
end

"""RRT-Connect tree storage for robot-arm joint-space warm starts."""
RobotArmRRTConnectTree(root) = RobotArmRRTConnectTree([Float64.(collect(root))], [0], [0.0])
