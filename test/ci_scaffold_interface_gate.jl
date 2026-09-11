const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

if !isdefined(@__MODULE__, :Resources)
    include(joinpath(REPO_ROOT, "experimental", "vehicle", "resources", "resources.jl"))
end
if !isdefined(@__MODULE__, :ConstellationNetwork)
    include(joinpath(REPO_ROOT, "experimental", "mission", "constellation", "network", "network.jl"))
end
if !isdefined(@__MODULE__, :Estimation)
    include(joinpath(REPO_ROOT, "experimental", "gnc", "estimation", "estimation.jl"))
end
if !isdefined(@__MODULE__, :LaserTerminal)
    include(joinpath(REPO_ROOT, "experimental", "vehicle", "laser_terminal", "laser_terminal.jl"))
end

function _expect_not_implemented(callable::Function, label::String)
    try
        callable()
    catch err
        err isa ErrorException || error("$label did not throw ErrorException: $(typeof(err))")
        msg = String(err.msg)
        startswith(msg, "Not implemented: ") || error("$label threw wrong message: $msg")
        return
    end
    error("$label did not throw")
end

let
    state = Resources.ResourceState()
    inputs = Resources.ResourceInputs()

    battery = Resources.BatteryResourceModel()
    solar = Resources.SolarArrayResourceModel()
    bus = Resources.PowerBusResourceModel()
    load = Resources.LoadResourceModel()

    _expect_not_implemented(() -> Resources.initialize_resource_state(battery), "battery initialize_resource_state")
    _expect_not_implemented(() -> Resources.step_resource!(battery, state, inputs), "battery step_resource!")
    _expect_not_implemented(() -> Resources.initialize_resource_state(solar), "solar initialize_resource_state")
    _expect_not_implemented(() -> Resources.step_resource!(solar, state, inputs), "solar step_resource!")
    _expect_not_implemented(() -> Resources.initialize_resource_state(bus), "power_bus initialize_resource_state")
    _expect_not_implemented(() -> Resources.step_resource!(bus, state, inputs), "power_bus step_resource!")
    _expect_not_implemented(() -> Resources.initialize_resource_state(load), "load initialize_resource_state")
    _expect_not_implemented(() -> Resources.step_resource!(load, state, inputs), "load step_resource!")
end

let
    topology = ConstellationNetwork.NetworkTopology()
    state = ConstellationNetwork.NetworkState()
    inputs = ConstellationNetwork.NetworkInputs()
    _expect_not_implemented(() -> ConstellationNetwork.step_network!(topology, state, inputs), "step_network!")
end

let
    estimator_state = Estimation.EstimatorState()
    measurement = Estimation.MeasurementPacket()
    _expect_not_implemented(() -> Estimation.estimate_state!(estimator_state, measurement), "estimate_state!")
end

let
    actuator = LaserTerminal.LaserTerminalActuator()
    command = LaserTerminal.LaserTerminalCommand()
    _expect_not_implemented(() -> LaserTerminal.apply_laser_terminal_command!(actuator, command, 0.0), "apply_laser_terminal_command!")
end

println("scaffold_interface_gate_ok")
