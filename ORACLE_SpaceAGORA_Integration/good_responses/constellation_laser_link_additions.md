# Constellation and Laser-Link Additions

## Main additions and updates

### `constellation.jl`

- **Purpose:** Defines the constellation model.
- **Simple explanation:** Adds the idea of a group of satellites and how they are connected by links.
- **Includes:** `ConstellationPattern`, `WalkerPattern`, `OraclePattern`, `FlowerPattern`, and link helpers such as `activate_link!` and `reset_active_links!`.
- **Plain words:** Models how satellites are arranged and which pairs can communicate or push on each other.

### `laser_link_effectors.jl`

- **Purpose:** Implements the laser force model.
- **Simple explanation:** Models a laser transmitting force from helper satellites to a target spacecraft.
- **Includes:** `OpenCavityLaserLinkModel`, link scheduling, force calculation, and impulse tracking.
- **Plain words:** Determines when helper satellites fire at the target and how much force is applied.

### `force_torque_models.jl`

- **Purpose:** Registers the laser model in the force/torque library.
- **Plain words:** Makes the laser model an available force in the simulation.

### `simulation_model.jl`

- **Purpose:** Exposes the constellation and laser-link types through the main simulation model API.
- **Plain words:** Makes the new features visible and usable by the simulation.

### `SpaceAGORA.jl`

- **Purpose:** Provides the public package entry point.
- **Plain words:** Exports the new constellation and laser-link names so users can call them directly.

## Supporting files

### `ConstellationDesign.md`

Explains the design and structure of the constellation link system.

### `oracle_laser_links.jl`

Example simulation that builds an ORACLE constellation, adds a laser model, and runs the simulation.

### `oracle_laser_plots.jl`

Creates plots for the example results, such as altitude, delta-v, and orbit geometry.

### `output`

Stores generated simulation results, plots, and other runtime artifacts from example runs.
