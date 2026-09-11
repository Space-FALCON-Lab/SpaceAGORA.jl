# Changed files summary: blank baseline vs. main project

## Overview

The folder [3_SpaceAGORA.jl-main_blank](3_SpaceAGORA.jl-main_blank) is the baseline/blank project copy, while [4_SpaceAGORA.jl-main](4_SpaceAGORA.jl-main) is the expanded project version with ORACLE/laser-link functionality added.

I compared the two directories directly. The main differences are the addition of the constellation model, the laser-link force model, and the related example/docs/output files.

## Main files added or changed

### 1. [4_SpaceAGORA.jl-main/src/simulation/constellation.jl](4_SpaceAGORA.jl-main/src/simulation/constellation.jl)

Purpose: defines the constellation model.

Simple explanation: this adds the idea of a group of satellites and how they are connected by links.

It includes:
- `ConstellationPattern`
- `WalkerPattern`
- `OraclePattern`
- `FlowerPattern`
- link activation helpers like `activate_link!`, `deactivate_link!`, and `reset_active_links!`

In plain words: it models how satellites are arranged and which pairs can communicate or push on each other.

### 2. [4_SpaceAGORA.jl-main/src/dynamics/coupled/force_torque_models/laser_link_effectors.jl](4_SpaceAGORA.jl-main/src/dynamics/coupled/force_torque_models/laser_link_effectors.jl)

Purpose: implements the laser force model.

Simple explanation: this is the physics of a laser transmitting force from helper satellites to a target spacecraft.

It adds:
- `OpenCavityLaserLinkModel`
- scheduling logic for when laser links are active
- force calculation and impulse tracking

In plain words: it says when helper satellites point a laser at the target, this is how much force is applied and when.

### 3. [4_SpaceAGORA.jl-main/src/dynamics/coupled/force_torque_models.jl](4_SpaceAGORA.jl-main/src/dynamics/coupled/force_torque_models.jl)

Purpose: registers the new laser model in the force/torque library.

Simple explanation: this file now includes the laser-link effector so the simulation can actually use it.

In plain words: make the new laser model part of the allowed forces in the simulation.

### 4. [4_SpaceAGORA.jl-main/src/core/simulation_model.jl](4_SpaceAGORA.jl-main/src/core/simulation_model.jl)

Purpose: exposes the new constellation and laser-link types through the main model API.

Simple explanation: this file pulls in the new modules and re-exports their public pieces.

In plain words: make the new features visible and usable from the main simulation model.

### 5. [4_SpaceAGORA.jl-main/src/SpaceAGORA.jl](4_SpaceAGORA.jl-main/src/SpaceAGORA.jl)

Purpose: public package entry point.

Simple explanation: it exports the new names so users can call them directly without internal module references.

In plain words: this is the package interface that exposes the laser and constellation features.

## Supporting files

### 6. [4_SpaceAGORA.jl-main/docs/ConstellationDesign.md](4_SpaceAGORA.jl-main/docs/ConstellationDesign.md)

Purpose: explains the design of the constellation model.

Simple explanation: design notes for how the link system is structured.

### 7. [4_SpaceAGORA.jl-main/examples/oracle_laser_links.jl](4_SpaceAGORA.jl-main/examples/oracle_laser_links.jl)

Purpose: example simulation for an ORACLE laser-link scenario.

Simple explanation: a sample run that builds an ORACLE constellation, adds a laser model, and simulates it.

### 8. [4_SpaceAGORA.jl-main/examples/oracle_laser_plots.jl](4_SpaceAGORA.jl-main/examples/oracle_laser_plots.jl)

Purpose: plotting helper for visualizing the example results.

Simple explanation: generates charts for altitude, delta-v, orbit geometry, and related metrics.

### 9. [4_SpaceAGORA.jl-main/output](4_SpaceAGORA.jl-main/output)

Purpose: output folder for generated simulation results and plots.

Simple explanation: runtime artifacts produced by the example runs.

## Bottom line

The main project is the evolved version with orbit-constellation logic plus a laser-force model for ORACLE-style helper satellites, while the blank project is the earlier baseline without those added features.
