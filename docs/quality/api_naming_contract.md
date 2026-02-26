# SpaceAGORA API Naming Contract

## Scope
This document defines the canonical simulation API naming baseline and file names.

## Canonical Public API Names
The canonical entrypoints are:

1. `execute_case(...)`
2. `execute_campaign(...)`
3. `execute_elements_case(...)`
4. `execute_analysis(...)`
5. `execute_orbital_elements_campaign(...)`
6. `execute_vgamma_campaign(...)`
7. `execute_ae_campaign(...)`

## Canonical File Names
Simulation entry files are:

1. `src/simulation/SimulationExecution.jl`
2. `src/simulation/SimulationElements.jl`

## Contract Rules
1. New public APIs must follow `execute_*` naming.
2. CI contracts and tests must reference canonical file names only.
3. Compatibility aliases are not reintroduced unless explicitly approved and documented.
