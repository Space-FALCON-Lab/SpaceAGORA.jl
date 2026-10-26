# Data

Onboard data storage: the data a spacecraft holds, produces and sends.

| File | Holds |
|---|---|
| `data_components.jl` | the records (owned by `Components`): `AbstractDataSource`, `DataStore` (with its run state: full, empty), `DataSource` |
| `data_storage.jl` | `DataStorage`: `add_data!`, the production hook `data_production_bps`, and the store's boundary and step-end updates |

The engine integrates the stored, produced and lost data
(`data_stored_bits`, `data_produced_bits`, `data_lost_bits`, in the state of
every spacecraft of a run with a store, zero for one without) in its data pass after every right-hand-side route
(`src/orchestration/simulation/engine/dynamics_rhs.jl`). The rates the data
links carry, limited at the stores, come from
`src/models/vehicle/interlinks/terminals.jl`. The store's flags switch at events
and step ends (`src/orchestration/simulation/callbacks/data_callbacks.jl`),
which also writes the storage columns.
