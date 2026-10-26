# Data

Onboard data storage: the data a spacecraft holds, produces and sends.

| File | Holds |
|---|---|
| `data_components.jl` | the records, owned by `Components`: `AbstractDataSource`, `DataStore` (with its run state: full, empty) and `DataSource` |
| `data_storage.jl` | `DataStorage`: `add_data!`, the production hook `data_production_bps`, and the store's updates at boundaries and step ends |

The engine integrates the stored, produced and lost data (`data_stored_bits`,
`data_produced_bits`, `data_lost_bits`) in a data pass after the right-hand
side (`src/orchestration/simulation/engine/dynamics_rhs.jl`). In a run with a
store, every spacecraft has these states, zero for one without a store. The
solver modes that split the right-hand side (`:gravity_backbone_split`,
`:split_imex`, `:multirate`) refuse a run with a data store. The rates the data
links carry, limited at the stores, come from
`src/models/vehicle/interlinks/terminals.jl`. The store's full and empty flags
switch at events and step ends in
`src/orchestration/simulation/callbacks/data_callbacks.jl`, which also writes
the storage columns.
