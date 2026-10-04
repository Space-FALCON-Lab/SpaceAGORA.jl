# SpaceAGORAHYPR

Compatibility entry point for the separately maintained HYPR package. Installing
this package installs HYPR and SpaceAGORA. `using SpaceAGORAHYPR` enables the
SpaceAGORA adapter and preserves the historical public configuration and result
types. The search implementation is owned by HYPR; this package contains no copy.

This is a local extraction candidate. External repository/release installation
will be pinned after independent review; no published release is claimed yet.

## Candidate compatibility

This 0.2.0 shim requires HYPR 0.1.0 and SpaceAGORA HYPRServices contract 1.0.0.
Companion 0.1 contained the old implementation and must not be loaded with this
pair. If a loading conflict occurs, start a fresh Julia process. Published source
resolution and precompiled installation remain release-integration requirements;
the local candidate is not yet a public installation route.
