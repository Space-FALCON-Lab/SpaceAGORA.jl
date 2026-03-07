# GNC Aerobraking Boundary Contract

## Canonical ownership
- `src/gnc/guidance/aerobraking/*` owns aerobraking strategy planning logic (`E-EDG`, `T-EDG`).
- `src/gnc/control/aerobraking/*` owns algorithm-agnostic tracking/execution only.
- `src/mission/operations/aerobraking_policy/*` owns strategy selection policy interfaces.

## Required boundaries
- Guidance files must not include control source files directly.
- Control files must not define targeting/switch-planning strategy solvers.
- Strategy selection must be typed and mission-owned (`AerobrakingStrategyKind`, selector contract).
- DRL integration point is a stub contract only (`DRLPolicyAdapterStub`), with explicit `Not implemented` behavior.

## Canonical interfaces
- `AerobrakingGuidanceInput`
- `AerobrakingGuidanceOutput`
- `AerobrakingControlCommand`
- `AbstractAerobrakingStrategy`
- `AerobrakingPolicyConfig`
- `AbstractAerobrakingPolicySelector`
- `AerobrakingStrategyKind` (`E_EDG`, `T_EDG`)
