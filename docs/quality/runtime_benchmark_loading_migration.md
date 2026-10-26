# Runtime benchmark package-route migration evidence

Status: implemented and under review. On 2026-10-26, the owner explicitly
accepted canonical behavior as a new benchmark baseline, with historical results
unchanged and the difference documented.
Published releases and historical benchmark records remain unchanged.

The maintained runtime entrypoint now uses package-loaded SpaceAGORA types and
the canonical GRAM extension. It removes the benchmark's independent model,
engine and native-adapter ownership. Case definitions, loading precedence,
worker project selection, sampling policies and output schemas remain in place.

## Verified scope

- Fresh full entrypoint import attaches the canonical extension and uses exact
  package module identities; the 167 synthetic loading-policy checks also pass.
- 29 real Earth/Mars native checks pass for construction, copying, serialization,
  production worker startup, transported sampling and collection accounting.
- A three-second non-atmospheric J2 benchmark retains saved columns/values and the
  complete returned trajectory exactly between old and new routes.
- Native samples are not numerically identical. At altitude 100 km, latitude
  0.2 rad, longitude 0.3 rad and elapsed time 30 s from 2020-01-01, Mars density
  changes from `3.8320816570617645e-8` to `3.8893379076285284e-8` kg/m³ (+1.4941%).
  Temperature changes from 124.06160849575598 to 124.56604462113526 K.
- Loading the canonical ephemeris hook while retaining the old independent
  adapter reproduces the new Mars sample exactly. Earth density and temperature
  likewise match; Earth wind still differs. The canonical constructor's native
  static-table priming is absent from the old adapter and addresses the known
  first-instance Earth wind issue. A fresh old-adapter process with the canonical hook and a throwaway Earth
  warm-up reproduces both new Earth and Mars samples exactly.

- Short three-second Earth and Mars runs with native atmosphere and drag each
  repeat exactly within their own route (12 assertions for each route). Their
  output schemas and time endpoints match (8 checks), while numerical values
  differ: the Mars final position components differ by up to 0.000521 m and
  velocity components by up to 0.000348 m/s in this bounded case. These short-run
  deltas cannot bound longer campaigns.

The canonical extension also supplies density-history resets and recognizes when
nominal wind sampling is history dependent. The old independent adapter does not
supply those same contracts. Native samples alone therefore do not qualify a full
atmospheric campaign or establish old/new timing comparability.

## Accepted adoption policy

Use the supported canonical route for maintained benchmarks and identify native
benchmark results as a new numerical baseline, recording source revision, wrapper,
native library/data and sampling policy. Keep historical outputs with their
original provenance. Do not rewrite them or call their timings interchangeable
with measurements from the canonical route. The owner accepted this baseline change on 2026-10-26. Required review, CI and
integration checks still apply; the decision does not certify full scientific
or performance acceptance.

Qualification used existing local native assets without editing them. Initial
harness failures from an incomplete asset directory, missing planet/kernel
initialization and undeclared imports in the legacy environment were corrected;
those failures are retained in the cleanup evidence packet. The package route's
native ownership tests do not certify the third-party binary's internal formulas.
