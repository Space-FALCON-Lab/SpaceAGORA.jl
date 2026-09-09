# Cartesian versus KS Propagator Comparison

## Scope

This benchmark compares two formulations of the same Earth-centered perturbed
two-body problem:

- Cartesian state: `[x,y,z,vx,vy,vz]`, with inverse-square gravity and J2
  evaluated explicitly.
- Kustaanheimo–Stiefel state: `[u₁,u₂,u₃,u₄,u₁′,u₂′,u₃′,u₄′,h_KS,t]`, with
  inverse-square gravity regularized through the energy state and J2 applied as
  a perturbation. The energy convention is `h_KS=-ε` and
  `ω_KS²=h_KS/2`.

Both candidates use classical fixed-step RK4 and exactly the same Earth
constants and analytic J2 acceleration. The Cartesian nominal step is physical
time. The KS nominal fictitious-time step is selected as `Δs = Δt_nominal/a`;
because `dt/ds = r`, KS automatically takes shorter physical-time steps near
perigee and longer ones near apogee. Final KS steps are shortened to terminate
at the same physical final time.

Accuracy is measured against a dense ninth-order `Vern9` Cartesian reference
with `reltol=1e-13` and `abstol=1e-9`. Runtime measurements exclude compilation,
use three repetitions, and report the median complete propagation time,
including state-history storage and KS-to-Cartesian conversion. The benchmark
source is
[`examples/AGORA_KS_Cartesian_Propagator_Comparison.jl`](../../../examples/AGORA_KS_Cartesian_Propagator_Comparison.jl).

## Test cases

| Case | Orbit | Inclination | Propagation interval |
|---|---:|---:|---:|
| Near-circular LEO | `a=7000 km`, `e=0.01` | 51.6° | 10 osculating periods |
| High-eccentricity Earth orbit | `rp=6678.137 km`, `ra=42164 km`, `e≈0.7265` | 63.4° | 3 osculating periods |

The initial point is perigee. Specific mechanical energy includes the J2
potential. The inertial angular-momentum z component, rather than total angular
momentum, is monitored because an axisymmetric J2 potential conserves the
former.

## Results

Measured with Julia 1.12.7 on a 64-bit `znver5` host. Timings are useful for
relative comparison on this host and should not be treated as portable absolute
performance claims.

| Case | Method | Step (s) | RHS calls | Runtime (ms) | Maximum position error (m) | Endpoint position error (m) | Maximum velocity error (m/s) | Relative energy drift | Relative hₙ drift | Allocation (MiB) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEO | Cartesian | 120 | 1944 | 0.190 | 24162.4 | 24162.4 | 26.4218 | 6.415e-5 | 3.196e-5 | 0.988 |
| LEO | KS | 120 | 1960 | 0.775 | 212.867 | 212.867 | 0.232089 | 1.001e-6 | 4.932e-7 | 2.351 |
| LEO | Cartesian | 60 | 3888 | 0.449 | 840.981 | 840.981 | 0.921875 | 2.001e-6 | 9.967e-7 | 1.983 |
| LEO | KS | 60 | 3904 | 1.585 | 9.87645 | 9.87645 | 0.0107691 | 3.131e-8 | 1.542e-8 | 4.688 |
| LEO | Cartesian | 30 | 7772 | 0.922 | 31.7308 | 31.7308 | 0.0349082 | 6.249e-8 | 3.113e-8 | 3.974 |
| LEO | KS | 30 | 7792 | 3.149 | 0.509776 | 0.509776 | 0.000555872 | 1.002e-9 | 4.819e-10 | 9.366 |
| High-e | Cartesian | 120 | 3800 | 0.492 | 79040.4 | 72299.3 | 69.8915 | 7.801e-5 | 4.305e-6 | 1.939 |
| High-e | KS | 120 | 3884 | 1.556 | 1.17672 | 1.08090 | 0.00103756 | 1.026e-9 | 1.292e-11 | 4.663 |
| High-e | Cartesian | 60 | 7600 | 0.997 | 3402.92 | 3116.13 | 3.01931 | 2.925e-6 | 1.352e-7 | 3.888 |
| High-e | KS | 60 | 7744 | 3.235 | 0.0730484 | 0.0671117 | 6.44098e-5 | 6.132e-11 | 4.070e-13 | 9.308 |
| High-e | Cartesian | 30 | 15196 | 1.850 | 164.563 | 150.774 | 0.146284 | 1.195e-7 | 4.220e-9 | 7.787 |
| High-e | KS | 30 | 15468 | 6.641 | 0.00440506 | 0.00404751 | 3.88389e-6 | 3.739e-12 | 2.023e-14 | 18.604 |

The independent tighter-reference check differed from the baseline reference by
`2.07e-6 m` and `2.29e-9 m/s` at the LEO endpoint, and by `1.25e-4 m` and
`8.92e-8 m/s` at the high-e endpoint. These differences are well below every
reported candidate error except the high-e KS 30-second result, where the
reference-position uncertainty is about 2.8% of the reported 4.41 mm error.

## Interpretation

At equal nominal step and nearly equal RHS-call count, KS is substantially more
accurate. At 60 seconds, maximum position error is about 85 times smaller in LEO
and about 46,600 times smaller in the high-eccentricity case. The advantage is
especially large for the high-e orbit because Sundman time transformation
concentrates work near perigee, where Cartesian dynamics vary most rapidly.

The present KS implementation is not faster per accepted step. Its ten-state
algebra, transformations, allocations, and history conversion make this
implementation approximately four to six times slower in these equal-step
tests. It also allocates roughly 2.4–2.6 times more memory. The meaningful
engineering trade is therefore accuracy per force evaluation, not raw cost per
step. In the high-e case, KS at a 120-second nominal step is still far more
accurate than Cartesian RK4 at 30 seconds while requiring about one quarter as
many RHS calls.

KS is most attractive for high-eccentricity motion, repeated close approaches,
long propagation where invariant drift matters, and algorithms needing smooth
regularized dynamics or state sensitivities. Cartesian coordinates remain
simpler for near-circular trajectories, direct coupling to SpaceAGORA force and
attitude states, adaptive production integrators, event handling in physical
time, and models where perturbation evaluation dominates coordinate overhead.

## Limitations and next decisions

This is a formulation benchmark, not yet a comparison of two fully integrated
SpaceAGORA simulation-engine backends. It uses J2 only; drag, third bodies,
spherical harmonics, SRP, thrust, discontinuities, and atmosphere-interface
callbacks require separate studies. Fixed RK4 isolates coordinate effects, but
SpaceAGORA's adaptive Cartesian solver may outperform fixed Cartesian RK4 for a
requested error tolerance. A production KS backend should therefore add an
adaptive fictitious-time solver, allocation-free kernels, physical-time dense
output, and event localization before it replaces the Cartesian backend in a
general simulation.
