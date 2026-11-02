# KS splitting work-precision

Compares the J2-only KS splitting propagator (`ks_splitting_propagate`, orders 2
and 4) with KS Tsit5, Cartesian Tsit5 and Cartesian DP8 over ten orbits of an
Earth 300 x 3000 km, i = 51.6 deg orbit.

```bash
julia --project=. --threads=1 benchmarks/studies/ks_splitting_work_precision/work_precision.jl
```

Outputs in `results/`:

- `raw_work_precision.csv`: per method and setting (steps per orbit `N` or
  tolerance), the median wall-clock time, accepted steps, force evaluations
  and maximum position error against a Cartesian Vern9 reference.
- `metadata.csv`: run time, orbit period and the reference error against a
  10x tighter Vern9 solution.

The splitting rows give the observed order directly: doubling `N` divides the
error by 4 (order 2) and 16 (order 4). Timings depend on the machine; compare
methods within one run.
