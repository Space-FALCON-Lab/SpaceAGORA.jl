# Constellation Design Parity Tests

This directory contains parity tests for comparing SpaceAGORA's constellation design outputs against CAPOConstellations reference outputs.

## Purpose

The parity tests ensure that the ported Stage 0, Stage 1, and Stage 2 functionality in SpaceAGORA produces results that match the original CAPOConstellations implementation within a floating-point tolerance of 1e-2.

## Test Structure

```
test/constellation_parity/
├── README.md                    # This file
├── parity_test_runner.jl        # Main test runner
└── cluster_data/                # Test data directory
    ├── cluster_1/
    │   ├── config.yaml          # YAML configuration for cluster 1
    │   ├── stage0_capo_results.jld2   # CAPO Stage 0 reference results
    │   ├── stage1_capo_results.jld2   # CAPO Stage 1 reference results
    │   └── stage2_capo_results.jld2   # CAPO Stage 2 reference results
    └── cluster_2/
        ├── config.yaml          # YAML configuration for cluster 2
        ├── stage0_capo_results.jld2   # CAPO Stage 0 reference results
        ├── stage1_capo_results.jld2   # CAPO Stage 1 reference results
        └── stage2_capo_results.jld2   # CAPO Stage 2 reference results
```

## Setting Up Test Data

To run the parity tests, you need to:

1. Copy the cluster configuration YAML files from CAPOConstellations to the `cluster_data/cluster_X/` directories
2. Run CAPOConstellations on clusters 1 and 2 to generate reference results
3. Save the reference results as JLD2 files in the appropriate directories

### Generating Reference Results

In CAPOConstellations, run:

```bash
# For cluster 1
julia --project src/debris_controllable_sim.jl config=<path_to_cluster1_config> mode=full

# Save the results from each stage to JLD2 files
```

## Running the Tests

### Run all parity tests for clusters 1 and 2:

```bash
julia --project test/constellation_parity/parity_test_runner.jl
```

### Run tests for specific clusters:

```bash
julia --project test/constellation_parity/parity_test_runner.jl 1
julia --project test/constellation_parity/parity_test_runner.jl 1 2
```

## Test Criteria

- **Floating-point arrays**: Maximum absolute difference < 1e-2
- **Integer/boolean arrays**: Exact match
- **Structured outputs**: Compare key fields with tolerance

## Current Status

The parity test framework is in place. However, the actual CAPO reference results need to be generated and copied from the CAPOConstellations repository before the tests can run successfully.

## Notes

- The current Stage 0, Stage 1, and Stage 2 implementations in SpaceAGORA are simplified stubs that provide the basic structure but do not yet implement the full mathematical models from CAPOConstellations.
- Full parity will require implementing the complete optimization logic, access matrix computation, and control verification algorithms.
- The default formulation (polyhedral set mode with constructive zonotope certificates) is used for comparison.
