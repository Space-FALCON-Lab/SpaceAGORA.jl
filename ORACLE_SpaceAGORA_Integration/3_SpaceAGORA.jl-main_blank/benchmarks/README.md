# Benchmarks Overview

This folder is for performance testing and analysis of the SpaceAGORA simulation.

## What it does

- Measures runtime performance for different execution modes such as serial, auto, threads, and process.
- Runs benchmark scenarios against orbital and simulation workloads.
- Collects timing data for the full benchmark run, per-orbit sweeps, and entry-duration sweeps.
- Saves raw CSV results and summarized CSV reports.
- Generates plots and markdown summaries for analysis and paper-style reporting.
- Compares how different configurations affect performance, especially parallelism and scaling.

## Main entry point

- [scripts/performance_paper_pipeline.jl](scripts/performance_paper_pipeline.jl)

This script defines the benchmark pipeline, parses settings, runs the benchmark cases, records outputs, and writes reports and plots.

## Study files

The [studies](studies) directory contains many focused performance tests, including:

- runtime analysis
- thread scaling
- constellation scaling
- parallelization benchmarks
- accuracy and surrogate comparisons
- telemetry and optimization tuning

## Output

The pipeline writes results to the repo output area under:

- output/performance/paper_pipeline

## In simple terms

This folder is not the main application code. It is the benchmarking and research tooling used to measure speed, scaling, and performance trade-offs for the simulation.
