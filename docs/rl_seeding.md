# RL Satellite Seeding

This module implements a DeepSet-based reinforcement learning agent to replace stochastic seeding in SpaceAGORA's constellation design pipeline.

## Overview

The RL seeding module uses a dual DeepSet architecture with PPO (Proximal Policy Optimization) to learn generalizable satellite placement strategies across multiple debris clusters.

## Architecture

### Dual DeepSet Policy Network

The policy network consists of:
- **Satellite encoder**: Per-satellite φ (6→hidden) → global ψ (hidden→embed)
- **Client encoder**: Per-client φ (6→hidden) → global ψ (hidden→embed)
- **Cross-attention**: Combines satellite and client embeddings
- **Policy head**: Outputs 6 orbital elements normalized to [-1, 1]

All layer sizes are configurable via hyperparameters.

### Environment

The `SatelliteSeedingEnv` implements:
- **Action space**: 6D continuous (orbital elements normalized to [-1, 1])
- **Observation space**: Satellite orbitals, client trajectories, orbital bounds, budget, certificate state
- **Reward**: Cost improvement + constraint violation penalties + feasibility bonus

### Cost Function

The cost function re-implements the CAPO certificate audit in SpaceAGORA:
- `unsafe_cost`: Cost from unsafe constraint violations
- `safe_cost`: Cost from safe constraint violations
- `pred_cost`: Cost from prediction errors
- `total_deficit`: Sum of all deficits
- `feasible`: Boolean indicating constellation feasibility

## Installation

The RL dependencies are already included in SpaceAGORA.jl/Project.toml:
- ReinforcementLearningCore
- ReinforcementLearningZoo
- DistributedReinforcementLearning
- Flux
- CUDA (optional, for GPU support)

To install dependencies:
```bash
cd SpaceAGORA.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## Data Setup

Debris cluster data should be placed in `data/debris_clusters/` as CSV files with the format:
```
a,e,inc,raan,arg_p,ta
7.126313709612876e6,0.0021,1.2925978253440085,2.5360785948783926,5.549105059778527,0.7360992321917623
...
```

Files should be named `clients_cluster_<id>.csv` for individual clusters.

## Training

### Basic Training

Run the training example:
```bash
julia --project=. examples/rl_constellation_training.jl
```

This will:
1. Load all debris clusters from `data/debris_clusters/`
2. Configure PPO hyperparameters
3. Setup TensorBoard logging
4. Train the policy network
5. Save checkpoints to `data/rl_models/`

### TensorBoard Monitoring

View training progress with TensorBoard:
```bash
tensorboard --logdir=data/rl_models/tensorboard
```

Logged metrics include:
- Episode rewards and lengths
- Policy and value losses
- Cost components (total_deficit, unsafe_cost, safe_cost, pred_cost)
- Satellite embeddings (φ output)
- Client embeddings (φ output)
- Cross-attention output
- Constellation state (n_sats, feasible)

### Configuration

Training is configured via the `PPOConfig` struct:

```julia
ppo_config = PPOConfig(
    # PPO parameters
    clip_range = 0.2,
    epochs = 10,
    batch_size = 64,
    gamma = 0.99,
    gae_lambda = 0.95,
    
    # Network parameters
    hidden_size = 32,
    embed_size = 16,
    learning_rate = 3e-4,
    
    # Training parameters
    max_episodes = 10000,
    max_steps_per_episode = 100,
    update_frequency = 2048,
    
    # Parallel training
    num_workers = 32,
    num_envs_per_worker = 4,
    
    # Logging
    tensorboard_log_dir = "data/rl_models/tensorboard",
    log_frequency = 10,
    
    # Reward weights
    unsafe_weight = 100.0,
    safe_weight = 10.0,
    pred_weight = 5.0,
    feasibility_threshold = 1e-6,
)
```

### GPU Support

The training automatically detects GPU availability:
- If CUDA is available, uses GPU for training
- If CUDA is unavailable, falls back to CPU with a warning

## Inference

### Running Inference

Use a trained model for satellite seeding:
```bash
julia --project=. examples/rl_constellation_inference.jl
```

This will:
1. Load a debris cluster
2. Load the trained policy
3. Run RL-based seeding
4. Compare with stochastic greedy baseline
5. Display cost comparison

### Programmatic Usage

```julia
using SpaceAGORA.ConstellationDesign

# Load cluster
config = load_cluster_from_csv("data/debris_clusters/clients_cluster_1.csv")

# Configure RL parameters
config["optimizer_params"] = Dict{String,Any}(
    "stage0_method" => "rl",
    "rl_config" => Dict{String,Any}(
        "rl_model_path" => "data/rl_models/latest_model.jld2",
        "max_sats" => 64,
    )
)

# Run RL seeding
result = run_rl_stage0_seeding(config)
constellation = result["constellation_orbitals"]
```

## Integration with Constellation Design

The RL seeding is integrated into the constellation design pipeline via the `stage0_method` configuration:

```yaml
optimizer_params:
  stage0_method: "rl"  # or "stochastic" or "stochastic_greedy"
  rl_config:
    rl_model_path: "data/rl_models/latest_model.jld2"
    max_sats: 64
    max_steps: 100
```

Run the pipeline:
```julia
using SpaceAGORA.ConstellationDesign

config_dict = ingest_yaml("config.yaml")
results = run_debris_controllable_sim(config_dict)
```

## API Reference

### Main Functions

- `load_training_scenarios()` - Load all debris cluster configurations
- `load_cluster_from_csv(path)` - Load a single cluster from CSV
- `load_superset_csv(ids)` - Load and merge multiple clusters
- `train_ppo(config, scenarios)` - Train PPO agent
- `run_rl_stage0_seeding(config)` - Run RL-based seeding
- `run_stochastic_greedy_seeding(config)` - Run stochastic greedy baseline
- `compute_stage0_cost(config, orbitals)` - Compute constellation cost

### Structures

- `DualDeepSetPolicy` - DeepSet policy network
- `SatelliteSeedingEnv` - RL environment
- `RLSatelliteSeedingObservation` - Observation structure
- `PPOConfig` - Training configuration

## Hardware Requirements

- **Minimum**: 32 CPU cores, 1 GPU (8GB+ VRAM) preferred
- **Fallback**: CPU-only training with warning
- **Storage**: 10GB+ for training logs and model checkpoints

## Success Metrics

- **Training convergence**: Reward improvement over episodes
- **Baseline comparison**: Match or exceed stochastic greedy performance
- **Transfer capability**: Zero-shot performance on new debris clusters
- **Scalability**: Training time scales with available cores
- **Feasibility rate**: Percentage of trained policies that produce feasible constellations

## Troubleshooting

### No training scenarios found
Ensure debris cluster CSV files are in `data/debris_clusters/` with the correct naming format.

### GPU not available
The training will fall back to CPU with a warning. For GPU training, install CUDA and the CUDA.jl package.

### Model not found
Ensure you have trained a model first using `rl_constellation_training.jl`. The default model path is `data/rl_models/latest_model.jld2`.

### Poor training convergence
Try adjusting hyperparameters:
- Increase `learning_rate` for faster learning
- Increase `hidden_size` and `embed_size` for more capacity
- Adjust reward weights to balance exploration vs. exploitation
