module SimulationModel

# --- Top-Level Dependencies ---
using LinearAlgebra
using StaticArrays
using CSV
using DataFrames
using Reexport

# --- Utils ---
# Assuming you have this util file
include("../utils/quaternion_utils.jl")

# --- Submodules ---
# We include the files, which define their own modules.
# We then @reexport their public APIs.

# 1. Core abstract types
include("abstract_types.jl")
@reexport using .AbstractTypes

# 2. Simple hardware data structs
include("components.jl")
@reexport using .Components

# 3. Main container structs (Link, Joint, Model)
include("model.jl")
@reexport using .Model

# 4. Functions for building the model (add_...!)
include("assembly.jl")
@reexport using .Assembly

# 5. Functions for rotations and frames
include("kinematics.jl")
@reexport using .Kinematics

# 6. Functions for calculating effectors (thrusters, etc.)
include("effectors.jl")
@reexport using .Effectors

# 7. Functions for analysis (get_COM, get_inertia, etc.)
include("analysis.jl")
@reexport using .Analysis

end # module SpacecraftSim