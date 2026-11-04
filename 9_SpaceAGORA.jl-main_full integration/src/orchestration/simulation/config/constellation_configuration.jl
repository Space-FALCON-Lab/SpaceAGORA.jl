# The existing spacecraft collection and selected dynamics effectors for one run.
# Included inside SpacecraftModels to preserve DynamicsModel type identity and public access.
# No additional collection wrapper or constellation-generation algorithm is introduced.

"""
Struct containing the information needed for the dynamics propagation of all spacecraft

Roots: Vector of root links (main bus or core bodies). This is the vector over which the threads will be parallelized, so each thread will handle the dynamics propagation of one root link and its associated sub-links.
DynamicEffectors: Tuple of dynamic effector models (gravity, drag, etc.) to be applied during the dynamics propagation.
"""
struct DynamicsModel{T_Effectors<:Tuple}
    spacecraft::Vector{SpacecraftModel} # Vector of root links (main bus or core bodies)
    dynamic_effectors::T_Effectors # Tuple of dynamic effector models (gravity, drag, etc.)

    function DynamicsModel(roots::Vector{SpacecraftModel}, dynamic_effectors::T_Effectors) where {T_Effectors<:Tuple}
        new{T_Effectors}(roots, dynamic_effectors)
    end
end
