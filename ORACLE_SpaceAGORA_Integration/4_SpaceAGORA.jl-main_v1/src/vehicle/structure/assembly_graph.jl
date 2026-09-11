"""
    traverse_bodies(model::SpacecraftModel, body::Link)

Traverses the spacecraft model starting from the given body and returns all
bodies connected to it via joints, as well as the index of the root body
for this assembly.
"""
function traverse_bodies(model::SpacecraftModel, body::Link)
    visited = Set{Link}() # Set to keep track of visited bodies
    root_index = 0 # Index of the root body, if needed
    queue = Link[body] # Initialize the queue with the starting body
    push!(visited, body) # Mark the initial body as visited
    while !isempty(queue)
        current_body = popfirst!(queue) # Get the next body in the queue
        if current_body.root
            root_index = findfirst(isequal(current_body), model.roots) # Find the index of the root body
        end
        # Add children bodies to the queue
        for joint in model.joints
            if joint.link1 == current_body && !in(joint.link2, visited)
                push!(queue, joint.link2) # Add the second link of the joint to the queue
                push!(visited, joint.link2) # Mark it as visited
            elseif joint.link2 == current_body && !in(joint.link1, visited)
                push!(queue, joint.link1) # Add the first link of the joint to the queue
                push!(visited, joint.link1) # Mark it as visited
            end
        end
    end
    @assert root_index != 0 "Root body not found in the model" # Ensure a root body was found
    return collect(visited), root_index # Return all visited bodies as a vector
end
