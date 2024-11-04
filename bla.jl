
module NestedNodes

export NestedNode

mutable struct Node
    value::Int
    next::Int
end

mutable struct Node2
    value::Int
    next::Int
end

mutable struct NestedNode
    Node::Node
    Node2::Node2
end

p = NestedNode(Node(1,2), Node2(3,4))

end
