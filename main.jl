using Plots
using Interpolations

import config



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
println(p)
println(p.Node.value)

fv = [2040, 1780, 1550, 1313, 1065, 850, 660, 495, 359, 238, 151, 115, 81, 55, 35, 19.5, 9.7, 4.3, 1.5, 0]
vf = [16000,15500, 15000, 14500,14000, 13500, 13000, 12500, 12000, 11500, 11000, 10750, 10500, 10250, 10000, 9750, 9500, 9250, 9000, 0]

fn = linear_interpolation(sort(vf), sort(fv)) # check interpolation

gr()
plot(vf, fv)
scatter!(vf, fn.(vf))

as = Dict(:a => 1, :b => 2, :c => 3)

typeof(as)

as = [1, 2, 3.0]

typeof(as)