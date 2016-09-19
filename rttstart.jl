# start = Float64[0, 0, 0]
# goal = Float64[10, 10, π]
start = Float64[0, 0]
goal = Float64[10, 10]
ε = 0.5 # meters
n_iter = 10000
using Graphs

g = graph(Int[], Graphs.IEdge[], is_directed=true)

# add vertices start and goal
nodes = Array{Float64}[start]
id_start = length(nodes)
add_vertex!(g, id_start) # start

function is_allowed(state::Array{Float64, 1})
    # geometric polygon obstacle [(4, 4), (4,8), (8,8), (8,4)]
    # then out of the boundaries (0,0), (10,10)
    if (4 < state[1] < 8) & (4 < state[2] < 8)
        return false
    elseif (10 < state[1]) | (state[1] < 0) | (10 < state[2]) | (state[2] < 0)
        return false
    end
    return true
end

function dist(p1::Array{Float64,1}, p2::Array{Float64,1})
    sqrt((p1[1]-p2[1])^2+(p1[2]-p2[2])^2)
end

function stop_from_to(p1::Array{Float64,1}, p2::Array{Float64,1})
    if dist(p1, p2)< ε
        return p2
    else
        theta = atan2(p2[2]-p1[2], p2[1]-p1[1])
        return Float64[p1[1] + ε*cos(theta), p1[2] + ε*sin(theta)]
    end
end

# RRT algorithm
for i in 1:n_iter
    r_v = rand(2).*Float64[10,10] # random vertex
    if !is_allowed(r_v)
        continue
    end
    c_v = id_start # closest vertex
    for p in g.vertices
        if dist(r_v, nodes[p]) < dist(nodes[c_v], nodes[p])
            c_v = p
        end
    end

    # create new vertex and edge
    n_v = stop_from_to(nodes[c_v], r_v)
    push!(nodes, n_v)
    new_id = length(nodes)
    add_vertex!(g, new_id) # start
    add_edge!(g, c_v, new_id)
end
