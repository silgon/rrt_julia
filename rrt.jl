start = Float64[1, 1]
goal = Float64[9, 9]
ε = 0.4 # meters
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

function transition(p1::Array{Float64,1}, p2::Array{Float64,1})
    if dist(p1, p2)< ε
        return p2
    else
        theta = atan2(p2[2]-p1[2], p2[1]-p1[1])
        return Float64[p1[1] + ε*cos(theta), p1[2] + ε*sin(theta)]
    end
end

# RRT algorithm
iter = 0
while iter < n_iter
    r_v = rand(2).*Float64[10,10] # random vertex
    c_v = id_start # closest vertex

    _, c_v = findmin(map(x->dist(r_v, nodes[x]), g.vertices))

    # create new vertex and edge
    n_v = transition(nodes[c_v], r_v)

    # check if allowed
    if !is_allowed(n_v)
        continue
    end

    push!(nodes, n_v)
    new_id = length(nodes)
    add_vertex!(g, new_id) # start
    add_edge!(g, c_v, new_id)
    iter += 1
end
println("Plotting data")
m_nodes = hcat(nodes...)' # matrix with value of nodes
using PyCall, PyPlot
@pyimport seaborn as sns
sns.set_style("whitegrid")
plt = sns.plt
plt[:figure](1,figsize=(8,8))
plt[:clf]()
for ed in g.edges
    plt[:plot](m_nodes[[ed.source, ed.target],1],
               m_nodes[[ed.source, ed.target],2], "b", linewidth=0.3)
end

_, g_v = findmin(map(x->dist(goal,nodes[x]), g.vertices))

# navigate from goal to beginning
c_v = g_v # current vertex
while true
    ie = in_edges(c_v, g)  # only has one parent (from the algorithm)
    if length(ie)==0
        break
    end
    plt[:plot](m_nodes[[ie[1].source, c_v],1], m_nodes[[ie[1].source, c_v],2], "r")
    c_v = ie[1].source
end
