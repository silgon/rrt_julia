start = Float64[3, 3, π]
goal = Float64[9, 9, 0]
# start = Float64[3, 3, π/4]
# goal = Float64[8, 2, -π/4]
# ε = 0.4 #
ε_d_max = 0.4
ε_d_min = 0.1
ε_α = pi/12
k_α = 1
n_iter = 20000
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

function boundedAngle(angle)  #[-pi,pi]
    atan2(sin(angle),cos(angle))
end

function dist(p1::Array{Float64,1}, p2::Array{Float64,1})
    sqrt((p1[1]-p2[1])^2+(p1[2]-p2[2])^2+(k_α*boundedAngle(p1[3]-p2[3]))^2)
end

function transition(p1::Array{Float64,1}, p2::Array{Float64,1})
    d = min(ε_d_max, sqrt((p1[1]-p2[1])^2+(p1[2]-p2[2])^2))
    d = max(ε_d_min, d)
    p_12_θ = atan2(p2[2]-p1[2], p2[1]-p1[1])-p1[3]
    α = min(ε_α, abs(p_12_θ))
    s_α = sign(p_12_θ)
    θ = min(ε_α, abs(p2[3]-p1[3]))
    s_θ = sign(p2[3]-p1[3])
    Float64[p1[1]+d*cos(p1[3]+s_α*α),
            p1[2]+d*sin(p1[3]+s_α*α),
            p1[3]+s_θ*θ]
end

# RRT algorithm
iter = 0
while iter < n_iter
    r_v = rand(3).*Float64[10, 10, 2*pi] # random vertex
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
plt = sns.plt
plt[:clf]()
for ed in g.edges
    plt[:plot](m_nodes[[ed.source, ed.target],1],
               m_nodes[[ed.source, ed.target],2], "b", linewidth=0.3)
end

# find closest node to goal
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
