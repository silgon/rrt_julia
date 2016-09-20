####### BEGINNING of COSTMAP ##########
using Distributions

# srand(5)
srand(7)
mvs = MvNormal[]
n_dims = 2
n_gmm = 30
for _ in 1:n_gmm
    tmp = rand(n_dims, n_dims)*10 # range [0, 10]
    Σ = tmp*tmp'
    mv = MvNormal(rand(n_dims)*10, Σ)
    push!(mvs,mv)
end
pre_w = rand(n_gmm)
w = pre_w/sum(pre_w)
m = MixtureModel(mvs, w)
# access the costmap with: pdf(m, [x, y])
##### END OF COSTMAP  ################

max_s = 0.000001 # max sampled

# start = Float64[1, 1]
start = rand(2)*10
ε = 0.4 # meters
n_iter = 2000
using Graphs

g = graph(Int[], Graphs.IEdge[], is_directed=true)

# add vertices start and goal
nodes = Array{Float64}[start]
id_start = length(nodes)
add_vertex!(g, id_start) # start

function is_allowed(state::Array{Float64, 1})
    # then out of the boundaries (0,0), (10,10)
    if (10 < state[1]) | (state[1] < 0) | (10 < state[2]) | (state[2] < 0)
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
iter = 0
while iter < n_iter
    r_v = rand(2).*Float64[10,10] # random vertex

    c_s = pdf(m, r_v)
    if c_s > max_s
        max_s = c_s
    end
    # flip a coin
    if rand() > (c_s/max_s)^3
        continue
    end
    # println(c_s,max_s,c_s/max_s)

    # nearest vertex to random position
    _, c_v = findmin(map(x->dist(r_v, nodes[x]), g.vertices))
    # create new vertex and edge
    n_v = stop_from_to(nodes[c_v], r_v)

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
               m_nodes[[ed.source, ed.target],2], "k", linewidth=0.5)
end



# goal is the maximum value found
_, g_v = findmax(map(x->pdf(m, nodes[x]), g.vertices))

c_v = g_v # current vertex
while true
    ie = in_edges(c_v, g)  # only has one parent (from the algorithm)
    if length(ie)==0
        break
    end
    plt[:plot](m_nodes[[ie[1].source, c_v],1], m_nodes[[ie[1].source, c_v],2], "r")
    c_v = ie[1].source
end



x_n= y_n= 50
x = linspace(0, 10, x_n)
y = linspace(0, 10, y_n)
X = repmat(x, 1, y_n)'
Y = repmat(y, 1, x_n)
Z = zeros(X)
for iter in eachindex(X)
    Z[iter] = pdf(m, [X[iter],Y[iter]])
end
plt[:contour](X, Y, Z, cmap="viridis")
# plt[:pcolormesh](X, Y, Z, cmap="winter")
plt[:colorbar]()
