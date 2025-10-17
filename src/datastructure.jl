"""
    Node(i::Int, x::Int, y::Int, q::Int, t::Int, h::int, v::Int)

A `Node` is a point on graph with index `i`, coordinates `(x,y)`, demand `q`, and
tail node index `t`, head node index `h`, and vehicle index `v` in the CVRP solution.
"""
mutable struct Node
    i::Int              # index
    x::Int              # abcissa
    y::Int              # ordinate
    q::Int              # demand
    t::Int              # tail node index
    h::Int              # head node index
    v::Int              # vehicle index
    Node(i, x, y, q) = new(i, x, y, q, i, i, 0)
end

"""
    Arc(t::Int, h::Int, l::Float64)

An `Arc` is a connection between a tail node with index `t` and a head node with 
index `h`, with a traversal cost `c`.
"""
struct Arc
    t::Int              # tail node index
    h::Int              # head node index
    c::Float64          # cost
end

"""
    Vehicle(i::Int, q::Int, s::Int, e::Int, n::Int, l::Int, x::Float64, y::Float64)

A `Vehicle` is a medium of delivery with index `i`, capacity `q`, and start node index `s`, 
end node index `e`, number of customers served `n`, total load `l` and route centroid 
coordinates `(x, y)` in the CVRP solution.
"""
mutable struct Vehicle
    i::Int              # index
    q::Int              # capacity
    s::Int              # start node index
    e::Int              # end node index
    n::Int              # number of customers served
    l::Int              # total load
    x::Float64          # centroid abcissa
    y::Float64          # centroid ordinate
    Vehicle(i, q) = new(i, q, 1, 1, 0, 0, 0., 0.)
end

"""
    Graph(N::Vector{Node}, A::Matrix{Arc}, V::Vector{Vehicle})

A `Graph` is collection of nodes `N`, arcs `A`, and vehicles `V`.
"""
struct Graph
    N::Vector{Node}     # nodes
    A::Matrix{Arc}      # arcs
    V::Vector{Vehicle}  # vehicles
end

"""
    Solution(G::Graph, c::Float64, p::Float64)

A `Solution` is a CVRP solution for graph `G` with cost `c` and penalty `p`.
"""
mutable struct Solution
    G::Graph            # graph
    c::Float64          # cost
    p::Float64          # penalty
end