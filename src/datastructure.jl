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

struct Arc
    i::Int              # tail node index
    j::Int              # head node index
    c::Float64          # cost
end

mutable struct Vehicle
    i::Int              # index
    q::Int              # capacity
    s::Int              # start node index
    e::Int              # end node index
    n::Int              # number of customers served
    l::Int              # load
    x::Float64          # centroid abcissa
    y::Float64          # centroid ordinate
    Vehicle(i, q) = new(i, q, 1, 1, 0, 0, 0., 0.)
end

mutable struct Solution
    N::Vector{Node}     # nodes
    A::Matrix{Arc}      # arcs
    V::Vector{Vehicle}  # vehicles
    c::Float64          # cost
    p::Float64          # penalty
    Solution(N, A, V) = new(N, A, V, 0., 0.)
end