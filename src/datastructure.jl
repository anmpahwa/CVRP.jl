abstract type Node end

mutable struct DepotNode <: Node
    i::Int              # index
    x::Float64          # ordinate
    y::Float64          # abcissa
    V::Vector{Vehicle}  # fleet
end


mutable struct CustomerNode <: Node
    i::Int              # index
    x::Float64          # ordinate
    y::Float64          # abcissa
    q::Float64          # demand
    t::Int              # tail node index
    h::Int              # head node index
    v::Int              # vehicle index
end

struct Arc
    i::Int              # tail node index
    j::Int              # head node index
    c::Float64          # cost
end

mutable struct Vehicle
    i::Int              # index
    q::Float64          # capacity
    s::Int              # start node index
    e::Int              # end node index
    n::Int              # number of customer nodes
    x::Float64          # centroid ordinate
    y::Float64          # centroid abcissa
end