"""
    isdepot(n::Node)

Returns `true` if index of node `n` is one, else returns `false`.
"""
@inline isdepot(n::Node) = isone(n.i)

"""
    iscustomer(n::Node)

Returns `false` if index of node `n` is one, else returns `true`.
"""
@inline iscustomer(n::Node) = !isdepot(n)

"""
    isopen(n::Node)

Returns `true` if vehicle index of node `n` is zero, else returns `false`.

Note: open refers to utilization status in the context of a depot node, 
and to service status for a customer node.
"""
@inline isopen(n::Node) = iszero(n.v)

"""
    isclose(n::Node)

Returns `false` if vehicle index of node `n` is zero, else returns `true`.

Note: close refers to utilization status in the context of a depot node, 
and to service status for a customer node.
"""
@inline isclose(n::Node) = !isopen(n)

"""
    isopt(v::Vehicle)

Returns `false` is vehicle `v` serves zero customers, else returns `true`.

Operational refers to utilization status.
"""
@inline isopt(v::Vehicle) = !iszero(v.n)

"""
    vectorize(s::Solution)

Returns a vector of customer node indices in the order of their visit in solution `s`.

Note: visits to/from depot node index is included in the vector.
"""
function vectorize(s::Solution)
    G = s.G
    Z = Int[]
    for v ∈ G.V
        push!(Z, 1)
        i = v.s
        for _ ∈ 1:v.n
            push!(Z, i)
            n = G.N[i]
            i = n.h
        end
        push!(Z, 1)
    end
    return Z
end

"""
    isfeasible(s::Solution)

Returns `true` if solution `s` is feasible, else `false`.

Note: infeasibility refers to capacity violations for CVRP.
"""
@inline function isfeasible(s::Solution)
    G = s.G
    for v ∈ G.V if v.l > v.q return false end end
    return true
end

"""
    f(s::Solution)

Returns the objective function value for solution `s` along with
any penalty for constraint violations, scaled appropriately.
"""
@inline f(s::Solution) = s.c + s.p * 10 ^ ceil(log10(s.c))

"""
    h(s::Solution)

Returns hash code for solution `s`.
"""
@inline h(s::Solution) = hash(vectorize(s))

"""
    Base.deepcopy_internal(G::Graph, dict::IdDict)

Creates a deep copy of Graph `G` using the provided `dict` to track already-copied objects.

Note: This method ensures that all mutable fields of the graph, such as the node set `N` and vehicle set `V`, are recursively deep-copied,
while the adjacency matrix `A` is copied as-is.
"""
@inline Base.deepcopy_internal(G::Graph, dict::IdDict) = Graph(Base.deepcopy_internal(G.N, dict), G.A, Base.deepcopy_internal(G.V, dict))

"""
        Base.deepcopy_internal(s::Solution, dict::IdDict)

Creates a deep copy of solution `s` using the provided `dict` to track already-copied objects.

Note: this method ensures that all mutable fields of the graph, such as the node set `N` and vehicle set `V`, are recursively deep-copied,
while the adjacency matrix `A` is copied as-is.
"""
@inline Base.deepcopy_internal(s::Solution, dict::IdDict) = Solution(Base.deepcopy_internal(s.G, dict), s.c, s.p)