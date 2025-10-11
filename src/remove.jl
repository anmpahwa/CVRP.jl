"""
    remove!([rng::AbstractRNG], k::Int, s::Solution, method::Symbol)

Returns solution removing `k` nodes from solution s using the given `method`.

Available methods include,
- Random Node Removal       : `:randomnode!`
- Random Arc Removal        : `:randomarc!`
- Random Vehicle Removal    : `:randomvehicle!`
- Related Node Removal      : `:relatednode!`
- Related Arc Removal       : `:relatedarc!`
- Related Vehicle Removal   : `:relatedvehicle!`
- Worst Node Removal        : `:worstnode!`
- Worst Arc Removal         : `:worstarc!`
- Worst Vehicle Removal     : `:worstvehicle!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
remove!(rng::AbstractRNG, k::Int, s::Solution, method::Symbol)::Solution = isdefined(CVRP, method) ? getfield(CVRP, method)(rng, k, s) : getfield(Main, method)(rng, k, s)
remove!(k::Int, s::Solution, method::Symbol) = remove!(Random.GLOBAL_RNG, k, s, method)

"""
    randomnode!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing exactly `k` nodes 
selected randomly.
"""
function randomnode!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # node indices
    I = eachindex(G.N)
    # set node weights: uniform
    W = zeros(Int, I)
    for i ∈ I
        n = G.N[i]
        if isdepot(n) continue end
        W[i] = 1
    end
    # loop: remove exactly k sampled nodes
    j = 0
    while j < k
        i = sample(rng, I, Weights(W))
        n = G.N[i]
        t = G.N[n.t]
        h = G.N[n.h]
        v = G.V[n.v]
        removenode!(n, t, h, v, s)
        W[i] = 0
        j += 1
    end
    # return solution
    return s
end

"""
    relatednode!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing exactly `k` nodes 
selected based on relatedness to a pivot node.
"""
function relatednode!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # node indices
    I = eachindex(G.N)
    # randomize a pivot node
    p = G.N[rand(rng, I)]
    # set node weights: relatedness
    W = zeros(Float64, I)
    for i ∈ I
        n = G.N[i]
        if isdepot(n) continue end
        d = abs(n.x - p.x) + abs(n.y - p.y)
        W[i] = 1 / (d + 1e-3)
    end
    # loop: remove exactly k sampled nodes
    j = 0
    while j < k
        i = sample(rng, I, Weights(W))
        n = G.N[i]
        t = G.N[n.t]
        h = G.N[n.h]
        v = G.V[n.v]
        removenode!(n, t, h, v, s)
        W[i] = 0.
        j += 1
    end
    # return solution
    return s
end

"""
    worstnode!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing exactly `k` nodes 
selected based on removal cost.
"""
function worstnode!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # node indices
    I = eachindex(G.N)
    # set node weights: cost
    W = zeros(Float64, I)
    for i ∈ I
        n = G.N[i]
        if isdepot(n) continue end
        t = G.N[n.t]
        h = G.N[n.h]
        W[i] = (G.A[t.i, n.i].c + G.A[n.i, h.i].c) - G.A[t.i, h.i].c
    end
    # loop: remove exactly k sampled nodes
    j = 0
    while j < k
        i = sample(rng, I, Weights(W))
        n = G.N[i]
        t = G.N[n.t]
        h = G.N[n.h]
        v = G.V[n.v]
        removenode!(n, t, h, v, s)
        W[i] = 0.
        j += 1
    end
    # return solution
    return s
end

"""
    randomarc!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes 
from arcs selected randomly.
"""
function randomarc!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # arc indices
    C = CartesianIndices(G.A)
    I = eachindex(G.A)
    # set arc weights: uniform
    W = zeros(Int, I)
    for i ∈ I
        a = G.A[C[i]]
        t = G.N[a.t]
        h = G.N[a.h] 
        W[i] = isequal(t.h, h.i) ? 1 : 0
    end
    # loop: until at least k nodes are removed
    j = 0
    while j < k
        # sample an arc
        i = sample(rng, I, Weights(W))
        a = G.A[C[i]]
        # remove tail node
        n = G.N[a.t]
        if iscustomer(n) && isclose(n)
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            j += 1
        end
        # remove head node
        n = G.N[a.h]
        if iscustomer(n) && isclose(n)
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            j += 1
        end
        # update arc weight
        W[i] = 0
    end
    # return solution
    return s
end

"""
    relatedarc!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes f
rom arcs selected based on relatedness to s pivot arc.
"""
function relatedarc!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # arc indices
    C = CartesianIndices(G.A)
    I = eachindex(G.A)
    # randomize a pivot arc
    p = G.A[rand(rng, C)]
    # set arc weights: relatedness
    W = zeros(Float64, I)
    for i ∈ I
        a = G.A[C[i]]
        t = G.N[a.t]
        h = G.N[a.h] 
        d = abs((G.N[a.t].x + G.N[a.h].x) - (G.N[p.t].x + G.N[p.h].x)) + abs((G.N[a.t].y + G.N[a.h].y) - (G.N[p.t].y + G.N[p.h].y))
        W[i] = isequal(t.h, h.i) ? (1 / (d + 1e-3)) : 0.
    end
    # loop: until at least k nodes are removed
    j = 0
    while j < k
        # sample an arc
        i = sample(rng, I, Weights(W))
        a = G.A[C[i]]
        # remove tail node
        n = G.N[a.t]
        if iscustomer(n) && isclose(n)
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            j += 1
        end
        # remove head node
        n = G.N[a.h]
        if iscustomer(n) && isclose(n)
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            j += 1
        end
        # update arc weight
        W[i] = 0.
    end
    # return solution
    return s
end

"""
    worstarc!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes 
from arcs selected based on removal cost.
"""
function worstarc!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # arc indices
    C = CartesianIndices(G.A)
    I = eachindex(G.A)
    # set arc weights: cost
    W = zeros(Float64, I)
    for i ∈ I
        a = G.A[C[i]]
        t = G.N[a.t]
        h = G.N[a.h] 
        W[i] = isequal(t.h, h.i) ? a.c : 0.
    end
    # loop: until at least k nodes are removed
    j = 0
    while j < k
        # sample an arc
        i = sample(rng, I, Weights(W))
        a = G.A[C[i]]
        # remove tail node
        n = G.N[a.t]
        if iscustomer(n) && isclose(n)
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            j += 1
        end
        # remove head node
        n = G.N[a.h]
        if iscustomer(n) && isclose(n)
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            j += 1
        end
        # update arc weight
        W[i] = 0.
    end
    # return solution
    return s
end

"""
    randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes 
from vehicles selected randomly.
"""
function randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # vehicle indices
    I = eachindex(G.V)
    # set vehicle weights: uniform
    W = ones(Int, I)
    # loop: until at least k nodes are removed
    j = 0
    while j < k
        # sample a vehicle
        i = sample(rng, I, Weights(W))
        v = G.V[i]
        # remove all associated nodes
        for _ ∈ 1:v.n
            n = G.N[v.s]
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            j += 1
        end
        # update vehicle weight
        W[i] = 0
    end
    # return solution
    return s
end

"""
    relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes 
from vehicles selected based on relatedness to a pivot vehicle.
"""
function relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # vehicle indices
    I = eachindex(G.V)
    # randomize a pivot vehicle
    p = G.V[rand(rng, I)]
    # set vehicle weights: relatedness
    W = zeros(Float64, I)
    for i ∈ I
        v = G.V[i]
        d = abs(v.x - p.x) + abs(v.y - p.y)
        W[i] = 1 / (d + 1e-3)
    end
    # loop: until at least k nodes are removed
    j = 0
    while j < k
        # sample a vehicle
        i = sample(rng, I, Weights(W))
        v = G.V[i]
        # remove all associated nodes
        for _ ∈ 1:v.n
            n = G.N[v.s]
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            j += 1
        end
        # update vehicle weight
        W[i] = 0.
    end
    # return solution
    return s
end

"""
    worstvehicle!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes 
from vehicles selected based on utilization.
"""
function worstvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # vehicle indices
    I = eachindex(G.V)
    # set vehicle weights: utilization
    W = zeros(Float64, I)
    for i ∈ I
        v = G.V[i]
        W[i] = 1 - v.l / v.q
    end
    # loop: until at least k nodes are removed
    j = 0
    while j < k
        # sample a vehicle
        i = sample(rng, I, Weights(W))
        v = G.V[i]
        # remove all associated nodes
        for _ ∈ 1:v.n
            n = G.N[v.s]
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            j += 1
        end
        # update vehicle weight
        W[i] = 0.
    end
    # return solution
    return s
end