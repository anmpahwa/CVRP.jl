function randomcustomer!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # set node weights: uniform
    W = [isone(i) ? 0 : 1 for i ∈ eachindex(G.N)]
    # sample node indices
    I = sample(rng, eachindex(G.N), Weights(W), k, replace=false)
    # loop: remove exactly k sampled nodes
    for i ∈ I
        n = G.N[i]
        t = G.N[n.t]
        h = G.N[n.h]
        v = G.V[n.v]
        removenode!(n, t, h, v, s)
    end
    # return solution
    return s
end

function relatedcustomer!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # randomize a pivot node
    p = rand(rng, 2:lastindex(G.N))
    # set node weights: relatedness
    W = [isone(i) ? 0. : relatedness(G.N[i], G.N[p]) for i ∈ eachindex(G.N)]
    # sample node indices
    I = sample(rng, eachindex(G.N), Weights(W), k, replace=false)
    # loop: remove exactly k sampled nodes
    for i ∈ I
        n = G.N[i]
        t = G.N[n.t]
        h = G.N[n.h]
        v = G.V[n.v]
        removenode!(n, t, h, v, s)
    end
    # return solution
    return s
end

function worstcustomer!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # set node weights: cost
    W = zeros(Float64, length(G.N))
    z = f(s)
    for i ∈ eachindex(G.N)
        if isone(i) continue end
        n = G.N[i]
        t = G.N[n.t]
        h = G.N[n.h]
        v = G.V[n.v]
        removenode!(n, t, h, v, s)
        z′ = f(s)
        W[i] = z - z′
        insertnode!(n, t, h, v, s)
    end
    # sample node indices
    I = sample(rng, eachindex(G.N), Weights(W), k, replace=false)
    # loop: remove exactly k sampled nodes
    for i ∈ I
        n = G.N[i]
        t = G.N[n.t]
        h = G.N[n.h]
        v = G.V[n.v]
        removenode!(n, t, h, v, s)
    end
    # return solution
    return s
end
    return s
end

function randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # set vehicle weights: uniform
    W = ones(Float64, length(G.V))
    # loop: until at least k nodes are removed
    m = 0
    while m ≤ k
        # sample a vehicle
        j = sample(rng, eachindex(G.V), Weights(W))
        v = G.V[j]
        # remove all associated nodes
        for _ ∈ 1:v.n
            i = v.s
            n = G.N[i]
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            m += 1
        end
        # update vehicle weight
        W[j] = 0.
    end
    # return solution
    return s
end


function relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # randomize a pivot vehicle
    p = sample(rng, eachindex(G.V))
    # set vehicle weights: relatedness
    W = [relatedness(G.V[i], G.V[p]) for i ∈ eachindex(G.V)]
    # loop: until at least k nodes are removed
    m = 0
    while m ≤ k
        # sample a vehicle
        j = sample(rng, eachindex(G.V), Weights(W))
        v = G.V[j]
        # remove all associated nodes
        for _ ∈ 1:v.n
            i = v.s
            n = G.N[i]
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            m += 1
        end
        # update vehicle weight
        W[j] = 0.
    end
    # return solution
    return s
end

function worstvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    # set vehicle weights: utilization
    W = [v.l / v.q for v ∈ G.V]
    # loop: until at least k nodes are removed
    m = 0
    while m ≤ k
        # sample a vehicle
        j = sample(rng, eachindex(G.V), Weights(W))
        v = G.V[j]
        # remove all associated nodes
        for _ ∈ 1:v.n
            i = v.s
            n = G.N[i]
            t = G.N[n.t]
            h = G.N[n.h]
            v = G.V[n.v]
            removenode!(n, t, h, v, s)
            m += 1
        end
        # update vehicle weight
        W[j] = 0.
    end
    # return solution
    return s
end