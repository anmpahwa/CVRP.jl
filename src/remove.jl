function randomcustomer!(rng::AbstractRNG, k::Int, s::Solution)
    # set node weights: uniform
    W = [isone(i) ? 0 : 1 for i ∈ eachindex(s.N)]
    # sample node indices
    I = sample(rng, k, eachindex(s.N), Weights(W))
    # loop: remove exactly k sampled nodes
    for i ∈ I
        n = s.N[i]
        t = s.N[n.t]
        h = s.N[n.h]
        v = s.V[n.v]
        removenode!(n, t, h, v, s)
    end
    return s
end

function relatedcustomer!(rng::AbstractRNG, k::Int, s::Solution)
    # randomize a pivot node
    j = rand(rng, 2:lastindex(s.N))
    # set node weights: relatedness
    W = [isone(i) ? 0. : relatedness(s.N[i], s.N[j]) for i ∈ eachindex(s.N)]
    # sample node indices
    I = sample(rng, k, eachindex(s.N), Weights(W))
    # loop: remove exactly k sampled nodes
    for i ∈ I
        n = s.N[i]
        t = s.N[n.t]
        h = s.N[n.h]
        v = s.V[n.v]
        removenode!(n, t, h, v, s)
    end
    return s
end

function worstcustomer!(rng::AbstractRNG, k::Int, s::Solution)
    # set node weights: cost
    W = zeros(Float64, length(s.N))
    z = f(s)
    for i ∈ eachindex(s.N)
        if isone(i) continue end
        n = s.N[i]
        t = s.N[n.t]
        h = s.N[n.h]
        v = s.V[n.v]
        removenode!(n, t, h, v, s)
        z′ = f(s)
        W[i] = z - z′
        insertnode!(n, p, q, v, s)
    end
    # sample node indices
    I = sample(rng, k, eachindex(s.N), Weights(W))
    # loop: remove exactly k sampled nodes
    for i ∈ I
        n = s.N[i]
        t = s.N[n.t]
        h = s.N[n.h]
        v = s.V[n.v]
        removenode!(n, t, h, v, s)
    end
    return s
end

function randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    # set vehicle weights: uniform
    W = ones(Float64, length(s.V))
    # loop: until at least k nodes are removed
    m = 0
    while true
        # sample a vehicle
        j = rand(rng, eachindex(s.V), Weights(W))
        v = s.V[j]
        i = v.s
        # remove all associated nodes
        for _ ∈ 1:v.n
            n = s.N[i]
            t = s.N[n.t]
            h = s.N[n.h]
            v = s.V[n.v]
            removenode!(n, t, h, v, s)
            m += 1
        end
        # update vehicle weight
        W[j] = 0.
        if m ≥ k break end
    end
    return s
end


function relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    # randomize a pivot vehicle
    j = rand(rng, eachindex(s.V))
    # set vehicle weights: relatedness
    W = [relatedness(s.V[i], s.V[j]) for i ∈ eachindex(s.V)]
    # loop: until at least k nodes are removed
    m = 0
    while true
        # sample a vehicle
        j = rand(rng, eachindex(s.V), Weights(W))
        v = s.V[j]
        i = v.s
        # remove all associated nodes
        for _ ∈ 1:v.n
            n = s.N[i]
            t = s.N[n.t]
            h = s.N[n.h]
            v = s.V[n.v]
            removenode!(n, t, h, v, s)
            m += 1
        end
        # update vehicle weight
        W[j] = 0.
        if m ≥ k break end
    end
    return s
end

function worstvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    # set vehicle weights: utilization
    W = [v.l / v.q for v ∈ s.V]
    # loop: until at least k nodes are removed
    m = 0
    while true
        # sample a vehicle
        j = rand(rng, eachindex(s.V), Weights(W))
        v = s.V[j]
        i = v.s
        # remove all associated nodes
        for _ ∈ 1:v.n
            n = s.N[i]
            t = s.N[n.t]
            h = s.N[n.h]
            v = s.V[n.v]
            removenode!(n, t, h, v, s)
            m += 1
        end
        # update vehicle weight
        W[j] = 0.
        if m ≥ k break end
    end
    return s
end