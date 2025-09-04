function randomcustomer!(rng::AbstractRNG, k::Int, s::Solution)
    W = [isone(i) ? 0 : 1 for i ∈ eachindex(s.N)]
    I = sample(rng, k, eachindex(s.N), Weights(W))
    for i ∈ I
        n = s.N[i]
        t = s.N[n.t]
        h = s.N[n.h]
        v = s.V[n.v]
        removecustomer!(s, n, t, h, v)
    end
    return s
end

function relatedcustomer!(rng::AbstractRNG, k::Int, s::Solution)
    j = rand(rng, 2:lastindex(s.N))
    W = [isone(i) ? 0. : relatedness(s.N[i], s.N[j]) for i ∈ eachindex(s.N)]
    I = sample(rng, k, eachindex(s.N), Weights(W))
    for i ∈ I
        n = s.N[i]
        t = s.N[n.t]
        h = s.N[n.h]
        v = s.V[n.v]
        removecustomer!(s, n, t, h, v)
    end
    return s
end

function worstcustomer!(rng::AbstractRNG, k::Int, s::Solution)
    W = zeros(Float64, length(s.N))
    z = f(s)
    for i ∈ eachindex(s.N)
        if isone(i) continue end
        n = s.N[i]
        t = s.N[n.t]
        h = s.N[n.h]
        v = s.V[n.v]
        removecustomer!(s, n, t, h, v)
        z′ = f(s)
        W[i] = z - z′
        insertcustomer!(s, n, p, q, v)
    end
    I = sample(rng, k, eachindex(s.N), Weights(W))
    for i ∈ I
        n = s.N[i]
        t = s.N[n.t]
        h = s.N[n.h]
        v = s.V[n.v]
        removecustomer!(s, n, t, h, v)
    end
    return s
end

function randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    W = ones(Float64, length(s.V))
    m = 0
    while true
        j = rand(rng, eachindex(s.V), Weights(W))
        v = s.V[j]
        i = v.s
        for _ ∈ 1:v.n
            n = s.N[i]
            t = s.N[n.t]
            h = s.N[n.h]
            v = s.V[n.v]
            removecustomer!(s, n, t, h, v)
            m += 1
        end
        W[j] = 0.
        if m ≥ k break end
    end
    return s
end


function relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    j = rand(rng, eachindex(s.V))
    W = [relatedness(s.V[i], s.V[j]) for i ∈ eachindex(s.V)]
    m = 0
    while true
        j = rand(rng, eachindex(s.V), Weights(W))
        v = s.V[j]
        i = v.s
        for _ ∈ 1:v.n
            n = s.N[i]
            t = s.N[n.t]
            h = s.N[n.h]
            v = s.V[n.v]
            removecustomer!(s, n, t, h, v)
            m += 1
        end
        W[j] = 0.
        if m ≥ k break end
    end
    return s
end

function worstvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    W = [v.l / v.q for v ∈ s.V]
    m = 0
    while true
        j = rand(rng, eachindex(s.V), Weights(W))
        v = s.V[j]
        i = v.s
        for _ ∈ 1:v.n
            n = s.N[i]
            t = s.N[n.t]
            h = s.N[n.h]
            v = s.V[n.v]
            removecustomer!(s, n, t, h, v)
            m += 1
        end
        W[j] = 0.
        if m ≥ k break end
    end
    return s
end