function randomcustomer!(rng::AbstractRNG, q::Int, s::Solution)
    W = [isone(i) ? 0 : 1 for i ∈ eachindex(s.N)]
    I = sample(rng, q, eachindex(s.N), Weights(W))
    for i ∈ I
        n = s.N[i]
        p = s.N[n.t]
        q = s.N[n.h]
        v = s.V[n.v]
        removecustomer!(s, n, p, q, v)
    end
    return s
end

function relatedcustomer!(rng::AbstractRNG, q::Int, s::Solution)
    k = rand(rng, 2:lastinex(s.N))
    W = [isone(i) ? 0. : relatedness(s.N[i], s.N[k]) for i ∈ eachindex(s.N)]
    I = sample(rng, q, eachindex(s.N), Weights(W))
    for i ∈ I
        n = s.N[i]
        p = s.N[n.t]
        q = s.N[n.h]
        v = s.V[n.v]
        removecustomer!(s, n, p, q, v)
    end
    return s
end

function worstcustomer!(ng::AbstractRNG, q::Int, s::Solution)
    W = zeros(Float64, length(s.N))
    z = f(s)
    for i ∈ eachindex(s.N)
        if isone(i) continue end
        n = s.N[i]
        p = s.N[n.t]
        q = s.N[n.h]
        v = s.V[n.v]
        removecustomer!(s, n, p, q, v)
        z′ = f(s)
        W[i] = z - z′
        insertcustomer!(s, n, p, q, v)
    end
    I = sample(rng, q, eachindex(s.N), Weights(W))
    for i ∈ I
        n = s.N[i]
        p = s.N[n.t]
        q = s.N[n.h]
        v = s.V[n.v]
        removecustomer!(s, n, p, q, v)
    end
    return s
end