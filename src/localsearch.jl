"""
    intramove!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after moving a randomly selected node to its best position 
in the same route, if the move results in a reduction in objective function value, 
repeating for `k` iterations.    
"""
function intramove!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # initialize
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample a random node
        n = sample(rng, N, Weights(W))
        # remove the node from its current position
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        # iterate through all positions in the route
        c = 0.
        p = (t.i, h.i, v.i)
        t = N[1]
        h = N[v.s]
        for _ ∈ 0:v.n
            # insert the node
            insertnode!(n, t, h, v, s)
            # evaluate the insertion cost
            z′ = f(s)
            c′ = z′ - z
            # reset the best insertion position for the node
            if c′ < c c, p = c′, (t.i, h.i, v.i) end
            # remove the node
            removenode!(n, t, h, v, s)
            t = h
            h = N[t.h]
        end
        # insert the node to its best positiion
        t = N[p[1]]
        h = N[p[2]]
        v = V[p[3]]
        insertnode!(n, t, h, v, s)
    end
    # return solution
    return s
end

"""
    intermove!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after moving a randomly selected node to its best position 
in a randomly selected route, if the move results in a reduction in objective 
function value, repeating for `k` iterations.    
"""
function intermove!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # initialize
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    R = [zeros(Float64, eachindex(V)) for _ ∈ eachindex(V)]
    for v ∈ V
        for u ∈ V
            if isequal(u,v) continue end
            x = abs(u.x - v.x)
            y = abs(u.y - v.y)
            d = x + y
            r = 1 / (d + 1e-3)
            R[u.i][v.i] = r
        end
    end
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample a random node
        n = sample(rng, N, Weights(W))
        # remove the node from its current position
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        # sample a random vehicle and iterate through all positions in it
        c = 0.
        p = (t.i, h.i, v.i)
        v = sample(rng, V, Weights(R[v.i]))
        t = N[1]
        h = N[v.s]
        for _ ∈ 0:v.n
            # insert the node
            insertnode!(n, t, h, v, s)
            # evaluate the insertion cost
            z′ = f(s)
            c′ = z′ - z
            # reset the best insertion position for the node
            if c′ < c c, p = c′, (t.i, h.i, v.i) end
            # remove the node
            removenode!(n, t, h, v, s)
            t = h
            h = N[t.h]
        end
        # insert the node to its best positiion
        t = N[p[1]]
        h = N[p[2]]
        v = V[p[3]]
        insertnode!(n, t, h, v, s)
    end
    # return solution
    return s
end

"""
    intraswap!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after swapping two randomly selected nodes from 
the same route if the swap results in a reduction in objective function 
value, repeating for `k` iterations.
"""
function intraswap!(rng::AbstractRNG, k::Int, s::Solution)
    # tₘ → m → hₘ and tₙ → n → hₙ
    # pre-initialize
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    # initialize
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    R = [zeros(Float64, eachindex(N)) for _ ∈ eachindex(N)]
    for n ∈ N
        for m ∈ N
            if isequal(n, m) continue end
            if !isequal(n.v, m.v) continue end
            d = A[n.i,m.i].c        
            w = W[n.i] * W[m.i]
            r = w / (d + 1e-3)
            R[m.i][n.i] = r
        end
    end
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample a random node
        n = sample(rng, N, Weights(W))
        # sample another node from the same route, weighted by relatedness
        m = sample(rng, N, Weights(R[n.i]))
        # swap the two nodes
        if isdepot(m) continue end
        if isequal(n,m) continue end
        tₙ = N[n.t]
        hₙ = N[n.h]
        vₙ = V[n.v]
        tₘ = N[m.t]
        hₘ = N[m.h]
        vₘ = V[m.v]
        # case 1: tₘ → m (tₙ) → n (hₘ) → hₙ
        if isequal(n, hₘ) # || isequal(m, tₙ)
            removenode!(n, tₙ, hₙ, vₙ, s)
            insertnode!(n, tₘ, tₙ, vₘ, s)
        # case 2: tₙ → n (tₘ) → m (hₙ) → hₘ
        elseif isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(n, tₙ, hₙ, vₙ, s)
            insertnode!(n, hₙ, hₘ, vₘ, s)
        # case 3: all other cases
        else
            removenode!(n, tₙ, hₙ, vₙ, s)
            removenode!(m, tₘ, hₘ, vₘ, s)
            insertnode!(n, tₘ, hₘ, vₘ, s)
            insertnode!(m, tₙ, hₙ, vₙ, s)
        end
        # evaluate the change in objective function value
        z′ = f(s)
        c  = z′ - z
        if c < 0 continue end
        # reswap the two nodes
        tₙ = N[n.t]
        hₙ = N[n.h]
        vₙ = V[n.v]
        tₘ = N[m.t]
        hₘ = N[m.h]
        vₘ = V[m.v]
        # case 1: tₘ → m (tₙ) → n (hₘ) → hₙ
        if isequal(n, hₘ) # || isequal(m, tₙ)
            removenode!(n, tₙ, hₙ, vₙ, s)
            insertnode!(n, tₘ, tₙ, vₘ, s)
        # case 2: tₙ → n (tₘ) → m (hₙ) → hₘ
        elseif isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(n, tₙ, hₙ, vₙ, s)
            insertnode!(n, hₙ, hₘ, vₘ, s)
        # case 3: all other cases
        else
            removenode!(n, tₙ, hₙ, vₙ, s)
            removenode!(m, tₘ, hₘ, vₘ, s)
            insertnode!(n, tₘ, hₘ, vₘ, s)
            insertnode!(m, tₙ, hₙ, vₙ, s)
        end
    end
    # return solution
    return s
end

"""
    interswap!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after swapping two randomly selected nodes from 
two different routes if the swap results in a reduction in objective 
function value, repeating for `k` iterations.
"""
function interswap!(rng::AbstractRNG, k::Int, s::Solution)
    # tₘ → m → hₘ and tₙ → n → hₙ
    # pre-initialize
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    # initialize
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    R = [zeros(Float64, eachindex(N)) for _ ∈ eachindex(N)]
    for n ∈ N
        for m ∈ N
            if isequal(n, m) continue end
            if isequal(n.v, m.v) continue end
            d = A[n.i,m.i].c        
            w = W[n.i] * W[m.i]
            r = w / (d + 1e-3)
            R[m.i][n.i] = r
        end
    end
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample a random node
        n = sample(rng, N, Weights(W))
        # sample another node from the same route, weighted by relatedness
        m = sample(rng, N, Weights(R[n.i]))
        # swap the two nodes
        if isdepot(m) continue end
        if isequal(n,m) continue end
        tₙ = N[n.t]
        hₙ = N[n.h]
        vₙ = V[n.v]
        tₘ = N[m.t]
        hₘ = N[m.h]
        vₘ = V[m.v]
        # case 1: tₘ → m (tₙ) → n (hₘ) → hₙ
        if isequal(n, hₘ) # || isequal(m, tₙ)
            removenode!(n, tₙ, hₙ, vₙ, s)
            insertnode!(n, tₘ, tₙ, vₘ, s)
        # case 2: tₙ → n (tₘ) → m (hₙ) → hₘ
        elseif isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(n, tₙ, hₙ, vₙ, s)
            insertnode!(n, hₙ, hₘ, vₘ, s)
        # case 3: all other cases
        else
            removenode!(n, tₙ, hₙ, vₙ, s)
            removenode!(m, tₘ, hₘ, vₘ, s)
            insertnode!(n, tₘ, hₘ, vₘ, s)
            insertnode!(m, tₙ, hₙ, vₙ, s)
        end
        # evaluate the change in objective function value
        z′ = f(s)
        c  = z′ - z
        if c < 0 continue end
        # reswap the two nodes
        tₙ = N[n.t]
        hₙ = N[n.h]
        vₙ = V[n.v]
        tₘ = N[m.t]
        hₘ = N[m.h]
        vₘ = V[m.v]
        # case 1: tₘ → m (tₙ) → n (hₘ) → hₙ
        if isequal(n, hₘ) # || isequal(m, tₙ)
            removenode!(n, tₙ, hₙ, vₙ, s)
            insertnode!(n, tₘ, tₙ, vₘ, s)
        # case 2: tₙ → n (tₘ) → m (hₙ) → hₘ
        elseif isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(n, tₙ, hₙ, vₙ, s)
            insertnode!(n, hₙ, hₘ, vₘ, s)
        # case 3: all other cases
        else
            removenode!(n, tₙ, hₙ, vₙ, s)
            removenode!(m, tₘ, hₘ, vₘ, s)
            insertnode!(n, tₘ, hₘ, vₘ, s)
            insertnode!(m, tₙ, hₙ, vₙ, s)
        end
    end
    # return solution
    return s
end

function intraopt!(rng::AbstractRNG, k::Int, s::Solution)
    # initialize
    G = s.G
    N = G.N
    V = G.V
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    R = [[(isequal(n, m) || !isequal(n.v, m.v)) ? 0. : relatedness(n, m) for m ∈ N] for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z  = f(s)
        # sample a random node
        n  = sample(rng, N, Weights(W))
        vₙ = V[n.v]
        # sample another node from the same route, weighted by relatedness
        m  = sample(rng, N, Weights(R[n.i]))
        vₘ = V[m.v]
        # perform operations
        if isdepot(m) continue end
        if isequal(n,m) continue end
        tₙ = N[n.i]
        hₙ = N[n.h]
        tₘ = N[m.t]
        hₘ = N[m.h]
        while !isequal(m, hₒ)
            removenode!(m, tₘ, hₘ, vₘ, s)
            insertnode!(m, tₙ, hₙ, vₙ, s)
            tₙ = m 
            hₙ = N[tₙ.h]
            m  = tₘ
            tₘ = N[m.t]
            hₘ = N[m.h]
        end
        # evaluate the change in objective function value
        z′ = f(s)
        c  = z′ - z
        if c < 0 continue end
        # reperform operations
        tₙ = N[n.i]
        hₙ = N[n.h]
        tₘ = N[m.t]
        hₘ = N[m.h]
        while !isequal(m, hₒ)
            removenode!(m, tₘ, hₘ, vₘ, s)
            insertnode!(m, tₙ, hₙ, vₙ, s)
            tₙ = m 
            hₙ = N[tₙ.h]
            m  = tₘ
            tₘ = N[m.t]
            hₘ = N[m.h]
        end
    end
    # return solution
    return s
end

function interopt!(rng::AbstractRNG, k::Int, s::Solution)
    # initialize
    G = s.G
    N = G.N
    V = G.V
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    R = [[(isequal(n, m) || isequal(n.v, m.v)) ? 0. : relatedness(n, m) for m ∈ N] for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z  = f(s)
        # sample a random node
        n  = sample(rng, N, Weights(W))
        vₙ = V[n.v]
        # sample another node from the same route, weighted by relatedness
        m  = sample(rng, N, Weights(R[n.i]))
        vₘ = V[m.v]
        # perform operationss
        if isdepot(m) continue end
        if isequal(n,m) continue end
        tₙ = N[n.i]
        hₙ = N[n.h]
        o  = N[m.i]
        tₘ = N[m.t]
        hₘ = N[m.h]
        while !isdepot(o)
            removenode!(o, tₘ, hₘ, vₘ, s)
            insertnode!(o, tₙ, hₙ, vₙ, s)
            tₙ = o 
            hₙ = N[tₙ.h]
            o  = hₘ
            tₘ = N[o.t]
            hₘ = N[o.h]
        end
        m  = hₙ
        tₘ = N[vₘ.e]
        hₘ = N[1]
        o  = hₙ
        tₙ = N[o.t]
        hₙ = N[o.h]
        while !isdepot(o)
            removenode!(o, tₙ, hₙ, vₙ, s)
            insertnode!(o, tₘ, hₘ, vₘ, s)
            tₘ = o
            hₘ = N[tₘ.h]
            o  = hₙ
            tₙ = N[o.t]
            hₙ = N[o.h]
        end
        # evaluate the change in objective function value
        z′ = f(s)
        c  = z′ - z
        if c < 0 continue end
        # reperform operations
        tₙ = N[n.i]
        hₙ = N[n.h]
        o  = N[m.i]
        tₘ = N[m.t]
        hₘ = N[m.h]
        while !isdepot(o)
            removenode!(o, tₘ, hₘ, vₘ, s)
            insertnode!(o, tₙ, hₙ, vₙ, s)
            tₙ = o 
            hₙ = N[tₙ.h]
            o  = hₘ
            tₘ = N[o.t]
            hₘ = N[o.h]
        end
        tₘ = N[vₘ.e]
        hₘ = N[1]
        o  = hₙ
        tₙ = N[o.t]
        hₙ = N[o.h]
        while !isdepot(o)
            removenode!(o, tₙ, hₙ, vₙ, s)
            insertnode!(o, tₘ, hₘ, vₘ, s)
            tₘ = o
            hₘ = N[tₘ.h]
            o  = hₙ
            tₙ = N[o.t]
            hₙ = N[o.h]
        end
    end
    # return solution
    return s
end