function intramove!(rng::AbstractRNG, k::Int, s::Solution)
    # initialize
    N = s.N
    V = s.V
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample a random node
        n = sample(rng, N, Weights(W))
        # remove the node from its current position
        tₙ= N[n.t]
        hₙ= N[n.h]
        vₙ= V[n.v]
        removenode!(n, tₙ, hₙ, vₙ, s)
        # iterate through all positions in the route
        v = vₙ
        c = 0.
        p = (tₙ.i, hₙ.i, vₙ.i)
        t = N[1]
        h = N[v.s]
        for _ ∈ 0:v.n
            # insert the node
            insertnode!(n, t, g, v, s)
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

function intermove!(rng::AbstractRNG, k::Int, s::Solution)
    # initialize
    N = s.N
    V = s.V
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    R = [[isequal(u, v) ? 0. : relatedness(u, v) for u ∈ V] for v ∈ V]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample a random node
        n = sample(rng, N, Weights(W))
        # remove the node from its current position
        tₙ= N[n.t]
        hₙ= N[n.h]
        vₙ= V[n.v]
        removenode!(n, tₙ, hₙ, vₙ, s)
        # sample a random vehicle and iterate through all positions in it
        v = sample(rng, V, Weights(R[vₙ.i]))
        c = 0.
        p = (tₙ.i, hₙ.i, vₙ.i)
        t = N[1]
        h = N[v.s]
        for _ ∈ 0:v.n
            # insert the node
            insertnode!(n, t, g, v, s)
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

function intraswap!(rng::AbstractRNG, k::Int, s::Solution)
    # tₘ → m → hₘ and tₙ → n → hₙ
    # initialize
    N = s.N
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    R = [[(isequal(n, m) || !isequal(n.v, m.v)) ? 0. : relatedness(n, m) for m ∈ N] for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample a random node
        n = sample(rng, N, Weights(W))
        # sample another node from the same route, weighted by relatedness
        m = sample(rng, N, Weights(R[n.i]))
        # swap the two nodes
        tₙ= N[n.t]
        hₙ= N[n.h]
        tₘ= N[m.t]
        hₘ= N[m.h]
        # case 1: tₘ → m (tₙ) → n (hₘ) → hₙ
        if isequal(n, hₘ) # || isequal(m, tₙ)
            removenode!(n, tₙ, hₙ, v, s)
            insertnode!(n, tₘ, tₙ, v, s)
        # case 2: tₙ → n (tₘ) → m (hₙ) → hₘ
        elseif isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(n, tₙ, hₙ, v, s)
            insertnode!(n, hₙ, hₘ, v, s)
        # case 3: all other cases
        else
            removenode!(n, tₙ, hₙ, v, s)
            removenode!(m, tₘ, hₘ, v, s)
            insertnode!(n, tₘ, hₘ, v, s)
            insertnode!(m, tₙ, hₙ, v, s)
        end
        # evaluate the change in objective function value
        z′ = f(s)
        c  = z′ - z
        if c < 0 continue end
        # reswap the two nodes
        tₙ= N[n.t]
        hₙ= N[n.h]
        tₘ= N[m.t]
        hₘ= N[m.h]
        # case 1: tₘ → m (tₙ) → n (hₘ) → hₙ
        if isequal(n, hₘ) # || isequal(m, tₙ)
            removenode!(n, tₙ, hₙ, v, s)
            insertnode!(n, tₘ, tₙ, v, s)
        # case 2: tₙ → n (tₘ) → m (hₙ) → hₘ
        elseif isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(n, tₙ, hₙ, v, s)
            insertnode!(n, hₙ, hₘ, v, s)
        # case 3: all other cases
        else
            removenode!(n, tₙ, hₙ, v, s)
            removenode!(m, tₘ, hₘ, v, s)
            insertnode!(n, tₘ, hₘ, v, s)
            insertnode!(m, tₙ, hₙ, v, s)
        end
    end
    # return solution
    return s
end

function interswap!(rng::AbstractRNG, k::Int, s::Solution)
    # tₘ → m → hₘ and tₙ → n → hₙ
    # initialize
    N = s.N
    W = [isone(i) ? 0 : 1 for i ∈ eachindex(s.N)]
    R = [[(isequal(n, m) || isequal(n.v, m.v)) ? 0 : relatedness(n, m) for m ∈ N] for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample a random node
        n = sample(rng, N, Weights(W))
        # sample another node from the same route, weighted by relatedness
        m = sample(rng, N, Weights(R[n.i]))
        # swap the two nodes
        if isequal(n,m) continue end
        tₙ= N[n.t]
        hₙ= N[n.h]
        tₘ= N[m.t]
        hₘ= N[m.h]
        # case 1: tₘ → m (tₙ) → n (hₘ) → hₙ
        if isequal(n, hₘ) # || isequal(m, tₙ)
            removenode!(n, tₙ, hₙ, v, s)
            insertnode!(n, tₘ, tₙ, v, s)
        # case 2: tₙ → n (tₘ) → m (hₙ) → hₘ
        elseif isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(n, tₙ, hₙ, v, s)
            insertnode!(n, hₙ, hₘ, v, s)
        # case 3: all other cases
        else
            removenode!(n, tₙ, hₙ, v, s)
            removenode!(m, tₘ, hₘ, v, s)
            insertnode!(n, tₘ, hₘ, v, s)
            insertnode!(m, tₙ, hₙ, v, s)
        end
        # evaluate the change in objective function value
        z′ = f(s)
        c  = z′ - z
        if c < 0 continue end
        # reswap the two nodes
        tₙ= N[n.t]
        hₙ= N[n.h]
        tₘ= N[m.t]
        hₘ= N[m.h]
        # case 1: tₘ → m (tₙ) → n (hₘ) → hₙ
        if isequal(n, hₘ) # || isequal(m, tₙ)
            removenode!(n, tₙ, hₙ, v, s)
            insertnode!(n, tₘ, tₙ, v, s)
        # case 2: tₙ → n (tₘ) → m (hₙ) → hₘ
        elseif isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(n, tₙ, hₙ, v, s)
            insertnode!(n, hₙ, hₘ, v, s)
        # case 3: all other cases
        else
            removenode!(n, tₙ, hₙ, v, s)
            removenode!(m, tₘ, hₘ, v, s)
            insertnode!(n, tₘ, hₘ, v, s)
            insertnode!(m, tₙ, hₙ, v, s)
        end
    end
    # return solution
    return s
end

function intraopt!(rng::AbstractRNG, k::Int, s::Solution)
    # initialize
    N = s.N
    V = s.V
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    R = [[(isequal(n, m) || !isequal(n.v, m.v)) ? 0. : relatedness(n, m) for m ∈ N] for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample a random node
        n = sample(rng, N, Weights(W))
        # sample another node from the same route, weighted by relatedness
        m = sample(rng, N, Weights(R[n.i]))
        # perform operations
        if isequal(n,m) continue end
        v = V[n.v]
        tₒ= n
        hₙ= N[n.h]
        tₘ= N[m.t]
        hₘ= N[m.h]
        while !isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(m, tₘ, hₘ, v, s)
            insertnode!(m, tₒ, hₙ, v, s)
            tₒ= m 
            m = tₘ
            hₙ= N[n.h]
            tₘ= N[m.t]
            hₘ= N[m.h]
        end
        # evaluate the change in objective function value
        z′ = f(s)
        c  = z′ - z
        if c < 0 continue end
        # reperform operations
        v = V[n.v]
        tₒ= n
        hₙ= N[n.h]
        tₘ= N[m.t]
        hₘ= N[m.h]
        while !isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(m, tₘ, hₘ, v, s)
            insertnode!(m, tₒ, hₙ, v, s)
            tₒ= m 
            m = tₘ
            hₙ= N[n.h]
            tₘ= N[m.t]
            hₘ= N[m.h]
        end
    end
end

function interopt!(rng::AbstractRNG, k::Int, s::Solution)
    # initialize
    N = s.N
    V = s.V
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    R = [[(isequal(n, m) || isequal(n.v, m.v)) ? 0. : relatedness(n, m) for m ∈ N] for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample a random node
        n = sample(rng, N, Weights(W))
        # sample another node from the same route, weighted by relatedness
        m = sample(rng, N, Weights(R[n.i]))
        # perform operations
        if isequal(n,m) continue end
        v = V[n.v]
        tₒ= n
        hₙ= N[n.h]
        tₘ= N[m.t]
        hₘ= N[m.h]
        while !isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(m, tₘ, hₘ, v, s)
            insertnode!(m, tₒ, hₙ, v, s)
            tₒ= m 
            m = tₘ
            hₙ= N[n.h]
            tₘ= N[m.t]
            hₘ= N[m.h]
        end
        # evaluate the change in objective function value
        z′ = f(s)
        c  = z′ - z
        if c < 0 continue end
        # reperform operations
        v = V[n.v]
        tₒ= n
        hₙ= N[n.h]
        tₘ= N[m.t]
        hₘ= N[m.h]
        while !isequal(n, tₘ) # || isequal(m, hₙ)
            removenode!(m, tₘ, hₘ, v, s)
            insertnode!(m, tₒ, hₙ, v, s)
            tₒ= m 
            m = tₘ
            hₙ= N[n.h]
            tₘ= N[m.t]
            hₘ= N[m.h]
        end
    end
end