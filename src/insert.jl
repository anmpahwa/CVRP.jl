function best!(rng::AbstractRNG, s::Solution; mode::Symbol)
    # initialize
    φ = isequal(mode, :ptb)
    N = s.N
    V = s.V
    L = [n for n ∈ N if iscustomer(n) && isopen(n)] # set of open nodes
    I = eachindex(L)
    J = eachindex(V)
    W = ones(Int, I)                                # W[i]: selection weight for node L[i]
    C = fill(Inf, I)                                # C[i]: best insertion cost of node L[i] 
    P = fill((0,0,0), I)                            # P[i]: best insertion position of node L[i]
    # loop: until all nodes are inserted
    for _ ∈ I
        z = f(s)
        # sample a node and iterate through all possible insertion positions in each route
        i = sample(rng, I, Weights(W))
        n = L[i]
        for (j,v) ∈ enumerate(V)
            t = N[1]
            h = N[v.s]
            for _ ∈ 0:v.n
                # insert the node
                insertnode!(n, t, h, v, s)
                # evaluate the insertion cost
                z′ = f(s) * (1 + φ * rand(rng, Uniform(-0.2, 0.2)))
                c  = z′ - z
                # reset the best insertion position for the node
                if c < C[i] C[i], P[i] = c, (t.i, h.i, v.i) end
                # remove the node
                removenode!(n, t, h, v, s)
                t = h
                h = N[t.h]
            end
        end
        # insert the node at its best position
        p = P[i]
        t = N[p[1]]
        h = N[p[2]]
        v = V[p[3]]
        insertnode!(n, t, h, v, s)
        # update node vectors
        W[i] = 0
        C[i] = Inf
        P[i] = (0,0,0)
    end
    # return solution
    return s
end

function greedy!(rng::AbstractRNG, s::Solution; mode::Symbol)
    # initialize
    φ = isequal(mode, :ptb)
    N = s.N
    V = s.V 
    L = [n for n ∈ N if iscustomer(n) && isopen(n)] # set of open nodes
    I = eachindex(L)
    J = eachindex(V)
    W = ones(Int, I)                                # W[i]: selection weight for node L[i]
    C = fill(Inf, I)                                # C[i]: best insertion cost of node L[i] 
    P = fill((0,0,0), I)                            # P[i]: best insertion position of node L[i]
    # loop: until all nodes are inserted
    for _ ∈ I
        z = f(s)
        for (i,n) ∈ enumerate(L)
            if isclose(n) continue end
            for (j,v) ∈ enumerate(V)
                t = N[1]
                h = N[v.s]
                for _ ∈ 0:v.n
                    # insert the node
                    insertnode!(n, t, h, v, s)
                    # evaluate the insertion cost
                    z′ = f(s) * (1 + φ * rand(rng, Uniform(-0.2, 0.2)))
                    c  = z′ - z
                    # reset the best insertion position for the node
                    if c < C[i] C[i], P[i] = c, (t.i, h.i, v.i) end
                    # remove the node
                    removenode!(n, t, h, v, s)
                    t = h
                    h = N[t.h]
                end
            end
        end
        # insert the node with the minimum insertion cost
        i = argmin(C)
        n = L[i]
        p = P[i]
        t = N[p[1]]
        h = N[p[2]]
        v = V[p[3]]
        insertnode!(n, t, h, v, s)
        # update node vectors
        W[i] = 0
        C[i] = Inf
        P[i] = (0,0,0)
    end
    # return solution
    return s
end

function regretk!(rng::AbstractRNG, s::Solution; k::Int)
    # initialize
    φ = isequal(mode, :ptb)
    N = s.N
    V = s.V 
    L = [n for n ∈ N if iscustomer(n) && isopen(n)] # set of open nodes
    I = eachindex(L)
    J = eachindex(V)
    W = ones(Int, I)                                # W[i]      : selection weight for node L[i]
    C = fill(Inf, I)                                # C[i]      : best insertion cost of node L[i] 
    Cₖ= fill(Inf, (I,k))                            # Cₖ[i,k]   : k-th least insertion cost of node L[i]
    R = fill(Inf, I)                                # R[i]      : regret cost of node L[i]
    P = fill((0,0,0), I)                            # P[i]      : best insertion position of node L[i]
    # loop: until all nodes are inserted
    for _ ∈ I
        z = f(s)
        for (i,n) ∈ enumerate(L)
            if isclose(n) continue end
            for (j,v) ∈ enumerate(V)
                t = N[1]
                h = N[v.s]
                for _ ∈ 0:v.n
                    # insert the node
                    insertnode!(n, t, h, v, s)
                    # evaluate the insertion cost
                    z′ = f(s) * (1 + φ * rand(rng, Uniform(-0.2, 0.2)))
                    c  = z′ - z
                    # reset the best insertion position for the node
                    if c < C[i] C[i], P[i] = c, (t.i, h.i, v.i) end
                    # revise k least insertion cost
                    for (k,cₖ) ∈ enumerate(Cₖ[i])
                        if c < cₖ
                            Cₖ[i,k] = c
                            c = cₖ
                        end
                    end
                    # remove the node
                    removenode!(n, t, h, v, s)
                    t = h
                    h = N[t.h]
                end
            end
            # compute the regret cost
            cₒ = Cₖ[i,1]
            for (k,cₖ) ∈ enumerate(Cₖ[i]) R[i] += cₖ - cₒ end
        end
        # insert the node with the maximum regret cost
        I̲ = findall(isequal.(R, maximum(R)))
        i = argmin(C[I̲])
        n = L[i]
        p = P[i]
        t = N[p[1]]
        h = N[p[2]]
        v = V[p[3]]
        insertnode!(n, t, h, v, s)
        # update node vectors
        W[i] = 0
        C[i] = Inf
        P[i] = (0,0,0)
    end
    # return solution
    return s        
end