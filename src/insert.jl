function best!(rng::AbstractRNG, s::Solution; mode::Symbol)
    # initialize
    φ = isequal(mode, :ptb)
    G = s.G
    N = G.N
    V = G.V
    L = [n for n ∈ N if iscustomer(n) && isopen(n)] # set of open nodes
    I = eachindex(L)
    J = eachindex(V)
    W = ones(Int, I)                                # W[i]: selection weight for node L[i]
    C = fill(Inf, (I,J))                            # C[i,j]: best insertion cost of node L[i] in vehicle route V[j]
    P = fill((0,0), (I,J))                          # P[i,j]: best insertion position on node L[i] in vehicle route V[j]
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
                if c < C[i,j] C[i,j], P[i,j] = c, (t.i, h.i) end
                # remove the node
                removenode!(n, t, h, v, s)
                t = h
                h = N[t.h]
            end
        end
        # insert the node at its best position
        j = argmin(C[i,:])
        p = P[i]
        t = N[p[1]]
        h = N[p[2]]
        v = V[j]
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
    G = s.G
    N = G.N
    V = G.V 
    L = [n for n ∈ N if iscustomer(n) && isopen(n)] # set of open nodes
    I = eachindex(L)
    J = eachindex(V)
    ϕ = ones(Int, J)                                # ϕ[j]  : binary weight for vehicle route V[j]
    C = fill(Inf, (I,J))                            # C[i,j]: best insertion cost of node L[i] in vehicle route V[j]
    P = fill((0,0), (I,J))                          # P[i,j]: best insertion position on node L[i] in vehicle route V[j]
    # loop: until all nodes are inserted
    for _ ∈ I
        z = f(s)
        for (i,n) ∈ enumerate(L)
            if isclose(n) continue end
            for (j,v) ∈ enumerate(V)
                if iszero(ϕ[j]) continue end
                t = N[1]
                h = N[v.s]
                for _ ∈ 0:v.n
                    # insert the node
                    insertnode!(n, t, h, v, s)
                    # evaluate the insertion cost
                    z′ = f(s) * (1 + φ * rand(rng, Uniform(-0.2, 0.2)))
                    c  = z′ - z
                    # reset the best insertion position for the node
                    if c < C[i,j] C[i,j], P[i,j] = c, (t.i, h.i) end
                    # remove the node
                    removenode!(n, t, h, v, s)
                    t = h
                    h = N[t.h]
                end
            end
        end
        # insert the node with the minimum insertion cost
        i,j = Tuple(argmin(C))
        p = P[i,j]
        n = L[i]
        t = N[p[1]]
        h = N[p[2]]
        v = V[j]
        insertnode!(n, t, h, v, s)
        # update node vectors
        ϕ .= 0
        ϕ[j] = 1
        C[i,:] .= Inf
        P[i,:] .= ((0, 0), )
        C[:,j] .= Inf
        P[:,j] .= ((0, 0), )
    end
    # return solution
    return s
end

function regretk!(rng::AbstractRNG, s::Solution; k::Int)
    # initialize
    φ = isequal(mode, :ptb)
    G = s.G
    N = G.N
    V = G.V 
    L = [n for n ∈ N if iscustomer(n) && isopen(n)] # set of open nodes
    I = eachindex(L)
    J = eachindex(V)
    ϕ = ones(Int, J)                                # ϕ[j]      : binary weight for vehicle route V[j]
    C = fill(Inf, (I,J))                            # C[i,j]    : best insertion cost of node L[i] in vehicle route V[j]
    P = fill((0,0), (I,J))                          # P[i,j]    : best insertion position on node L[i] in vehicle route V[j]
    Cₖ= fill(Inf, (I,k))                            # Cₖ[i,k]   : k-th least insertion cost of node L[i]
    R = fill(Inf, I)                                # R[i]      : regret cost of node L[i]
    # loop: until all nodes are inserted
    for _ ∈ I
        z = f(s)
        for (i,n) ∈ enumerate(L)
            if isclose(n) continue end
            for (j,v) ∈ enumerate(V)
                if iszero(ϕ[j]) continue end
                t = N[1]
                h = N[v.s]
                for _ ∈ 0:v.n
                    # insert the node
                    insertnode!(n, t, h, v, s)
                    # evaluate the insertion cost
                    z′ = f(s) * (1 + φ * rand(rng, Uniform(-0.2, 0.2)))
                    c  = z′ - z
                    # reset the best insertion position for the node
                    if c < C[i,j] C[i,j], P[i,j] = c, (t.i, h.i) end
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
        i,j = Tuple(argmin(C[I̲,:]))
        p = P[i,j]
        n = L[i]
        t = N[p[1]]
        h = N[p[2]]
        v = V[j]
        insertnode!(n, t, h, v, s)
        # update node vectors
        ϕ .= 0
        ϕ[j] = 1
        C[i,:] .= Inf
        P[i,:] .= ((0, 0), )
        C[:,j] .= Inf
        P[:,j] .= ((0, 0), )
    end
    # return solution
    return s        
end