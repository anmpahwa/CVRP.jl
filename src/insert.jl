"""
    insert!([rng::AbstractRNG], s::Solution, method::Function)

Returns solution `s` after inserting open customer nodes to the solution using the given `method`.

Available methods include,
- Best Precise Insertion   : `bestprecise!`
- Best Perturb Insertion   : `bestperturb!`
- Greedy Precise Insertion : `greedyprecise!`
- Greedy Perturb Insertion : `greedyperturb!`
- Regret₂ Precise Insertion: `regret2precise!`
- Regret₂ Perturb Insertion: `regret2perturb!`
- Regret₃ Precise Insertion: `regret3precise!`
- Regret₃ Perturb Insertion: `regret3perturb!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
insert!(rng::AbstractRNG, s::Solution, method::Function)::Solution = method(rng, s)
insert!(s::Solution, method::Function) = insert!(Random.GLOBAL_RNG, s, method)

"""
    best!(rng::AbstractRNG, s::Solution; mode::Symbol)

Returns solution `s` after inserting all open customer nodes into the solution `s`.
In each iteration, the algorithm selects a random open node and places it at its best 
position in the solution according to the specified mode (`:precise` or `:perturb`).
"""
function best!(rng::AbstractRNG, s::Solution; mode::Symbol)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # initialize
    φ = isequal(mode, :perturb)
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
        p = P[i,j]
        t = N[p[1]]
        h = N[p[2]]
        v = V[j]
        insertnode!(n, t, h, v, s)
        # update node vectors
        W[i] = 0
        C[i,:] .= Inf
        P[i,:] .= ((0,0),)
    end
    # return solution
    return s
end
"""
    bestprecise!(rng::AbstractRNG, s::Solution)

Returns solution `s` after inserting all open customer nodes into the solution `s`.
In each iteration, the algorithm selects a random open node and places it at its best 
position in the solution, evaluating insertion costs with precision.
"""
bestprecise!(rng::AbstractRNG, s::Solution) = best!(rng, s; mode=:precise)
"""
    bestperturb(rng::AbstractRNG, s::Solution)

Returns solution `s` after inserting all open customer nodes into the solution `s`.
In each iteration, the algorithm selects a random open node and places it at its best 
position in the solution, evaluating insertion costs with perturbation.
"""
bestperturb!(rng::AbstractRNG, s::Solution) = best!(rng, s; mode=:perturb)

"""
    greedy!(rng::AbstractRNG, s::Solution; mode::Symbol)

Returns solution `s` after inserting all open customer nodes into the solution `s`.
In each iteration, the algorithm selects the node with the minimum insertion cost
and places it at its best position in the solution according to the specified mode 
(`:precise` or `:perturb`).
"""
function greedy!(rng::AbstractRNG, s::Solution; mode::Symbol)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # initialize
    φ = isequal(mode, :perturb)
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
        P[i,:] .= ((0,0),)
        C[:,j] .= Inf
        P[:,j] .= ((0,0),)
    end
    # return solution
    return s
end
"""
    greedyprecise!(rng::AbstractRNG, s::Solution)
Returns solution `s` after inserting all open customer nodes into the solution `s`.
In each iteration, the algorithm selects the node with the minimum insertion cost
and places it at its best position in the solution, evaluating insertion costs with 
precision.
"""
greedyprecise!(rng::AbstractRNG, s::Solution) = greedy!(rng, s; mode=:precise)
"""
    greedyperturb!(rng::AbstractRNG, s::Solution)
Returns solution `s` after inserting all open customer nodes into the solution `s`.
In each iteration, the algorithm selects the node with the minimum insertion cost
and places it at its best position in the solution, evaluating insertion costs with 
perturbation.
"""
greedyperturb!(rng::AbstractRNG, s::Solution) = greedy!(rng, s; mode=:perturb)

"""
    regretk!(rng::AbstractRNG, s::Solution; k::Int, mode::Symbol)

Returns solution `s` after inserting all open customer nodes into the solution `s`.
In each iteration, the algorithm selects the node with the maximum regret cost
and places it at its best position in the solution according to the specified mode
(`:precise` or `:perturb`). The regret cost is calculated based on the difference 
between the best insertion cost and the k-th best insertion cost.
"""
function regretk!(rng::AbstractRNG, s::Solution; k::Int, mode::Symbol)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # initialize
    φ = isequal(mode, :perturb)
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
        P[i,:] .= ((0,0),)
        C[:,j] .= Inf
        P[:,j] .= ((0,0),)
    end
    # return solution
    return s        
end
"""
    regret2precise!(rng::AbstractRNG, s::Solution)

Returns solution `s` after inserting all open customer nodes into the solution `s`.
In each iteration, the algorithm selects the node with the maximum regret cost
and places it at its best position, evaluating insertion cost with precision. 
The regret cost is calculated based on the difference between the best insertion 
cost and the 2nd best insertion cost.
"""
regret2precise!(rng::AbstractRNG, s::Solution) = regretk!(rng, s; k=2, mode=:precise)
"""
    regret2perturb!(rng::AbstractRNG, s::Solution)

Returns solution `s` after inserting all open customer nodes into the solution `s`.
In each iteration, the algorithm selects the node with the maximum regret cost
and places it at its best position, evaluating insertion cost with perturbation. 
The regret cost is calculated based on the difference between the best insertion 
cost and the 2nd best insertion cost.
"""
regret2perturb!(rng::AbstractRNG, s::Solution) = regretk!(rng, s; k=2, mode=:perturb)
"""
    regret3precise!(rng::AbstractRNG, s::Solution)

Returns solution `s` after inserting all open customer nodes into the solution `s`.
In each iteration, the algorithm selects the node with the maximum regret cost
and places it at its best position, evaluating insertion cost with precision. 
The regret cost is calculated based on the difference between the best insertion 
cost and the 3rd best insertion cost.
"""
regret3precise!(rng::AbstractRNG, s::Solution) = regretk!(rng, s; k=3, mode=:precise)
"""
    regret2perturb!(rng::AbstractRNG, s::Solution)

Returns solution `s` after inserting all open customer nodes into the solution `s`.
In each iteration, the algorithm selects the node with the maximum regret cost
and places it at its best position, evaluating insertion cost with perturbation. 
The regret cost is calculated based on the difference between the best insertion 
cost and the 3rd best insertion cost.
"""
regret3perturb!(rng::AbstractRNG, s::Solution) = regretk!(rng, s; k=3, mode=:perturb)