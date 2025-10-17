"""
    localsearch!([rng::AbstractRNG], k::Int, s::Solution, method::Function)

Returns solution `s` after performing local seach on the solution using given `method` for `k` iterations.

Available methods include,
- intra-move    : `intramove!`
- inter-move    : `intermove!`
- intra-swap    : `intraswap!`
- inter-swap    : `interswap!`
- intra-opt     : `intraopt!`
- inter-opt     : `interopt!`

Optionally specify a random number generator `rng` as the first argument (defaults to `Random.GLOBAL_RNG`).
"""
localsearch!(rng::AbstractRNG, k::Int, s::Solution, method::Function)::Solution = method(rng, k, s)
localsearch!(k::Int, s::Solution, method::Function) = localsearch!(Random.GLOBAL_RNG, k, s, method)

"""
    move!(rng::AbstractRNG, k::Int, s::Solution; scope::Symbol)

Returns solution `s` after moving a randomly selected node to its best position 
in the specified solution `scope`, if the move results in a reduction in objective 
function value, repeating for `k` iterations.
"""
function move!(rng::AbstractRNG, k::Int, s::Solution; scope::Symbol)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # initialize
    Wₙ = [isdepot(n) ? 0 : 1 for n ∈ N]
    Wᵥ = [[isequal(n.v, v.i) ? isequal(scope, :intra) : isequal(scope, :inter) for v ∈ V] for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample a random node
        n = sample(rng, N, Weights(Wₙ))
        # remove the node from its current position
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        # sample a random vehicle and iterate through all positions in it
        c = 0.
        p = (t.i, h.i, v.i)
        v = sample(rng, V, Weights(Wᵥ[n.i]))
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
    intramove!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after moving a randomly selected node to its best position 
in the same route, if the move results in a reduction in objective function value, 
repeating for `k` iterations.    
"""
intramove!(rng::AbstractRNG, k::Int, s::Solution) = move!(rng, k, s; scope=:intra)
"""
    intermove!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after moving a randomly selected node to its best position 
in a randomly selected route, if the move results in a reduction in objective 
function value, repeating for `k` iterations.    
"""
intermove!(rng::AbstractRNG, k::Int, s::Solution) = move!(rng, k, s; scope=:inter)

"""
    intraswap!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after swapping two randomly selected nodes in the specified 
solution `scope`, if the swap results in a reduction in objective function 
value, repeating for `k` iterations.
"""
function swap!(rng::AbstractRNG, k::Int, s::Solution; scope::Symbol)
    # tₘ → m → hₘ and tₙ → n → hₙ
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # initialize
    Wₙ = [isdepot(n) ? 0 : 1 for n ∈ N]
    Wₘ = [[isequal(n.v, m.v) ? isequal(scope, :intra) : isequal(scope, :inter) for m ∈ N] for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample nodes
        n = sample(rng, N, Weights(Wₙ))
        m = sample(rng, N, Weights(Wₘ[n.i] .* Wₙ))
        # swap the two nodes
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
    intraswap!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after swapping two randomly selected nodes from 
the same route, if the swap results in a reduction in objective function 
value, repeating for `k` iterations.
"""
intraswap!(rng::AbstractRNG, k::Int, s::Solution) = swap!(rng, k, s; scope=:intra)
"""
    interswap!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after swapping two randomly selected nodes from 
two different routes, if the swap results in a reduction in objective 
function value, repeating for `k` iterations.
"""
interswap!(rng::AbstractRNG, k::Int, s::Solution) = swap!(rng, k, s; scope=:inter)

"""
    intraopt!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after iteratively taking 2 arcs from the same route 
and reconfiguring them (total possible reconfigurations 2₂-1 = 3) if the 
reconfiguration results in a reduction in objective function value, repeating 
for `k` iterations.
"""
function intraopt!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # initialize
    Wₙ = [isdepot(n) ? 0 : 1 for n ∈ N]
    Wₘ = [[isequal(n.v, m.v) ? (isequal(n, m) ? 0 : 1) : 0 for m ∈ N] for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample nodes
        n = sample(rng, N, Weights(Wₙ))
        m = sample(rng, N, Weights(Wₘ[n.i]))
        # order correction
        p = N[n.t]
        while !isdepot(p)
            if isequal(p, m)
                m = n
                n = p
                break
            end
            p = N[p.t]
        end
        # skip invalid configurations
        if isdepot(n) continue end
        if isdepot(m) continue end
        if isequal(n, m) continue end
        if !isequal(n.v, m.v) continue end
        # perform operations
        # fetch vehicle
        vₒ = V[n.v]
        # set pivot node
        p  = N[m.i]
        # fetch tail and head node of new position
        tₙ = N[n.i]
        hₙ = N[n.h]
        # fetch tail and head node of original position
        tₒ = N[p.t]
        hₒ = N[p.h]
        # move nodes
        b = N[n.h]
        while !isequal(p.i, b.i)
            removenode!(p, tₒ, hₒ, vₒ, s)
            insertnode!(p, tₙ, hₙ, vₒ, s)
            tₙ = N[p.i]
            hₙ = N[p.h]
            p  = N[tₒ.i]
            tₒ = N[p.t]
            hₒ = N[p.h]
        end
        m = p
        # evaluate the change in objective function value
        z′ = f(s)
        c  = z′ - z
        if c < 0 continue end
        # reperform operations
        # set pivot node
        p  = N[m.i]
        # fetch tail and head node of new position
        tₙ = N[n.i]
        hₙ = N[n.h]
        # fetch tail and head node of original position
        tₒ = N[p.t]
        hₒ = N[p.h]
        # move nodes
        b = N[n.h]
        while !isequal(p.i, b.i)
            removenode!(p, tₒ, hₒ, vₒ, s)
            insertnode!(p, tₙ, hₙ, vₒ, s)
            tₙ = N[p.i]
            hₙ = N[p.h]
            p  = N[tₒ.i]
            tₒ = N[p.t]
            hₒ = N[p.h]
        end
        m = p
    end
    # return solution
    return s
end

"""
    interopt!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after iteratively taking 2 arcs from different routes 
and reconfiguring them (total possible reconfigurations 2₂-1 = 3) if the 
reconfiguration results in a reduction in objective function value, repeating 
for `k` iterations.
"""
function interopt!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # initialize
    Wₙ = [isdepot(n) ? 0 : 1 for n ∈ N]
    Wₘ = [[isequal(n.v, m.v) ? 0 : 1 for m ∈ N] for n ∈ N]
    # iterate
    for _ ∈ 1:k
        # compute cost of the current solution
        z = f(s)
        # sample nodes
        n = sample(rng, N, Weights(Wₙ))
        m = sample(rng, N, Weights(Wₘ[n.i]))
        # skip invalid configurations
        if isdepot(n) continue end
        if isdepot(m) continue end
        if isequal(n.v, m.v) continue end
        # perform operations
        # fetch vehicles
        vₙ = V[n.v]
        vₘ = V[m.v]
        # fetch tail and head nodes in vehicle n
        tₙ = N[n.i]
        hₙ = N[n.h]
        # set pivot node
        p  = N[m.i]
        # fetch tail and head nodes in vehicle m
        tₘ = N[p.t]
        hₘ = N[p.h]
        # move nodes from vehicle m to n
        while !isdepot(p)
            removenode!(p, tₘ, hₘ, vₘ, s)
            insertnode!(p, tₙ, hₙ, vₙ, s)
            # update nodes
            tₙ = p
            hₙ = N[tₙ.h]
            p  = hₘ
            tₘ = N[p.t]
            hₘ = N[p.h]
        end
        m = hₙ
        # fetch tail and head nodes in vehicle m
        tₘ = N[vₘ.e]
        hₘ = N[1]
        # set pivot node
        p  = hₙ
        # fetch tail and head nodes in vehicle n
        tₙ = N[p.t]
        hₙ = N[p.h]
        # move nodes from vehicle n to m
        while !isdepot(p)
            removenode!(p, tₙ, hₙ, vₙ, s)
            insertnode!(p, tₘ, hₘ, vₘ, s)
            # update nodes
            tₘ = p
            hₘ = N[tₘ.h]
            p  = hₙ
            tₙ = N[p.t]
            hₙ = N[p.h]
        end
        # evaluate the change in objective function value
        z′ = f(s)
        c  = z′ - z
        if c < 0 continue end
        # reperform operations
        # fetch tail and head nodes in vehicle n
        tₙ = N[n.i]
        hₙ = N[n.h]
        # set pivot node
        p  = N[m.i]
        # fetch tail and head nodes in vehicle m
        tₘ = N[p.t]
        hₘ = N[p.h]
        # move nodes from vehicle m to n
        while !isdepot(p)
            removenode!(p, tₘ, hₘ, vₘ, s)
            insertnode!(p, tₙ, hₙ, vₙ, s)
            # update nodes
            tₙ = p
            hₙ = N[tₙ.h]
            p  = hₘ
            tₘ = N[p.t]
            hₘ = N[p.h]
        end
        m = hₙ
        # fetch tail and head nodes in vehicle m
        tₘ = N[vₘ.e]
        hₘ = N[1]
        # set pivot node
        p  = hₙ
        # fetch tail and head nodes in vehicle n
        tₙ = N[p.t]
        hₙ = N[p.h]
        # move nodes from vehicle n to m
        while !isdepot(p)
            removenode!(p, tₙ, hₙ, vₙ, s)
            insertnode!(p, tₘ, hₘ, vₘ, s)
            # update nodes
            tₘ = p
            hₘ = N[tₘ.h]
            p  = hₙ
            tₙ = N[p.t]
            hₙ = N[p.h]
        end
    end
    # return solution
    return s
end