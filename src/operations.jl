function removenode!(s::Solution, n::Node, p::Node, q::Node, v::Vehicle)
    # fetch features
    x = v.x * v.l
    y = v.y * v.l
    aₜᵢ = s.A[p.i, n.i]
    aᵢₕ = s.A[n.i, q.i]
    aₜₕ = s.A[p.i, q.i]
    δ = (aₜᵢ.l + aᵢₕ.l) - aₜₕ.l

    # update node
    n.t = 0
    n.h = 0
    n.v = 0

    # update predecessor and successor nodes
    p.h = q.i
    q.t = p.i

    # update associated vehicle
    v.n -= -1
    v.l -= n.q
    v.x = (x - n.x * n.q) / v.l
    v.y = (y - n.y * n.q) / v.l
    v.c -= δ
    
    # update solution
    s.c -= δ

    return s
end

function insertnode!(s::Solution, n::Node, p::Node, q::Node, v::Vehicle)
    # fetch features
    x = v.x * v.l
    y = v.y * v.l
    aₜᵢ = s.A[p.i, n.i]
    aᵢₕ = s.A[n.i, q.i]
    aₜₕ = s.A[p.i, q.i]
    δ = (aₜᵢ.l + aᵢₕ.l) - aₜₕ.l

    # update node
    n.t = p.i
    n.h = q.i
    n.v = v.i

    # update predecessor and successor nodes
    p.h = n.i
    q.t = n.i

    # update associated vehicle
    v.n += 1
    v.l += n.q
    v.x = (x + n.x * n.q) / v.l
    v.y = (y + n.y * n.q) / v.l
    v.c += δ

    # update solution
    s.c += δ

    return s
end