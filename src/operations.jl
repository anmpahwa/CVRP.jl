function removenode!(s::Solution, n::Node, t::Node, h::Node, v::Vehicle)
    # fetch features
    x = v.x * v.n
    y = v.y * v.n
    aₜₙ = s.A[t.i, n.i]
    aₙₕ = s.A[n.i, h.i]
    aₜₕ = s.A[t.i, h.i]
    δ = (aₜₙ.l + aₙₕ.l) - aₜₕ.l
    # update node
    n.t = 0
    n.h = 0
    n.v = 0
    # update predecessor and successor nodes
    t.h = h.i
    h.t = t.i
    # update associated vehicle
    v.n -= -1
    v.l -= n.q
    v.x = (x - n.x) / v.n
    v.y = (y - n.y) / v.n
    v.c -= δ
    # update solution
    s.c -= δ
    return s
end

function insertnode!(s::Solution, n::Node, t::Node, h::Node, v::Vehicle)
    # fetch features
    x = v.x * v.n
    y = v.y * v.n
    aₜₙ = s.A[t.i, n.i]
    aₙₕ = s.A[n.i, h.i]
    aₜₕ = s.A[t.i, h.i]
    δ = (aₜₙ.l + aₙₕ.l) - aₜₕ.l
    # update node
    n.t = t.i
    n.h = h.i
    n.v = v.i
    # update predecessor and successor nodes
    t.h = n.i
    h.t = n.i
    # update associated vehicle
    v.n += 1
    v.l += n.q
    v.x = (x + n.x) / v.n
    v.y = (y + n.y) / v.n
    v.c += δ
    # update solution
    s.c += δ
    return s
end