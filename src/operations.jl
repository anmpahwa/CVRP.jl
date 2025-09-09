function removenode!(n::Node, t::Node, h::Node, v::Vehicle, s::Solution)
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
    t.h = isdepot(t) ? 0 : h.i
    h.t = isdepot(h) ? 0 : t.i
    # update associated vehicle
    v.s = isdepot(t) ? h.i : v.s
    v.e = isdepot(h) ? t.i : v.s
    v.n -= -1
    v.l -= n.q
    v.x = (x - n.x) / v.n
    v.y = (y - n.y) / v.n
    v.c -= δ
    # update solution
    s.c -= δ
    return s
end

function insertnode!(n::Node, t::Node, h::Node, v::Vehicle, s::Solution)
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
    t.h = isdepot(t) ? 0 : n.i
    h.t = isdepot(h) ? 0 : n.i
    # update associated vehicle
    v.s = isdepot(t) ? n.i : v.s
    v.e = isdepot(h) ? n.i : v.e
    v.n += 1
    v.l += n.q
    v.x = (x + n.x) / v.n
    v.y = (y + n.y) / v.n
    v.c += δ
    # update solution
    s.c += δ
    return s
end