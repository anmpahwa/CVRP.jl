function removenode!(n::Node, t::Node, h::Node, v::Vehicle, s::Solution)
    # fetch features
    G = s.G
    x = v.x * v.n
    y = v.y * v.n
    aₜₙ = G.A[t.i, n.i]
    aₙₕ = G.A[n.i, h.i]
    aₜₕ = G.A[t.i, h.i]
    # remove penalty
    s.p -= (v.l > v.q) ? (v.l - v.q) : 0
    # update node
    n.t = 0
    n.h = 0
    n.v = 0
    # update predecessor and successor nodes
    t.h = isdepot(t) ? t.h : h.i
    h.t = isdepot(h) ? h.t : t.i
    # update associated vehicle
    v.s = isdepot(t) ? h.i : v.s
    v.e = isdepot(h) ? t.i : v.e
    v.n -= 1
    v.l -= n.q
    v.x = iszero(v.n) ? 0. : (x - n.x) / v.n
    v.y = iszero(v.n) ? 0. : (y - n.y) / v.n
    # update solution cost
    s.c -= (aₜₙ.c + aₙₕ.c) - aₜₕ.c
    # add penalty
    s.p += (v.l > v.q) ? (v.l - v.q) : 0
    # return solution
    return s
end

function insertnode!(n::Node, t::Node, h::Node, v::Vehicle, s::Solution)
    # fetch features
    G = s.G
    x = v.x * v.n
    y = v.y * v.n
    aₜₙ = G.A[t.i, n.i]
    aₙₕ = G.A[n.i, h.i]
    aₜₕ = G.A[t.i, h.i]
    # remove penalty
    s.p -= (v.l > v.q) ? (v.l - v.q) : 0
    # update node
    n.t = t.i
    n.h = h.i
    n.v = v.i
    # update predecessor and successor nodes
    t.h = isdepot(t) ? t.h : n.i
    h.t = isdepot(h) ? h.t : n.i
    # update associated vehicle
    v.s = isdepot(t) ? n.i : v.s
    v.e = isdepot(h) ? n.i : v.e
    v.n += 1
    v.l += n.q
    v.x = (x + n.x) / v.n
    v.y = (y + n.y) / v.n
    # update solution cost
    s.c += (aₜₙ.c + aₙₕ.c) - aₜₕ.c
    # add penalty
    s.p += (v.l > v.q) ? (v.l - v.q) : 0
    # return solution
    return s
end