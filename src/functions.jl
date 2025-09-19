isdepot(n::Node) = isone(n.i)

iscustomer(n::Node) = !isdepot(n)

isopen(n::Node) = iszero(n.v)

isclose(n::Node) = !isopen(n)

function relatedness(m::Node, n::Node)
    z = sqrt((m.x - n.x)^2 + (m.y - n.y)^2) * abs(m.q - n.q)
    r = 1 / (z + 1)
    return r
end

isopt(v::Vehicle) = !iszero(v.n)

function relatedness(u::Vehicle, v::Vehicle)
    z = sqrt((u.x - v.x)^2 + (u.y - v.y)^2) * abs(u.l - v.l) * abs(u.n - v.n)
    r = 1 / (z + 1)
    return r
end

function vectorize(s::Solution)
    Z = Int[]
    for v ∈ s.V
        push!(Z, 1)
        i = v.s
        for _ ∈ 1:v.n
            push!(Z, i)
            n = s.N[i]
            i = n.h
        end
        push!(Z, 1)
    end
    return Z
end

function isfeasible(s::Solution)
    for v ∈ s.V if v.l > v.q return false end end
    return true
end

f(s::Solution) = s.c + s.p * 10 ^ ceil(log10(s.c))

h(s::Solution) = hash(vectorize(s))