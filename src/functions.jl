function vectorize(s::Solution)
    Z = Int[]
    for v ∈ s.V
        push!(Z, 1)
        i = v.s
        for k ∈ 1:v.n
            push!(Z, s.N[i])
            n = s.N[i]
            i = n.h
        end
        push!(Z, 1)
    end
    return Z
end

function relatedness(m::Node, n::Node)
    z = sqrt((m.x - n.x)^2 + (m.y - n.y)^2) * abs(m.q - n.q)
    r = 1 / (z + 1)
    return r
end

function relatedness(u::Vehicle, v::Vehicle)
    z = sqrt((u.x - v.x)^2 + (u.y - v.y)^2) * abs(u.l - v.l) * abs(u.n - v.n)
    r = 1 / (z + 1)
    return r
end

f(s::Solution) = s.c
h(s::Solution) = hash(vectorize(s))