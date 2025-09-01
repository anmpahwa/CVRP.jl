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