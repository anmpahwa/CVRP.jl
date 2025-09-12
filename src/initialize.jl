function build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # read instance file
    df = CSV.read("$dir/$instance.vrp", DataFrame, silencewarnings=true)
    # fetch key indices
    k₁ = findfirst(contains("DIMENSION"), df[:,1])::Int
    k₂ = findfirst(contains("CAPACITY"), df[:,1])::Int
    k₃ = findfirst(contains("NODE_COORD_SECTION"), df[:,1])::Int
    k₄ = findfirst(contains("DEMAND_SECTION"), df[:,1])::Int
    # fetch nodes
    n = parse(Int, df[k₁,2])
    N = Vector{Node}(undef, n)
    for i ∈ 1:n
        x = parse(Int, split(df[k₃+i,1])[2])
        y = parse(Int, split(df[k₃+i,1])[3])
        q = parse(Int, split(df[k₄+i,1])[2])
        N[i] = Node(i, x, y, q)
    end
    # create arcs
    A = Matrix{Arc}(undef, n, n)
    for i ∈ 1:n
        xᵢ = N[i].x
        yᵢ = N[i].y 
        for j ∈ 1:n
            xⱼ = N[j].x
            yⱼ = N[j].y 
            c  = ((xⱼ-xᵢ)^2 + (yⱼ-yᵢ)^2)^0.5
            a  = Arc(i, j, c)
            A[i,j] = a
        end
    end   
    # fetch vehicles
    q = parse(Int, df[k₂,2])
    V = Vector{Vehicle}(undef, n-1)
    for i ∈ 1:(n-1) V[i] = Vehicle(i, q) end
    # create graph
    G = (N, A, V)
    return G
end


function initialize(instance::String, dir=joinpath(dirname(@__DIR__), "instances"))
    # pre-initialize
    G = build(instance; dir=dir)
    s = Solution(G...)
    N = s.N
    V = s.V
    # initialize
    K = eachindex(V)
    C = fill(-Inf, (K,K))       # C[i,j]: Savings from merging route for vehicle V[i] and V[j]
    ϕ = ones(Int, K)            # ϕ[i]  : binary weight for vehicle V[i]
    m = 0.
    for k ∈ K
        n = N[k+1]
        t = N[1]
        h = N[1]
        v = V[k]
        m += N[k+1].q / V[begin].q
        insertnode!(n, t, h, v, s)
    end
    while sum(isopt.(V)) > ceil(Int, m)
        z = f(s)
        # iterate through each vehicle pair
        for i ∈ K
            vᵢ = V[i]
            if !isopt(vᵢ) continue end
            for j ∈ K
                vⱼ = V[j]
                k  = vⱼ.n
                if !isopt(vⱼ) continue end
                if isequal(i,j) continue end
                if iszero(ϕ[i]) & iszero(ϕ[j]) continue end
                # merge vehicle route V[j] into V[i]
                for _ ∈ 1:k
                    d = N[1]
                    n = N[vⱼ.s]
                    t = N[vᵢ.e]
                    h = N[n.h]
                    removenode!(n, d, h, vⱼ, s)
                    insertnode!(n, t, d, vᵢ, s)
                end
                # compute savings
                z⁻ = f(s)
                Δ  = z - z⁻
                C[i,j] = Δ
                # unmerge vehicle routes
                for _ ∈ 1:k
                    d = N[1]
                    n = N[vᵢ.e]
                    t = N[n.t]
                    h = N[vⱼ.s]
                    removenode!(n, t, d, vᵢ, s)
                    insertnode!(n, d, h, vⱼ, s)
                end
            end
        end
        # merge vehicle routes with highest savings
        i,j = Tuple(argmax(C))
        vᵢ = V[i]
        vⱼ = V[j]
        k  = vⱼ.n
        for _ ∈ 1:k
            d = N[1]
            n = N[vⱼ.s]
            t = N[vᵢ.e]
            h = N[n.h]
            removenode!(n, d, h, vⱼ, s)
            insertnode!(n, t, d, vᵢ, s)
        end
        # update vectors
        C[i, :] .= -Inf
        C[:, i] .= -Inf
        C[j, :] .= -Inf
        C[:, j] .= -Inf
        ϕ .= 0
        ϕ[i] = 1
        ϕ[j] = 1
    end
    # remove non-operational vehicles
    filter!(isopt, V)
    # reset indices
    K = eachindex(V)
    for k ∈ K
        v = V[k]
        v.i = k
        n = N[v.s]
        for _ ∈ 1:v.n
            n.v = k
            n = N[n.h]
        end
    end
    # return solution
    return s
end
