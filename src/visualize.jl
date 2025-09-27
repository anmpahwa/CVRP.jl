function visualize(instance::String; backend=gr)
    backend()
    G = build(instance)
    N = G.N
    fig = plot(legend=:none)
    K = lastindex(N)
    X = zeros(Float64, K)       # abcissa
    Y = zeros(Float64, K)       # ordinate
    C = fill("color", K)        # color
    S = zeros(Int, K)           # size
    M = fill(:shape, K)         # marker
    for (k,n) ∈ enumerate(N)
        X[k] = n.x
        Y[k] = n.y
        C[k] = isdepot(n) ? "#b4464b" : "#d1e0ec"
        S[k] = isdepot(n) ? 6 : 5
        M[k] = isdepot(n) ? :rect : :circle
    end
    scatter!(X, Y, color=C, markersize=S, markershape=M, markerstrokewidth=0)
    return fig
end

function visualize(s::Solution; backend=gr)
    backend()
    fig = plot(legend=:none)
    # nodes
    G = s.G
    N = G.N
    K = lastindex(N)
    X = zeros(Float64, K)       # abcissa
    Y = zeros(Float64, K)       # ordinate
    C = fill("color", K)        # color
    S = zeros(Int, K)           # size
    M = fill(:shape, K)         # marker
    for (k,n) ∈ enumerate(N)
        X[k] = n.x
        Y[k] = n.y
        C[k] = isdepot(n) ? "#b4464b" : isopen(n) ? "#d1e0ec" : "#4682b4"
        S[k] = isdepot(n) ? 6 : 5
        M[k] = isdepot(n) ? :rect : :circle
    end
    scatter!(X, Y, color=C, markersize=S, markershape=M, markerstrokewidth=0)
    # arcs
    Z = vectorize(s)
    K = lastindex(Z)
    X = [N[Z[k]].x for k ∈ 1:K]    # abcissa    
    Y = [N[Z[k]].y for k ∈ 1:K]    # ordinate
    plot!(X, Y, color="#23415a")
    # figure
    return fig
end