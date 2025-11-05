"""
    visualize(instance::String; backend=gr)

Returns plot visualizing the `instance` using the given `backend` library.

Note: the plot illustrates the customer nodes and the depot node.
"""
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

"""
    visualize(s::Solution; backend=gr)

Returns plot visualizing the solution `s` using the given `backend` library.

Note: the plot illustrates the customer nodes, depot node, and arcs traversed
by the vehicle in the solution.
"""
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

"""
    pltcnv(X::Vector{UInt64}, Z::Vector{Float64}; backend=gr)

Plots convergence using objective function evaluations vector `Z`. 
Uses given `backend` to plot (defaults to `gr`).
"""
function pltcnv(X::Vector{UInt64},Z::Vector{Float64}; backend=gr)
    backend()
    fig = plot(legend=:none)
    xₒ = X[1]
    zₒ = Z[1]
    z⃰  = minimum(Z)
    I  = eachindex(Z)
    Y₁ = []             # improvements
    Y₂ = zeros(I)       # convergence
    for i ∈ I
        x = X[i]
        z = Z[i]
        if !isequal(x, xₒ) && isless(z, zₒ)
            xₒ = x
            zₒ = z
            push!(Y₁, i)
        end
        Y₂[i] = (z/z⃰ - 1) * 100
    end
    vline!(Y₁, color=:black, linewidth=0.25)
    plot!(I, Y₂, xlabel="iterations", ylabel="deviation from the best (%)", color=:steelblue)
    # figure
    return fig
end