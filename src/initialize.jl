function build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    
    df = CSV.read("$dir/$instance.vrp", DataFrame, silencewarnings=true)
    
    k₁ = findfirst(contains("DIMENSION"), df[:,1])
    k₂ = findfirst(contains("CAPACITY"), df[:,1])
    k₃ = findfirst(contains("NODE_COORD_SECTION"), df[:,1])
    k₄ = findfirst(contains("DEMAND_SECTION"), df[:,1])

    # Nodes
    n = parse(Int, df[k₁,2])
    N = Vector{Node}(undef, n)
    for i ∈ 1:n
        x = parse(Int, split(df[k₃+i,1])[2])
        y = parse(Int, split(df[k₃+i,1])[3])
        q = parse(Int, split(df[k₄+i,1])[2])
        N[i] = Node(i, x, y, q)
    end
    
    # Arcs
    A = Matrix{Arc}(undef, n, n)
    for i ∈ 1:n
        xᵢ = N[i].x
        yᵢ = N[i].y 
        for j ∈ 1:n
            xⱼ = N[j].x
            yⱼ = N[j].y 
            l  = ((xⱼ-xᵢ)^2 + (yⱼ-yᵢ)^2)^0.5
            a  = Arc(i, j, l)
            A[i,j] = a
        end
    end
            
    # Vehicles
    d = 0
    for i ∈ 1:n d += N[i].q end
    q = parse(Int, df[k₂,2])
    m = d ÷ q + 1
    V = Vector{Vehicle}(undef, m)
    for i ∈ 1:m V[i] = Vehicle(i, q) end

    # Graph
    G = (N, A, V)
    return G
end


function initialize(rng::AbstractRNG, instance::String; dir=joinpath(dirname(@__DIR__), "instances"))

end
