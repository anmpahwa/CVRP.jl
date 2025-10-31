using CVRP
using Revise
using Random
using CPUTime
using DataFrames

let
    # Define instances
    instances = ["X-n101-k25", "X-n200-k36", "X-n308-k13"]
    # Define random number generator seeds
    seeds = [1010, 1106, 1509, 1604, 1905, 2104, 2412, 2703, 2710, 2807]
    # Dataframes to store solution quality and run time
    dfᶠ = DataFrame([instances, zeros(length(instances)), [zeros(length(instances)) for _ ∈ seeds]...], ["instance", "initial", ["$seed" for seed ∈ seeds]...])
    dfᵗ = DataFrame([instances, zeros(length(instances)), [zeros(length(instances)) for _ ∈ seeds]...], ["instance", "initial", ["$seed" for seed ∈ seeds]...])
    for (i,instance) ∈ enumerate(instances)
        # Visualize instance
        display(visualize("benchmark/$instance"))
        # Define inital solution method and build the initial solution
        t₁ = @elapsed s₁ = initialize("benchmark/$instance"; method=:dynamic);
        # Visualize initial solution
        display(visualize(s₁))
        # Fetch solution characteristics
        println("\nINSTANCE: $instance")
        println("Initial Solution:")
        println("   Run Time: $(round(t₁, digits=3)) seconds")
        println("   Feasiblity: $(isfeasible(s₁)) | $(round(s₁.p, digits=3))")
        println("   Objective Function: $(round(f(s₁), digits=3))")
        # Store results
        dfᶠ[i,2] = f(s₁)
        dfᵗ[i,2] = t₁
        # Save results
        CSV.write("objective_function.csv", dfᶠ)
        CSV.write("run_time.csv", dfᵗ)
        for (j,seed) ∈ enumerate(seeds)
            println("\nOptimal Solution | seed: $seed")
            rng = MersenneTwister(seed);
            # Define ALNS parameters
            x = max(100, lastindex(s₁.G.N))
            χ = ALNSparameters(
                j   =   50                      ,
                k   =   5                       ,
                n   =   x                       ,
                m   =   100x                    ,
                Ψᵣ  =   (
                            randomnode!         ,
                            randomarc!          ,
                            randomvehicle!      ,
                            relatednode!        ,
                            relatedarc!         ,
                            relatedvehicle!     ,
                            worstnode!          ,
                            worstarc!           ,
                            worstvehicle!
                        )                       ,
                Ψᵢ  =   (
                            bestprecise!        ,
                            bestperturb!        ,
                            greedyprecise!      ,
                            greedyperturb!      ,
                            regret2precise!     ,
                            regret2perturb!     ,
                            regret3precise!     ,
                            regret3perturb!
                        )                       ,
                Ψₗ  =   (
                            intramove!          ,
                            intermove!          ,
                            intraswap!          ,
                            interswap!          ,
                            intraopt!           ,
                            interopt!
                        )                       ,
                σ₁  =   15.0                    ,
                σ₂  =   10.0                    ,
                σ₃  =   3.0                     ,
                μ̲   =   0.05                    ,
                e̲   =   2                       ,
                μ̅   =   0.2                     ,
                e̅   =   30                      ,
                ω̅   =   0.05                    ,
                τ̅   =   0.5                     ,
                ω̲   =   0.01                    ,
                τ̲   =   0.01                    ,
                θ   =   0.9993                  ,
                ρ   =   0.1
            );
            # Run ALNS and fetch best solution
            t₂ = @elapsed s₂ = ALNS(rng, χ, s₁);
            # Visualize best solution
            display(visualize(s₂))
            # Fetch solution characteristics
            println("   Run Time: $(round(t₂, digits=3)) seconds")
            println("   Feasiblity: $(isfeasible(s₂)) | $(round(s₂.p, digits=3))")
            println("   Objective Function: $(round(f(s₂), digits=3))")
            # Store results
            dfᶠ[i,j+2] = f(s₂)
            dfᵗ[i,j+2] = t₂
            # Save results
            CSV.write("objective_function.csv", dfᶠ)
            CSV.write("run_time.csv", dfᵗ)
        end
        println(dfᶠ)
        println(dfᵗ)
    end
    return
end 