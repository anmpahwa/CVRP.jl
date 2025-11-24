using CVRP
using CSV
using Revise
using Random
using CPUTime
using DataFrames

let
    # Define instances
    I   = ["X-n101-k25", "X-n106-k14", "X-n110-k13", "X-n115-k10", "X-n120-k6", "X-n125-k30", "X-n129-k18", "X-n134-k13", "X-n139-k10", "X-n143-k7", "X-n148-k46", "X-n153-k22", "X-n157-k13", "X-n162-k11", "X-n167-k10", "X-n172-k51", "X-n176-k26", "X-n181-k23", "X-n186-k15", "X-n190-k8", "X-n195-k51"]
    II  = ["X-n200-k36", "X-n204-k19", "X-n209-k16", "X-n214-k11", "X-n219-k73", "X-n223-k34", "X-n228-k23", "X-n233-k16", "X-n237-k14", "X-n242-k48", "X-n247-k50", "X-n251-k28", "X-n256-k16", "X-n261-k13", "X-n266-k58", "X-n270-k35", "X-n275-k28", "X-n280-k17", "X-n284-k15", "X-n289-k60", "X-n294-k50", "X-n298-k31"]
    III = ["X-n303-k21", "X-n308-k13", "X-n313-k71", "X-n317-k53", "X-n322-k28", "X-n327-k20", "X-n331-k15", "X-n336-k84", "X-n344-k43", "X-n351-k40", "X-n359-k29", "X-n367-k17", "X-n376-k94", "X-n384-k52", "X-n393-k38"]
    IV  = ["X-n401-k29", "X-n411-k19", "X-n420-k130", "X-n429-k61", "X-n439-k37", "X-n449-k29", "X-n459-k26", "X-n469-k138", "X-n480-k70", "X-n491-k59"]
    V   = ["X-n502-k39", "X-n513-k21", "X-n524-k153", "X-n536-k96", "X-n548-k50", "X-n561-k42", "X-n573-k30", "X-n586-k159", "X-n599-k92"]
    VI  = ["X-n613-k62", "X-n627-k43", "X-n641-k35", "X-n655-k131", "X-n670-k130", "X-n685-k75"]
    VII = ["X-n701-k44", "X-n716-k35", "X-n733-k159", "X-n749-k98", "X-n766-k71", "X-n783-k48"]
    VIII= ["X-n801-k40", "X-n819-k171", "X-n837-k142", "X-n856-k95", "X-n876-k59", "X-n895-k37"]
    IX  = ["X-n916-k207", "X-n936-k151", "X-n957-k87", "X-n979-k58", "X-n1001-k43"]
    instances = [I..., II..., III..., IV...]
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
                j   =   100                     ,
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
                ω̅   =   0.01                    ,
                τ̅   =   0.5                     ,
                ω̲   =   0.001                   ,
                τ̲   =   0.05                    ,
                θ   =   exp(log(0.001 / 0.01 * log(1 / 0.5) / log(1 / 0.05)) / (0.50 * 100 * x))                  ,
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