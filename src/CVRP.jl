module CVRP

using CSV
using Plots
using Random
using Revise
using DataFrames

include("datastructure.jl")
include("functions.jl")
include("initialize.jl")
include("operations.jl")
include("remove.jl")
include("insert.jl")
include("localsearch.jl")
include("parameters.jl")
include("ALNS.jl")
include("visualize.jl")

export  build, 
        vectorize, f, 
        ALNSparameters, ALNS, 
        visualize

end
