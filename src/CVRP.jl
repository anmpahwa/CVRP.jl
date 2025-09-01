module CVRP

using CSV
using Plots
using Random
using Revise
using DataFrames

include("datastructure.jl")
include("functions.jl")
include("initialize.jl")
include("visualize.jl")

export build, initialize, vectorize, visualize

end
