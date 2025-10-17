"""
    ALNSparameters

Optimization parameters for Adaptive Large Neighborhood Search (ALNS).

- j     :   Number of segments in the ALNS
- k     :   Number of segments to reset ALNS
- n     :   Number of iterations in an ALNS segment
- m     :   Number of iterations in a local search
- Ψᵣ    :   Tuple of removal operators
- Ψᵢ    :   Tuple of insertion operators
- Ψₗ    :   Tuple of local search operators
- σ₁    :   Score for a new best solution
- σ₂    :   Score for a new better solution
- σ₃    :   Score for a new worse but accepted solution
- μ̲     :   Minimum removal fraction
- e̲     :   Minimum nodes removed
- μ̅     :   Maximum removal fraction
- e̅     :   Maximum nodes removed
- ω̅     :   Initial temperature deviation parameter
- τ̅     :   Initial temperature probability parameter
- ω̲     :   Final temperature deviation parameter
- τ̲     :   Final temperature probability parameter
- θ     :   Cooling rate
- ρ     :   Reaction factor
"""
Base.@kwdef struct ALNSparameters{R<:Tuple, I<:Tuple, L<:Tuple}
    j::Int
    k::Int
    n::Int
    m::Int
    Ψᵣ::R
    Ψᵢ::I
    Ψₗ::L
    σ₁::Float64
    σ₂::Float64
    σ₃::Float64
    μ̲::Float64
    e̲::Int
    μ̅::Float64
    e̅::Int
    ω̅::Float64
    τ̅::Float64
    ω̲::Float64
    τ̲::Float64
    θ::Float64
    ρ::Float64
end