"""
Struct to track convergence of the SCF procedure
"""
mutable struct ScfConvergence
    error_norm::Float64          # Frobenius error norm of the Pulay error
    energy_change::Dict{String, Float64}   # Change in the energy terms
    is_converged::Bool           # Is SCF converged
end
