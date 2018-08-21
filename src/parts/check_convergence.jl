"""
Struct to track convergence of the SCF procedure
"""
struct ScfConvergence
    error_norm::Float64          # Frobenius error norm of the Pulay error
    energy_change::Dict{String, Float64}   # Change in the energy terms
    is_converged::Bool           # Is SCF converged
end

"""
Check whether an SCF is converged
"""
function check_convergence(olditerate::ScfIterState, newiterate::ScfIterState;
                           max_error_norm=5e-7,
                           max_energy_total_change=1.25e-07,
                           max_energy_1e_change=5e-5,
                           kwargs...)
    error_norm = norm(reshape(newiterate.error_pulay, :))
    energy_change = Dict{String, Float64}()
    for (key, oval) in olditerate.energies
        nval = newiterate.energies[key]
        val = nval - oval
        energy_change[key] = val
    end

    converged = true
    if error_norm >= max_error_norm
        converged = false
    end
    if energy_change["total"] >= max_energy_total_change
        converged = false
    end
    if energy_change["1e"] >= max_energy_1e_change
        converged = false
    end
    return ScfConvergence(error_norm, energy_change, converged)
end
