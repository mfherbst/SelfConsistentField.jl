"""Abstract type for SCF convergence accelerators"""
abstract type Accelerator end

"""
Use the knowledge of the previous iterates to refine
the current iterate for the next iteration
"""
function compute_next_iterate(acc::Union{Accelerator, Nothing}, iterate::ScfIterState)
    return iterate
end

"""
Postprocess an iterate, removing any implicit effect of the
accelerator (e.g. a level shifting in the eigenvalues or so.
"""
function postprocess_iterate(acc::Union{Accelerator, Nothing}, iterate::ScfIterState)
    return iterate
end

