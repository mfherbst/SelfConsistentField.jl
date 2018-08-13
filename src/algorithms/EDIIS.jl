"""
    EDIIS
"""
mutable struct EDIIS <: Accelerator
    state::Tuple{DiisState, DiisState}
    sync_spins::Bool
    conditioning_threshold::Float64
    coefficient_threshold::Float64

    function EDIIS(problem::ScfProblem; n_diis_size = 5, sync_spins = true, conditioning_threshold = 1e-14, coefficient_threshold = 1e-6, kwargs...)
        stateα = DiisState(n_diis_size)
        stateβ = DiisState(n_diis_size)
        new((stateα, stateβ), sync_spins, conditioning_threshold, coefficient_threshold)
    end
end

function diis_solve_coefficients(::EDIIS, A::AbstractArray, threshold::Float64)
    
    return c, 0
end

"""
    Linear System Matrix for the EDIIS accelerator.
    This is a symmetric matrix containing error overlaps B
    and ones in the form
    
    A = B  1
        1† 0
"""
function diis_build_matrix(::EDiis, state::DiisState)
    @assert state.n_diis_size > 0
    @assert state.iterate.length > 0


    return Symmetric(A)
end
