"""
    Computes next iterate using DIIS
"""
function compute_next_iterate(acc::cDIIS, iterstate::ScfIterState)
    # Push iterate and error to state
    map(σ -> push_iterate!(acc.state[σ], spin(iterstate.fock, σ), spin(iterstate.error_pulay, σ)), spinloop(iterstate.fock))

    # Check if the number of known fock and error matrices is equal for both
    # spins before doing anything
    history_size = acc.state[1].iterate.length
    for σ in spinloop(iterstate.fock)
        @assert acc.state[σ].iterate.length == history_size
        @assert acc.state[σ].error.length == history_size
    end

    # To save memory we store only new_iterate once and pass views of it to the
    # computation routines that write directly into the view.
    new_iterate = zeros(size(iterstate.fock))

    # Defining anonymous functions with given arguments improves readability later on.
    matrix(σ) = diis_build_matrix(acc, acc.state[σ])
    coefficients(A) = diis_solve_coefficients(acc, A, acc.conditioning_threshold)
    compute(c, σ) = compute_next_iterate_matrix!(acc.state[σ], c, spin(new_iterate, σ), acc.coefficient_threshold)

    # If sync_spins is enabled, we need to calculate the coefficients using the
    # merged matrix. This also means we need to remove the same number of
    # matrices from both histrories.
    if acc.sync_spins & (spincount(iterstate.fock) == 2)
        A = merge_matrices(acc, matrix(1), matrix(2))
        c, matrixpurgecount = coefficients(A)
        map(σ -> compute(c, σ), spinloop(iterstate.fock))
        map(σ -> purge_from_state(acc.state[σ], matrixpurgecount), spinloop(iterstate.fock))
    else
        # If we are calculating the spins separately, each spin has its own coefficients.
        for σ in spinloop(iterstate.fock)
            c, matrixpurgecount = coefficients(matrix(σ))
            compute(c, σ)
            purge_from_state(acc.state[σ], matrixpurgecount)
        end
    end
    return update_iterate_matrix(iterstate, new_iterate)
end

"""
    Computes a new matrix to be used as Fock Matrix in next iteration
    The function writes the result into the passed argument 'iterate'
"""
function compute_next_iterate_matrix!(state::DiisState, c::AbstractArray, iterate::SubArray, coefficient_threshold::Float64)
    # add very small coefficients to the largest one but always use the most
    # recent iterate matrix regardless of the coefficient value
    mask = map(x -> norm(x) > coefficient_threshold, c)
    mask[1] = true
    c[argmax(c)] += sum(c[ .! mask])

    # Construct new Fock Matrix using obtained coefficients
    # and write it to the given iterate matrix. We assume, that
    # iterate is a matrix of zeros.
    for i in eachindex(c)[mask]
        iterate .+= c[i] * state.iterate[i]
    end
end

