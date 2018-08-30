"""
    Computes next iterate using DIIS
"""
# TODO rename acc to algorithm
function compute_next_iterate(acc::Union{cDIIS,EDIIS}, state::ScfIterState)
    iterate = get_iterate_matrix(state)

    # Push iterate and error to history
    map(σ -> push_iterate!(acc, acc.history[σ], spin(state.fock, σ), spin(state.error_pulay, σ), spin(state.density, σ), state.energies), spinloop(iterate))

    # To save memory we store only new_iterate once and pass views of it to the
    # computation routines that write directly into the view.
    new_iterate = zeros(size(iterate))

    # Defining anonymous functions with given arguments improves readability later on.
    matrix(σ) = diis_build_matrix(acc, acc.history[σ])
    coefficients(A) = diis_solve_coefficients(acc, A)
    compute(c, σ) = compute_next_iterate_matrix!(acc.history[σ], c, spin(new_iterate, σ), acc.coefficient_threshold)

    # If sync_spins is enabled, we need to calculate the coefficients using the
    # merged matrix. This also means we need to remove the same number of
    # matrices from both histrories.
    if (:sync_spins ∈ fieldnames(typeof(acc))) ? (acc.sync_spins & (spincount(iterate) == 2)) : false
        A = merge_matrices(acc, matrix(1), matrix(2))
        c, matrixpurgecount = coefficients(A)
        map(σ -> compute(c, σ), spinloop(iterate))
        map(σ -> purge_from_history(acc, acc.history[σ], matrixpurgecount + 1), spinloop(iterate)) # +1, since we have already pushed a new matrix!
    else
        # If we are calculating the spins separately, each spin has its own coefficients.
        for σ in spinloop(iterate)
            c, matrixpurgecount = coefficients(matrix(σ))
            compute(c, σ)
            purge_from_history(acc, acc.history[σ], matrixpurgecount + 1) # +1, since we have already pushed a new matrix!
        end
    end
    return update_iterate_matrix(state, new_iterate)
end

"""
    Computes a new matrix to be used as Fock Matrix in next iteration
    The function writes the result into the passed argument 'iterate'
"""
function compute_next_iterate_matrix!(history::DiisHistory, c::AbstractArray, iterate::SubArray, coefficient_threshold::Float64)
    # add very small coefficients to the largest one but always use the most
    # recent iterate matrix regardless of the coefficient value
    mask = map(x -> norm(x) >= coefficient_threshold, c)
    mask[1] = true
    c[argmax(c)] += sum(c[ .! mask])

    # Construct new Fock Matrix using obtained coefficients
    # and write it to the given iterate matrix. We assume, that
    # iterate is a matrix of zeros.
    for i in eachindex(c)[mask]
        iterate .+= c[i] * history.iterate[i]
    end
end

