#
# Direct inversion in the iterative subspace
#

"""
    cDIIS
"""
struct cDIIS <: Algorithm
    n_diis_size::Int
    sync_spins::Bool
    conditioning_threshold::Float64
    coefficient_threshold::Float64
    history::Tuple{DiisHistory, DiisHistory}
end

function cDIIS(problem::ScfProblem, state::ScfIterState, lg::Logger;
               n_diis_size = 5, sync_spins = true,
               conditioning_threshold = 1e-14, coefficient_threshold = 10e-6,
               params...)

    log!(lg, "setting up cDIIS", :info, :cdiis, :setup)
    historyα = DiisHistory(n_diis_size, need_error = true)
    if spincount(get_iterate_matrix(state)) == 2
        history = (historyα, DiisHistory(n_diis_size, need_error = true))
    else
        history = (historyα, historyα)
    end

    cDIIS(n_diis_size, sync_spins, conditioning_threshold, coefficient_threshold, history)
end

function iterate(cdiis::cDIIS, rp::SubReport)
    lg = Logger(rp)

    # Push iterate and error to history
    log!(lg, "pushing new iterate to history", :debug, :diis)
    history = push_iterate(cdiis.history, rp.state)


    iterate_computation = compute_next_iterate
    diis_matrix_formula(i,j) = sum(σ -> tr(history[σ].error[i]' * history[σ].error[j]), spinloop(rp.state))

    if !cdiis.sync_spins && (spincount(iterate) == 2)
        iterate_computation = compute_next_iterate_separating_spins
        diis_matrix_formula = (diis_matrix_formula, diis_matrix_formula)
    end

    iterate, new_history, matrixpurgecount = iterate_computation(get_iterate_matrix(rp.state),
                                                history, diis_matrix_formula,
                                                compute_cdiis_coefficients, lg,
                                                coefficient_threshold = cdiis.coefficient_threshold,
                                                conditioning_threshold = cdiis.conditioning_threshold)

    for σ in spinloop(iterate)
        # The next iteraton removes automatically removes one entry from the
        # history. This means we need to purge an additional matrix to actually
        # have a removal effect. Hence +1
        if matrixpurgecount[σ] > 0
            was_at_capacity = length(cdiis.history[σ].iterate) == cdiis.history[σ].iterate.capacity
            compensation = was_at_capacity ? 1 : 0
            purge_from_history!(new_history[σ], matrixpurgecount[σ] + compensation)
        end
    end

    state = update_iterate_matrix(rp.state, iterate)
    new_cdiis = cDIIS(cdiis.n_diis_size, cdiis.sync_spins,
                      cdiis.conditioning_threshold, cdiis.coefficient_threshold, new_history)

    return new_cdiis, new_subreport(new_cdiis, state, lg, rp)
end

"""
    Solves the linear system after removing small eigenvalues to improve
    consistency.
"""
function compute_cdiis_coefficients(B::AbstractArray, lg::Logger; conditioning_threshold::Float64, params...)
    A = build_diis_linear_system_matrix(B)
    # Right hand side of the equation
    rhs = cdiis_build_rhs(size(A, 1))

    # calculate the eigenvalues of A and select sufficiently large eigenvalues
    λ, U = eigen(A)
    mask = map(x -> norm(x) > conditioning_threshold, λ)

    if !all(mask)
        log!(lg, @sprintf("Removing %i of %i eigenvalues from DIIS linear system.", count(.! mask), length(mask)),
             :cdiis, :info)
    end

    # if all eigenvalues are under the threshold, we cannot calculate sane
    # coefficients. The current fock matrix should be used without
    # modifications.
    if all( .! mask)
        log!(lg, "All eigenvalues are under the threshold! Skipping iterate modification…", :cdiis, :info)
        c = zeros(size(A, 1) - 1)
        c[1] = 1
        return c
    end

    # Obtain the solution of the linear system A * c = rhs using the inverse
    # matrix of A constructed using the above decomposition
    c = U[:,mask] * Diagonal(1 ./ λ[mask]) * U[:,mask]' * rhs

    # Note that c has size (n_diis_size + 1) since the last element is the
    # lagrange multiplier corresponding to the constraint
    # \sum_{i=1}^n_diis_size c_i = 1 We need to remove this element!
    return c[1:length(c) - 1], count(.! mask)
end

"""
    Right hand side of the cDIIS linear system.
    This is a vector of size m+1 containing only ones
    except for the last element which is zero.
"""
function cdiis_build_rhs(vectorsize::Int)
    rhs = zeros(vectorsize)
    rhs[end] = 1
    return rhs
end
