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

function cDIIS(problem::ScfProblem, iterate::Iterate, lg::Logger;
               n_diis_size = 5, sync_spins = true,
               conditioning_threshold = 1e-14, coefficient_threshold = 10e-6,
               params...)

    log!(lg, "setting up cDIIS", :info, :cdiis, :setup)
    historyα = DiisHistory(n_diis_size, need_error = true)
    if !sync_spins && spincount(iterate.fock) == 2
        history = (historyα, DiisHistory(n_diis_size, need_error = true))
    else
        history = (historyα, historyα)
    end

    cDIIS(n_diis_size, sync_spins, conditioning_threshold, coefficient_threshold, history)
end

function notify(cdiis::cDIIS, rp::StepState)
    history = push_iterate(cdiis.history, rp.iterate)
    diis_matrix_formula(i,j) = sum(σ -> tr(history[σ].error[i]' * history[σ].error[j]), spinloop(rp.iterate))

    if !cdiis.sync_spins && spincount(rp.iterate) == 2
        historyα_with_diis_matrix_entries = compute_diis_matrix(diis_matrix_formula, history[1], false)
        historyβ_with_diis_matrix_entries = compute_diis_matrix(diis_matrix_formula, history[2], false)
        new_history = (historyα_with_diis_matrix_entries, historyβ_with_diis_matrix_entries)
    else
        historyα_with_diis_matrix_entries = compute_diis_matrix(diis_matrix_formula, history[1], false)
        new_history = (historyα_with_diis_matrix_entries, history[2])
    end

    new_cdiis = cDIIS(cdiis.n_diis_size, cdiis.sync_spins, cdiis.conditioning_threshold, cdiis.coefficient_threshold, new_history)
    return new_cdiis, StepState(new_cdiis, rp)
end

function iterate(cdiis::cDIIS, rp::StepState)
    lg = Logger(rp)

    # Push iterate and error to history
    log!(lg, "pushing new iterate to history", :debug, :diis)
    history = push_iterate(cdiis.history, rp.iterate)


    fock_computation = compute_next_fock
    diis_matrix_formula(i,j) = sum(σ -> tr(history[σ].error[i]' * history[σ].error[j]), spinloop(rp.iterate))

    if !cdiis.sync_spins && (spincount(fock) == 2)
        fock_computation = compute_next_fock_separating_spins
        diis_matrix_formula = (diis_matrix_formula, diis_matrix_formula)
    end

    fock, new_history, matrixpurgecount = fock_computation(rp.iterate.fock,
                                                history, diis_matrix_formula,
                                                compute_cdiis_coefficients, lg,
                                                coefficient_threshold = cdiis.coefficient_threshold,
                                                conditioning_threshold = cdiis.conditioning_threshold)

    for σ in spinloop(fock)
        # The next iteraton removes automatically removes one entry from the
        # history. This means we need to purge an additional matrix to actually
        # have a removal effect. Hence +1
        if matrixpurgecount[σ] > 0
            was_at_capacity = length(cdiis.history[σ].fock) == cdiis.history[σ].fock.capacity
            compensation = was_at_capacity ? 1 : 0
            purge_from_history!(new_history[σ], matrixpurgecount[σ] + compensation)
        end
    end

    iterate = update_fock_matrix(rp.iterate, fock)
    new_cdiis = cDIIS(cdiis.n_diis_size, cdiis.sync_spins,
                      cdiis.conditioning_threshold, cdiis.coefficient_threshold, new_history)

    return new_cdiis, StepState(new_cdiis, iterate, lg, rp)
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
        log!(lg, "All eigenvalues are under the threshold! Skipping fock modification…", :cdiis, :info)
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
