"""
    EDIIS
"""
mutable struct EDIIS <: Algorithm
    n_diis_size::Int
    coefficient_threshold::Float64
    history::Tuple{DiisHistory, DiisHistory}
end

function EDIIS(problem::ScfProblem, state::ScfIterState, lg::Logger;
               n_diis_size = 5, coefficient_threshold = 10e-6,
               params...)

    log!(lg, "setting up EDIIS", :info, :ediis, :setup)
    historyα = DiisHistory(n_diis_size, need_density = true, need_energies = true)
    if spincount(get_iterate_matrix(state)) == 2
        history = (historyα, DiisHistory(n_diis_size, need_density = true, need_energies = true))
    else
        history = (historyα, historyα)
    end

    EDIIS(n_diis_size, coefficient_threshold, history)
end

function notify(ediis::EDIIS, rp::StepState)
    history = push_iterate(ediis.history, rp.state)
    diis_matrix_formula(i,j) = sum(σ -> tr((history[σ].iterate[i] - history[σ].iterate[j]) *
                                  (history[σ].density[i] - history[σ].density[j])),
                                   spinloop(rp.state))

    historyα_with_diis_matrix_entries = compute_diis_matrix(diis_matrix_formula, history[1], false)
    new_history = (historyα_with_diis_matrix_entries, history[2])

    new_ediis = EDIIS(ediis.n_diis_size, ediis.coefficient_threshold, new_history)
    return new_ediis, StepState(new_ediis, rp)
end

function iterate(ediis::EDIIS, rp::StepState)
    lg = Logger(rp)

    # Push iterate and error to history
    log!(lg, "pushing new iterate to history", :debug, :diis)
    history = push_iterate(ediis.history, rp.state)


    diis_matrix_formula(i,j) = sum(σ -> tr((history[σ].iterate[i] - history[σ].iterate[j]) *
                                  (history[σ].density[i] - history[σ].density[j])),
                                   spinloop(rp.state))

    iterate, new_history, matrixpurgecount = compute_next_iterate(get_iterate_matrix(rp.state),
                                                history, diis_matrix_formula,
                                                compute_ediis_coefficients, lg,
                                                coefficient_threshold = ediis.coefficient_threshold,
                                                energies = history[1].energies)

    for σ in spinloop(iterate)
        # The next iteraton removes automatically removes one entry from the
        # history. This means we need to purge an additional matrix to actually
        # have a removal effect. Hence +1
        if matrixpurgecount[σ] > 0
            was_at_capacity = length(ediis.history[σ]) == ediis.history[σ].capacity
            compensation = was_at_capacity ? 1 : 0
            purge_from_history!(new_history[σ], matrixpurgecount[σ] + compensation)
        end
    end

    state = update_iterate_matrix(rp.state, iterate)
    new_ediis = EDIIS(ediis.n_diis_size, ediis.coefficient_threshold, new_history)

    return new_ediis, StepState(new_ediis, state, lg, rp)
end

"""
    Calculates coefficients
"""
function compute_ediis_coefficients(B::AbstractArray, lg::Logger; energies::CircularBuffer{Dict}, params...)
    m = size(B,1)
    # TODO put the energies somewhere else
    E = map(energies -> energies["total"], energies)

    function f(x)
        c = x.^2/sum(x.^2)
        E'*c - 1/2 * c'*B*c
    end

    guess = ones(m) / m
    guess[1] = 1

    options = Optim.Options(x_tol = 10e-5)
    res = optimize(f, guess, BFGS(), options; inplace = false, autodiff = :forward)
    x = Optim.minimizer(res)
    c = x.^2/sum(x.^2)

    # If number of iterations in optimization is zero, reuse old matrix
    if Optim.iterations(res) == 0
        c = zeros(m)
        c[1] = 1
    end

    return c, 0
end
