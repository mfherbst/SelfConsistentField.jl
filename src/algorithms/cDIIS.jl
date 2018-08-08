#
# Direct inversion in the iterative subspace
#

"""
    Container type for the state of one spin type for cDIIS
"""
mutable struct cDIISstate
    fock::CircularBuffer
    error::CircularBuffer
    errorOverlaps::CircularBuffer
end

"""
    cDIIS
"""
mutable struct cDIIS <: Accelerator
    n_diis_size::Int
    state::Tuple{cDIISstate, cDIISstate}
    conditioning_threshold::Float64
    coefficient_threshold::Float64

    function cDIIS(problem::ScfProblem; n_diis_size = 5, conditioning_threshold = 1e-14, coefficient_threshold = 1e-6, kwargs...)
        stateα = cDIISstate(CircularBuffer{AbstractArray}(n_diis_size),
            CircularBuffer{AbstractArray}(n_diis_size),
            CircularBuffer{AbstractArray}(sum(1:n_diis_size -1))
        )
        stateβ = deepcopy(stateα)
        new(n_diis_size, (stateα, stateβ), conditioning_threshold, coefficient_threshold)
    end
end

"""
    Helper function.
    pushes current fock and error matrices to states of both spin types
"""
function push_iterate_to_state!(iterate::ScfIterState, states::Tuple)
    for i in 1:size(iterate.fock, 3)
        pushfirst!(states[i].fock,  view(iterate.fock, :,:,i))
        pushfirst!(states[i].error, view(iterate.error_pulay, :,:,i))
    end
end

"""
    Computes next iterate using cDIIS
"""
function compute_next_iterate(acc::cDIIS, iterate::ScfIterState)
    # Check if the number of known fock and error matrices is equal for both
    # spins before doing anything
    history_size = acc.state[1].fock.length
    for i in 1:size(iterate.fock, 3)
        @assert acc.state[i].fock.length == history_size
        @assert acc.state[i].error.length == history_size
    end

    # Store the current fock and error matrices first
    push_iterate_to_state!(iterate, acc.state)

    # Calculate the new fock matrix for each spin type separately
    # and write them in the result matrix directly to save memory
    new_iterate_fock = zeros(size(iterate.fock))
    for i in 1:size(iterate.fock, 3)
        compute_next_iterate_fock!(acc.state[i], acc.n_diis_size, view(new_iterate_fock, :,:,i), acc.conditioning_threshold, acc.coefficient_threshold)
    end

    return update_iterate_matrix(iterate, new_iterate_fock)
end

"""
    Computes a new matrix to be used as Fock Matrix in next iteration
    The function writes the result into the passed argument 'fock'
"""
function compute_next_iterate_fock!(state::cDIISstate, n_diis_size::Int, fock::SubArray, conditioning_threshold::Float64, coefficient_threshold::Float64)
    @assert n_diis_size > 0
    @assert state.fock.length > 0

    # The Fock Matrix in the next Iteration is a linear combination of
    # previous Fock Matrices.
    #
    # The linear system has dimension m <= n_diis_size, since in the
    # beginning of the iteration we do not have the full number of
    # previous Fock Matrices yet.

    m = min(n_diis_size, state.fock.length)

    # Build the linear system we need to solve in order to obtain
    # the neccessary coefficients c_i.

    A = diis_build_matrix(n_diis_size, m, state)
    rhs = diis_build_rhs(n_diis_size, m, state)
    (c, bad_condition) = diis_solve_coefficients(A, rhs, conditioning_threshold)

    # In some cases the matrix so badly conditioned, that we cannot produce a
    # sane new iterate. In this case we need to reuse the newest fock matrix
    if bad_condition
        fock .= state.fock[1]
        return
    end

    # add very small coefficients to the largest one but always use the most
    # recent fock matrix regardless of the coefficient value
    mask = map(x -> norm(x) > coefficient_threshold, c)
    mask[1] = true
    c[argmax(c)] += sum(c[ .! mask])

    # Construct new Fock Matrix using obtained coefficients
    # and write it to the given fock matrix. We assume, that
    # fock is a matrix of zeros.
    for i in eachindex(c)[mask]
        fock .+= c[i] * state.fock[i]
    end
end

"""
    Solves the linear system after removing small eigenvalues to improve
    consistency.

    Returns the vector c and and a boolean value representing if the matrix A
    is so badly conditioned, that the previous fock matrix should be used.
"""
function diis_solve_coefficients(A, rhs, threshold)
    println(cond(A))
    # calculate the eigenvalues of A and select sufficiently large eigenvalues
    λ, U = eigen(A)
    mask = map(x -> norm(x) > threshold, λ)

    # if all eigenvalues are under the threshold, return and indicate, that the
    # previous fock matrix should be used without modifications.
    if all( .! mask)
        println("All eigenvalues are under the threshold! Skipping iterate modification…")
        return (nothing, true)
    end

    # Obtain the solution of the linear system A * c = rhs using the inverse
    # matrix of A constructed using the above decomposition
    c = U[:,mask] * Diagonal(1 ./ λ[mask]) * U[:,mask]' * rhs

    # Warning: Note that c has size (n_diis_size + 1) since
    #          the last element is the lagrange multiplier
    #          corresponding to the constraint
    #          \sum_{i=1}^n_diis_size c_i = 1
    #          We need to remove this element!

    return (c[1:length(c) - 1], false)
end

"""
    Linear System Matrix for the cDIIM accelerator.
    This is a hermitian matrix containing error overlaps B
    and ones in the form
    
    B  1
    1† 0
"""
function diis_build_matrix(n_diis_size::Int, m::Int, state::cDIISstate)
    A = zeros(m +1,m +1)

    # Since the Matrix A is Hermitian, we only have to calculate
    # the upper triagonal and can use the Julia function 'Hermitian'
    # to fill the lower triagonal of the matrix. This also allowes
    # Julia to use optimized algorithems for Hermitian matrices.

    # We can reuse most of the matrix B from the last iteration,
    # since the values do not change and only have to calculate
    # the first lign of B.
    # We accomplish this using a circular buffer of size n_diis_size,
    # which holds circular buffers of size n_diis_size to store already
    # calculated elements.

    # Fill the first row with newly calculated values and cache them
    # in a newly created Circular Buffer
    newValues = CircularBuffer{Any}(n_diis_size)
    for j in 1:m
        A[1,j] = tr(state.error[1]' * state.error[j])
        push!(newValues, A[1,j])
    end
    # Since we want to use this buffer as the 2nd row of A in the next
    # iteration we need the following layout of the buffer
    #   0 a1 a2 … am
    # so we need to push a 0 at the beginning
    pushfirst!(newValues, 0)

    # The last element of each row has to be 1, see above for a
    # detailed explanation.
    A[1, m + 1] = 1

    # Now fill the rest of the lines using cached values,
    # push a '0' on each buffer to prepare it for the next iteration
    # and set the last element of each row to 1.
    for i in 2:m
        A[i,1:m] = state.errorOverlaps[i-1][1:m]
        pushfirst!(state.errorOverlaps[i-1], 0)
        A[i, m + 1] = 1
    end

    # Push newly calculated row to the row buffer.
    pushfirst!(state.errorOverlaps, newValues)

    return Symmetric(A)
end

"""
    Right hand side of the cDIIS linear system.
    This is a vector of size m+1 containing only ones
    except for the last element which is zero.
"""
function diis_build_rhs(n_diis_size::Int, m::Int, state::cDIISstate)
    rhs = zeros(m + 1)
    rhs[m + 1] = 1
    return rhs
end
