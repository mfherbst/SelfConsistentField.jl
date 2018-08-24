#
# Direct inversion in the iterative subspace
#

"""
    Container type for the state of one spin type for cDIIS
"""
mutable struct DiisState
    iterate::CircularBuffer
    error::CircularBuffer

    # errorOverlaps is a Circular Buffer containing already calculated rows of
    # the next iterate matrix. Each row is also stored as a Circular Buffer.
    errorOverlaps::CircularBuffer
    n_diis_size::Int

    function DiisState(n_diis_size::Int)
        new(CircularBuffer{AbstractArray}(n_diis_size),
            CircularBuffer{AbstractArray}(n_diis_size),
            CircularBuffer{AbstractArray}(n_diis_size),
            n_diis_size
        )
    end
end


"""
    cDIIS
"""
mutable struct cDIIS <: Accelerator
    state::Tuple{DiisState, DiisState}
    sync_spins::Bool
    conditioning_threshold::Float64
    coefficient_threshold::Float64

    function cDIIS(problem::ScfProblem; n_diis_size = 5, sync_spins = true, conditioning_threshold = 1e-14, coefficient_threshold = 1e-6, kwargs...)
        stateα = DiisState(n_diis_size)
        stateβ = DiisState(n_diis_size)
        new((stateα, stateβ), sync_spins, conditioning_threshold, coefficient_threshold)
    end
end

"""
    Helper function.
    pushes current iterate and error matrices to states of both spin types
"""
function push_iterate!(state::DiisState, iterate::AbstractArray, error::Union{AbstractArray,Nothing} = nothing)
    pushfirst!(state.iterate,  iterate)

    # Push difference to previous iterate if no error given
    pushfirst!(state.error,
               error != nothing ? error : iterate - state.iterate[1])
end

"""
    Helper functions to get views for specific spins
"""
function spin(obj::AbstractArray, dim::Int)
    view(obj, ntuple(x -> Colon(), ndims(obj) - 1)..., dim)
end

function spincount(obj::AbstractArray)
    size(obj, ndims(obj))
end

function spinloop(obj::Union{AbstractArray, Accelerator})
    typeof(obj) == Accelerator ?
        (1:spincount(obj.state.iterate)) :
        (1:spincount(obj))
end

"""
    Computes next iterate using cDIIS
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
    matrix(σ) = diis_build_matrix(acc.state[σ])
    coefficients(A) = diis_solve_coefficients(A, acc.conditioning_threshold)
    compute(c, σ) = compute_next_iterate_matrix!(acc.state[σ], c, spin(new_iterate, σ), acc.coefficient_threshold)

    # If sync_spins is enabled, we need to calculate the coefficients using the
    # merged matrix. This also means we need to remove the same number of
    # matrices from both histrories.
    if acc.sync_spins & (spincount(iterstate.fock) == 2)
        A = merge_matrices(matrix(1), matrix(2))
        c, matrixpurgecount = coefficients(A)
        map(σ -> compute(c, σ), spinloop(iterstate.fock))
        map(σ -> purge_matrix_from_state(acc.state[σ], matrixpurgecount), spinloop(iterstate.fock))
    else
        # If we are calculating the spins separately, each spin has its own coefficients.
        for σ in spinloop(iterstate.fock)
            c, matrixpurgecount = coefficients(matrix(σ))
            compute(c, σ)
            purge_matrix_from_state(acc.state[σ], matrixpurgecount)
        end
    end
    return update_iterate_matrix(iterstate, new_iterate)
end

"""
    When synchronizing spins the resulting DIIS matrices have to be added
    together but the constraint must be kept as is.
"""
function merge_matrices(A1::AbstractArray, A2::AbstractArray)
    view(A1, 1:size(A1, 1) - 1, 1:size(A1, 2) - 1) .+ view(A2, 1:size(A2, 1) - 1, 1:size(A2, 2) - 1)
    return A1
end

function purge_matrix_from_state(state::DiisState, count::Int)
    for i in 1:2*count
        pop!(state.iterate)
        pop!(state.error)
        pop!(state.errorOverlaps)
    end
end

"""
    Computes a new matrix to be used as Fock Matrix in next iteration
    The function writes the result into the passed argument 'fock'
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

"""
    Solves the linear system after removing small eigenvalues to improve
    consistency.

    Returns the vector c and and a boolean value representing if the matrix A
    is so badly conditioned, that the previous fock matrix should be used.
"""
function diis_solve_coefficients(A::AbstractArray, threshold::Float64)
    # Right hand side of the equation
    rhs = diis_build_rhs(size(A, 1))

    # calculate the eigenvalues of A and select sufficiently large eigenvalues
    λ, U = eigen(A)
    mask = map(x -> norm(x) > threshold, λ)

    if !all(mask)
        println("   Removing ", count(.! mask), " of ", length(mask), " eigenvalues from DIIS linear system.")
    end

    # if all eigenvalues are under the threshold, we cannot calculate sane
    # coefficients. The current fock matrix should be used without
    # modifications.
    if all( .! mask)
        println("All eigenvalues are under the threshold! Skipping iterate modification…")
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
    Linear System Matrix for the cDIIM accelerator.
    This is a hermitian matrix containing error overlaps B
    and ones in the form
    
    A = B  1
        1† 0
"""
function diis_build_matrix(state::DiisState)
    @assert state.n_diis_size > 0
    @assert state.iterate.length > 0

    # The Fock Matrix in the next Iteration is a linear combination of
    # previous Fock Matrices.
    #
    # The linear system has dimension m <= state.n_diis_size, since in the
    # beginning of the iteration we do not have the full number of
    # previous Fock Matrices yet.

    m = min(state.n_diis_size, length(state.iterate))

    A = zeros(m +1,m +1)

    # Since the Matrix A is symmetric, we only have to calculate
    # the upper triagonal and can use the Julia object 'Symmetric'
    # to fill the lower triagonal of the matrix. This also allows
    # Julia to use optimized algorithems for symmetric matrices.
    #
    # The matrix A is filled in the following form:
    #
    # B[1,1] B[1,2] B[1,3] … B[1,m] 1
    # 0      B[2,2] B[2,3] … B[2,m] 1
    # 0      0      B[3,3] … B[3,m] 1
    # ⋮      ⋮      ⋮      … ⋮      ⋮
    # 0      0      0      … B[m,m] 1
    # 0      0      0      … 0      0
    #
    # Note, that the values A[1,1] … A[1,m-1] will become the values
    # of A[2,2] … A[2,m] in the next iteration. This means, that we effectively
    # only have to calculate the first row and can reuse all consecutive ones.
    #
    # We would like to store the rows in such a way, that the storage variable
    # can be mapped to a row immediately. Since the values are shifted to the
    # right in each iteration and the last element is discarded it makes sense
    # to use a Circular Buffer to hold each row.
    #
    # Since the last row is also discarded after the new row is calculated we
    # use a Circular Buffer again to store the rows.

    # Fill the first row with newly calculated values and cache them
    # in a newly created Circular Buffer
    newValues = CircularBuffer{Any}(state.n_diis_size)
    map(j -> push!(newValues,
                   tr(state.error[1]' * state.error[j])), 1:m)
    fill!(newValues, 0)

    # Push newly calculated row to the row buffer.
    pushfirst!(state.errorOverlaps, newValues)

    # The last element of each row of A has to be 1. After calling Symmetric(A)
    # the copy of these 1s in the bottom row of A defines the constraint
    #   sum(c) == 1.
    A[1, m + 1] = 1

    # Now fill all rows with cached values,
    # push a '0' on each buffer to prepare it for the next iteration
    # and set the last element of each row to 1.
    for i in 1:m
        A[i,1:m] = state.errorOverlaps[i][1:m]
        A[i, m + 1] = 1

        # Since we want to use this buffer as the 2nd row of A in the next
        # iteration we need the following layout of the buffer
        #   0 A[1,1] A[1,2] … A[1,m-1]
        # so we need to push a 0 to the beginning
        pushfirst!(state.errorOverlaps[i], 0)
    end

    return Symmetric(A)
end

"""
    Right hand side of the cDIIS linear system.
    This is a vector of size m+1 containing only ones
    except for the last element which is zero.
"""
function diis_build_rhs(vectorsize::Int)
    rhs = zeros(vectorsize)
    rhs[end] = 1
    return rhs
end
