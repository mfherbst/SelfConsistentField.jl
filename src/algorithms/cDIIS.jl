#
# Direct inversion in the iterative subspace
#

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

function needs_error(::cDIIS)
    true
end

function needs_density(::cDIIS)
    false
end

"""
    Solves the linear system after removing small eigenvalues to improve
    consistency.
"""
function diis_solve_coefficients(cdiis::cDIIS, A::AbstractArray)
    # Right hand side of the equation
    rhs = cdiis_build_rhs(size(A, 1))

    # calculate the eigenvalues of A and select sufficiently large eigenvalues
    λ, U = eigen(A)
    mask = map(x -> norm(x) > cdiis.conditioning_threshold, λ)

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
    When synchronizing spins the resulting DIIS matrices have to be added
    together but the constraint must be kept as is.
"""
function merge_diis_matrices_spin_blocks(::cDIIS, A1::AbstractArray, A2::AbstractArray)
    view(A1, 1:size(A1, 1) - 1, 1:size(A1, 2) - 1) .+ view(A2, 1:size(A2, 1) - 1, 1:size(A2, 2) - 1)
    return A1
end

"""
    Linear System Matrix for the cDIIS accelerator.
    This is a hermitian matrix containing error overlaps B
    and ones in the form
    
    A = B  1
        1† 0
"""
function diis_build_matrix(::cDIIS, state::DiisState)
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

    # Push newly calculated row to the row buffer. We use the iterationstate to
    # store these.
    pushfirst!(state.iterationstate, newValues)

    # The last element of each row of A has to be 1. After calling Symmetric(A)
    # the copy of these 1s in the bottom row of A defines the constraint
    #   sum(c) == 1.
    A[1, m + 1] = 1

    # Now fill all rows with cached values,
    # push a '0' on each buffer to prepare it for the next iteration
    # and set the last element of each row to 1.
    for i in 1:m
        A[i,1:m] = state.iterationstate[i][1:m]
        A[i, m + 1] = 1

        # Since we want to use this buffer as the 2nd row of A in the next
        # iteration we need the following layout of the buffer
        #   0 A[1,1] A[1,2] … A[1,m-1]
        # so we need to push a 0 to the beginning
        pushfirst!(state.iterationstate[i], 0)
    end

    return Symmetric(A)
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
