"""
    DIIS Matrix.
"""
function compute_diis_matrix(diis_matrix_formula::Function, history::DiisHistory)
    @assert history.n_diis_size > 0
    @assert length(history.iterate) > 0

    # B has dimension m <= history.n_diis_size, since in the
    # beginning of the iteration we do not have the full number of
    # previous iterates yet.

    m = min(history.n_diis_size, length(history.iterate))

    B = zeros(m,m)

    # Since the Matrix B is symmetric, we only have to calculate
    # the upper triagonal and can use the Julia object 'Symmetric'
    # to fill the lower triagonal of the matrix. This also allows
    # Julia to use optimized algorithems for symmetric matrices.
    #
    # The matrix B is filled in the following form:
    #
    # B[1,1] B[1,2] B[1,3] … B[1,m]
    # 0      B[2,2] B[2,3] … B[2,m]
    # 0      0      B[3,3] … B[3,m]
    # ⋮      ⋮      ⋮      … ⋮     
    # 0      0      0      … B[m,m]
    #
    # Note, that the values B[1,1] … B[1,m-1] will become the values
    # of B[2,2] … B[2,m] in the next iteration. This means, that we effectively
    # only have to calculate the first row and can reuse all consecutive ones.
    #
    # We would like to store the rows in such a way, that the storage variable
    # can be mapped to a row immediately. Since the values are shifted to the
    # right in each iteration and the last element is discarded it makes sense
    # to use a Circular Buffer to hold each row.
    #
    # Since the last row is also discarded after the new row is calculated we
    # use a Circular Buffer again to store the rows.

    # Rename diis_matrix_formula for better readability.
    # The entry in the ith row and jth column of B is b(i,j)
    b = diis_matrix_formula

    # Fill the first row with newly calculated values and cache them
    # in a newly created Circular Buffer
    newValues = CircularBuffer{Any}(history.n_diis_size)
    map(j -> push!(newValues, b(1,j)), 1:m)
    fill!(newValues, 0)

    # Push newly calculated row to the row buffer. We use the iterationhistory to
    # store these.
    new_history = copy(history)
    pushfirst!(new_history.iterationhistory, newValues)

    # Now fill all rows with cached values,
    # push a '0' on each buffer to prepare it for the next iteration
    for i in 1:m
        B[i,1:m] = new_history.iterationhistory[i][1:m]

        # Since we want to use this buffer as the 2nd row of A in the next
        # iteration we need the following layout of the buffer
        #   0 A[1,1] A[1,2] … A[1,m-1]
        # so we need to push a 0 to the beginning
        pushfirst!(new_history.iterationhistory[i], 0)
    end

    return Symmetric(B), new_history
end

"""
    Linear System Matrix for the cDIIS accelerator.
    This is a hermitian matrix containing error overlaps B
    and ones in the form
    
    A = B  1
        1† 0
"""
function build_diis_linear_system_matrix(B::AbstractArray)
    @assert size(B, 1) == size(B, 2)
    m = size(B, 1)
    # The matrix A is filled in the following form:
    #
    # B[1,1] B[1,2] B[1,3] … B[1,m] 1
    # 0      B[2,2] B[2,3] … B[2,m] 1
    # 0      0      B[3,3] … B[3,m] 1
    # ⋮      ⋮      ⋮      … ⋮      ⋮
    # 0      0      0      … B[m,m] 1
    # 1      1      1      … 1      0
    #
    A = ones(m+1, m+1)
    A[end,end] = 0
    A[1:size(B, 1), 1:size(B, 2)] .= B
    return Symmetric(A)
end

function compute_next_iterate(iterate::AbstractArray, history::Tuple{DiisHistory,DiisHistory}, diis_matrix_formula::Function, compute_diis_coefficients::Function, lg::Logger; params...)
    # To save memory we store only new_iterate once and pass views of it to the
    # computation routines that write directly into the view.
    new_iterate = zeros(size(iterate))
    purgecount = 0

    # Only store precaclulated values for the DIIS matrix in the first spinblock history.
    diis_matrix, new_historyα = compute_diis_matrix(diis_matrix_formula, history[1])
    c, purgecount = compute_diis_coefficients(diis_matrix, lg; params...)
    new_history = (new_historyα, history[2])
    for σ in spinloop(iterate)
        apply_diis_coefficients!(new_history[σ], c, spin(new_iterate, σ); params...)
    end

    return new_iterate, new_history, (purgecount, purgecount)
end

function compute_next_iterate_separating_spins(iterate::AbstractArray, history::Tuple{DiisHistory,DiisHistory}, diis_matrix_spinblock_formulas::Tuple{Function,Function}, compute_diis_coefficients::Function, lg::Logger; params...)
    new_iterate = zeros(size(iterate))
    new_history = Vector{DiisHistory}
    purgecount = (0,0)

    # If we are calculating the spins separately, each spin has its own coefficients.
    for σ in spinloop(iterate)
        diis_matrix, new_history[σ] = compute_diis_matrix(diis_matrix_formula, history[σ])
        c, purgecount[σ] = compute_diis_coefficients(diis_matrix; params...)
        apply_diis_coefficients!(new_history[σ], c, spin(new_iterate, σ); params...)
    end
    return new_iterate, collect(new_history), purgecount
end

"""
    Computes new iterate and stores it in the passed argument "iterate"
"""
function apply_diis_coefficients!(history::DiisHistory, c::AbstractArray, iterate::SubArray; coefficient_threshold::Float64, params...)
    @assert length(c) > 0
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
