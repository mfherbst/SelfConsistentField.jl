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

"""
    Calculates coefficients
"""
function diis_solve_coefficients(::EDIIS, B::AbstractArray; energybuffer::AbstractArray, kwargs...)
    m = size(B,1)
    E = map(energies -> energies["total"], energybuffer)

    function f(x)
        c = x.^2/sum(x.^2)
        E'*c - 1/2 * c'*B*c
    end

    function gradf(x)
        c = x.^2/sum(x.^2)
        ((Diagonal(x * sum(x.^2)) - x.^2*x').*(2/(sum(x.^2)^2)))*(E - B*c)
        #Diagonal((x * sum(x.^2) - x.^3).*(2/(sum(x.^2)^2)))*(E - 1/2 *(B*c+diag(B).*c))
    end

    res = optimize(f, gradf, ones(m), BFGS(); inplace = false)
    c = Optim.minimizer(res)

    return c, 0
end

"""
    Linear System Matrix for the EDIIS accelerator.
"""
function diis_build_matrix(::EDIIS, state::DiisState)
    @assert state.n_diis_size > 0
    @assert state.iterate.length > 0

    # B has dimension m <= state.n_diis_size, since in the
    # beginning of the iteration we do not have the full number of
    # previous iterates yet.

    m = min(state.n_diis_size, length(state.iterate))

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

    # Definition of matrix elements
    b(i,j) = tr(state.iterate[i] * state.density[i]) - tr(state.iterate[i] * state.density[j]) - tr(state.iterate[j] * state.iterate[i]) + tr(state.iterate[j] * state.iterate[j])

    # Fill the first row with newly calculated values and cache them
    # in a newly created Circular Buffer
    newValues = CircularBuffer{Any}(state.n_diis_size)
    map(j -> push!(newValues, b(1,j)), 1:m)
    fill!(newValues, 0)

    # Push newly calculated row to the row buffer. We use the iterationstate to
    # store these.
    pushfirst!(state.iterationstate, newValues)

    # Now fill all rows with cached values,
    # push a '0' on each buffer to prepare it for the next iteration
    for i in 1:m
        B[i,1:m] = state.iterationstate[i][1:m]

        # Since we want to use this buffer as the 2nd row of A in the next
        # iteration we need the following layout of the buffer
        #   0 A[1,1] A[1,2] … A[1,m-1]
        # so we need to push a 0 to the beginning
        pushfirst!(state.iterationstate[i], 0)
    end

    return Symmetric(B)
end

function needs_error(::EDIIS)
    false
end

function needs_density(::EDIIS)
    true
end

"""
    When synchronizing spins the resulting DIIS matrices have to be added
    together but the constraint must be kept as is.
"""
function merge_matrices(::EDIIS, B1::AbstractArray, B2::AbstractArray)
    B1 .+ B2
end
