"""
    EDIIS
"""
mutable struct EDIIS <: Algorithm
    n_diis_size::Union{Missing, Int}
    coefficient_threshold::Union{Missing,Float64}
    state::Tuple{DiisState, DiisState}

    function EDIIS(;n_diis_size::Union{Missing, Int} = missing, coefficient_threshold::Union{Missing, Float64} = missing)
        new(n_diis_size, coefficient_threshold)
    end
end

function initialize(ediis::EDIIS, problem::ScfProblem, iterstate::ScfIterState, defaults::Defaults)
    # TODO needs to become a separate function using reflection
    ediis.n_diis_size = ismissing(ediis.n_diis_size) & :n_diis_size ∈ defaults ? defaults[:n_diis_size] : 5
    ediis.coefficient_threshold = ismissing(ediis.coefficient_threshold) & :coefficient_threshold ∈ defaults ? defaults[:coefficient_threshold] : 10e-6

    stateα = DiisState(n_diis_size)
    ediis.state = spincount(get_iterate_matrix(iterstate)) == 2 ? (stateα, DiisState(n_diis_size)) : (stateα)
end

function iterate(ediis::EDIIS, rp::SubReport)
    rp.state = compute_next_iterate(ediis, rp.source.state)
    rp
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
    b(i,j) = tr((state.iterate[i] - state.iterate[j]) * (state.density[i] - state.density[j]))

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
    B1 + B2
end
