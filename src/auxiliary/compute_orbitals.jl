"""
Diagonalise the Fock matrix and return a new eigensolution
"""
function compute_orbitals(problem::ScfProblem, fock::AbstractArray; kwargs...)
    n_orb = problem.n_orb
    n_bas, _, n_spin = size(fock)
    @assert size(fock) == (n_bas, n_bas, n_spin) "fock should be a square matrix"

    # For each spin block diagonalise
    orbcoeff = Array{eltype(fock)}(undef, n_bas, n_orb, n_spin)
    orben = Array{eltype(fock)}(undef, n_orb, n_spin)
    for s in 1:n_spin
        focks = view(fock, :, :, s)
        eval, evec = eigen(Hermitian(focks), Hermitian(problem.overlap))
        orben[:, s] = eval[1:n_orb]
        orbcoeff[:, :, s] = evec[:, 1:n_orb]
    end
    return orben, orbcoeff
end
