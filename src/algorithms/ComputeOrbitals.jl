struct ComputeOrbitals <: Algorithm end

"""
Diagonalise the Fock matrix and return a new eigensolution
"""
function iterate(::ComputeOrbitals, rp::SubReport)
    lg = Logger(rp)
    log!(lg, "Computing density", :debug, :computeorbitals)

    orben, orbcoeff = compute_orbitals(rp.problem, rp.state.fock)
    log!(lg, "New orben", orben, :debug, :computeorbitals)
    log!(lg, "New orbcoeff", orbcoeff, :debug, :computeorbitals)

    state = FockIterState(rp.state.fock, rp.state.error_pulay, rp.state.energies, orbcoeff, orben, rp.state.density)
    return ComputeOrbitals(), new_subreport(ComputeOrbitals(), state, lg, rp)
end

# Separate function to allow orbitals to be computed in guesses
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
