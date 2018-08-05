"""
Diagonalise the Fock matrix and return a new eigensolution
"""
function compute_orbitals(problem::ScfProblem, fock::AbstractArray; kwargs...)
	n_orbs = problem.n_orbs
	n_bas, _, n_spin = size(fock)
	assert(size(fock) == (n_bas, n_bas, n_spin))

	# For each spin block diagonalise
	orbcoeff = Array{eltype(fock)}(n_bas, n_orbs, n_spin)
	orben = Array{eltype(fock)}(n_orbs, n_spin)
	for s in 1:n_spin
		focks = view(fock, :, :, s)
		eval, evec = eig(Hermitian(focks), Hermitian(problem.overlap))
		orben[:, s] = eval[1:n_orbs]
		orbcoeff[:, :, s] = evec[:, 1:n_orbs]
	end
	return orben, orbcoeff
end
