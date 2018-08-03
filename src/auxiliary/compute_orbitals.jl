"""
Diagonalise the Fock matrix and return a new eigensolution
"""
function compute_orbitals(problem::ScfProblem, fock::AbstractArray; kwargs...)
	# Obtain new mo energies and coefficients
	# by diagonalising the fock matrix and overlap
	# matrix on input.
	n_orbs = problem.n_orbs
	orben, orbcoeff = eig(Hermitian(fock), Hermitian(problem.overlap))
	return orben[1:n_orbs], orbcoeff[:, 1:n_orbs]
end
