"""
Obtain an Hcore guess for the provided problem
and the atoms described by atom numbers and coords.

The coords array should be an array of size (n_atoms, 3)

The returned density is only the spatial part, i.e.
of size (n_bas, n_bas) and makes no difference between alpha and beta components.
"""
function compute_guess_hcore(problem::ScfProblem,
		             coords::Array{Float64, 2},
			     atom_numbers::Array{Float64, 1}; kwargs...)
	n_bas, _ = size(problem.Hcore)
	assert(size(problem.Hcore) == (n_bas, n_bas))
	Fcore = reshape(problem.Hcore, n_bas, n_bas, 1)

	# Compute the orbitals by diagonalising Hcore
	# Compute the resulting density and return
	_, orbcoeff = compute_orbitals(problem, Fcore; kwargs...)
	return compute_density(problem, orbcoeff)
end

