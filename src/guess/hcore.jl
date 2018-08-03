"""
Obtain an Hcore guess for the provided problem
and the atoms described by atom numbers and coords.

The coords array should be an array of size (n_atoms, 3)
"""
function compute_guess_hcore(problem::ScfProblem,
		             coords::Array{Float64, 2},
			     atom_numbers::Array{Float64, 1}; kwargs...)
	assert(size(coords)[2] == 3)
	assert(get_n_spin(problem) == 1)
	# Build Hcore guess
	orben, orbcoeff = compute_orbitals(problem, problem.Hcore; kwargs...)
	return compute_density(problem, orbcoeff)
end

