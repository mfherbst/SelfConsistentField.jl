"""
Compute the new density from the new orbital coefficients
"""
function compute_density(problem::ScfProblem, orbcoeff::Array)
	# TODO Later we can use a more clever object for this
	n_spins = get_n_spin(problem)
	assert(n_spins == 1)
	C = orbcoeff[:, 1:problem.n_occ]
	return C * C'
end
