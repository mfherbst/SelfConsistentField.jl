using TensorOperations

"""
Compute the new density from the new orbital coefficients
"""
function compute_density(problem::ScfProblem, orbcoeff::Array)
	n_bas, n_orbs, n_spin = size(orbcoeff)

	# TODO Later we can use a more clever object for this
	density = Array{eltype(orbcoeff)}(n_bas, n_bas, n_spin)

	for s in 1:n_spin
		n_occs = problem.n_occ[s]
		assert(n_occs < n_orbs)

		Cs = view(orbcoeff, :, 1:n_occs, s)
		density[:, :, s] = Cs * Cs'
	end

	return density
end
