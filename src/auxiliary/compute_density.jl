using TensorOperations

"""
Compute the new density from the new orbital coefficients
"""
function compute_density(problem::ScfProblem, orbcoeff::Array)
	n_occ = problem.n_occ
	n_bas, n_orbs, n_spin = size(orbcoeff)
	assert(n_spin == 1 || n_spin == 2)
	assert(n_occ[1] < n_orbs && n_occ[2] < n_orbs)

	# Alpha and beta occupied coefficients
	# Notice that for the beta coefficients n_spin is used.
	# This makes sure that the same block is used for alpha and beta
	# if there only is one spin and else that the beta block is used.
	Cocc = (view(orbcoeff, :, 1:n_occ[1], 1),
		view(orbcoeff, :, 1:n_occ[2], n_spin))

	# Usually the number of spin blocks in orbcoeff and the number
	# of densities to be computed matches,
	# unless only a single Fock matrix block was diagonalised,
	# but two different occupancies are present in alpha and beta
	# in this case we still need to compute two densities
	n_densities = n_spin
	if n_occ[1] != n_occ[2]
		n_densities = 2
	end

	# TODO Later we can use a more clever object for this
	density = Array{eltype(orbcoeff)}(n_bas, n_bas, n_densities)
	for s in 1:n_densities
		density[:, :, s] = Cocc[s] * Cocc[s]'
	end
	return density
end
