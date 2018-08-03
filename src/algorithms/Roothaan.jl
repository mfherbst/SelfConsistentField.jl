"""
Perform one step of the Roothaan algorithm if the Fock matrix is iterated
"""
function roothan_step(problem::ScfProblem, iterate::FockIterState; kwargs...)
	orben, orbcoeff = compute_orbitals(problem, iterate.fock; kwargs...)
	density = compute_density(problem, orbcoeff)
	fock, error, energies = compute_fock(problem, density; kwargs...)
	return FockIterState(fock, error, energies, orbcoeff, orben)
end

