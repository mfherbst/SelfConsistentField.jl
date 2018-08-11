"""
Obtain a Loewdin guess for the provided problem
and the atoms described by atom numbers and coords.

The coords array should be an array of size (n_atoms, 3)
"""
function compute_guess_loewdin(problem::ScfProblem,
                               coords::Array{Float64, 2},
                               atom_numbers::Array{Float64, 1}; kwargs...)
    orben, orbcoeff = eig(Hermitian(problem.overlap))
    orben = 1 / sqrt(orben)

    # Build Hcore guess
    orben, orbcoeff = compute_orbitals(problem, problem.Hcore; kwargs...)
    return compute_density(problem, orbcoeff)
end


