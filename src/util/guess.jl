"""
Obtain an Hcore guess for the provided problem
and the atoms described by atom numbers and coords.

The coords array should be an array of size (n_atoms, 3)

The returned density is of the form expected for the ScfProblem,
i.e. only has a single spin component for restricted-closed SCF,
but two for unrestricted SCF or restricted-open-shell SCF.
"""
function compute_guess_hcore(problem::ScfProblem,
                             coords::Array{Float64, 2},
                             atom_numbers::Array{T, 1}; kwargs...) where T <: Real
    n_bas, _ = size(problem.h_core)
    @assert size(problem.h_core) == (n_bas, n_bas) "Hcore should be a square matrix."
    Fcore = reshape(problem.h_core, n_bas, n_bas, 1)

    # Compute the orbitals by diagonalising Hcore
    # Compute the resulting density and return
    _, orbcoeff = compute_orbitals(problem, Fcore; kwargs...)
    return compute_density(problem, orbcoeff)
end


"""
Obtain a Loewdin guess for the provided problem
and the atoms described by atom numbers and coords.

The coords array should be an array of size (n_atoms, 3)
"""
function compute_guess_loewdin(problem::ScfProblem,
                               coords::Array{Float64, 2},
                               atom_numbers::Array{T, 1}; kwargs...) where T <: Real
    # TODO Actually we only need max(problem.n_occ) eigenvectors here
    λ, orbcoeff = eigen(Symmetric(problem.overlap))

    # Orthonormalise coefficients
    orbcoeff = orbcoeff * Diagonal(1 ./ sqrt.(λ))

    # Permute coefficients such that largest eigenvalues of overlap
    # are first
    orbcoeff = reverse(orbcoeff, dims=2)

    orbcoeff = cat(orbcoeff, dims=3)
    return compute_density(problem, orbcoeff)
end


"""
Obtain a totally random guess density for the provided problem
and the atoms described by atom numbers and coords.

The coords array should be an array of size (n_atoms, 3)

The returned density is of the form expected for the ScfProblem,
i.e. only has a single spin component for restricted-closed SCF,
but two for unrestricted SCF or restricted-open-shell SCF.
"""
function compute_guess_random(problem::ScfProblem,
                              coords::Array{Float64, 2},
                              atom_numbers::Array{T, 1}; kwargs...) where T <: Real
    # TODO Actually we only need max(problem.n_occ) vectors here
    S = problem.overlap

    λ, orbcoeff = eigen(Symmetric(randn(size(S)) + S), Symmetric(S))

    # Permute coefficients such that largest eigenvalues of overlap
    # are first
    orbcoeff = reverse(orbcoeff, dims=2)

    orbcoeff = cat(orbcoeff, dims=3)
    return compute_density(problem, orbcoeff)
end


"""
Obtain a guess density for the provided problem
and the atoms described by atom numbers and coords.

The coords array should be an array of size (n_atoms, 3)

The returned density is of the form expected for the ScfProblem,
i.e. only has a single spin component for restricted-closed SCF,
but two for unrestricted SCF or restricted-open-shell SCF.

Valid guess methods are
    - loewdin
    - hcore (default)
    - random
"""
function compute_guess(problem::ScfProblem,
                       coords::Array{Float64, 2},
                       atom_numbers::Array{T, 1}; method="hcore", kwargs...) where T <: Real
    avail_methods = Dict(
         "hcore" => compute_guess_hcore,
         "loewdin" => compute_guess_loewdin,
         "random" => compute_guess_random,
    )

    if ! (method in keys(avail_methods))
        error("Method $method is not an available guess method, " *
              "try one of $( join(keys(avail_methods), " ") )")
    end
    return avail_methods[method](problem, coords, atom_numbers; kwargs...)
end
