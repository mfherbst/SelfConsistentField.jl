"""
Compute the new density from the new orbital coefficients
"""
function compute_density(problem::ScfProblem, orbcoeff::Array)
    n_elec = problem.n_elec
    n_bas, n_orb, n_spin = size(orbcoeff)
    @assert n_spin == 1 || n_spin == 2 "Only one or two spin components allowed"
    @assert n_elec[1] < n_orb && n_elec[2] < n_orb begin
        "n_elec of alpha and beta should be less than n_orb"
    end


    occupancies = (zeros(eltype(orbcoeff), n_orb),
                   zeros(eltype(orbcoeff), n_orb))

    # Hard-coded oxygen case
    occupancies[1][1] = 1.0; occupancies[2][1] = 1.0  # 1s
    occupancies[1][2] = 1.0; occupancies[2][2] = 1.0  # 2s
    occupancies[1][3] = 1.0; occupancies[2][3] = 1/3  # 2p
    occupancies[1][4] = 1.0; occupancies[2][4] = 1/3  # 2p
    occupancies[1][5] = 1.0; occupancies[2][5] = 1/3  # 2p

    #occupancies = ((occupancies[1] + occupancies[2]) / 2, )

    # Alpha and beta occupied coefficients
    # Notice that for the beta coefficients n_spin is used.
    # This makes sure that the same block is used for alpha and beta
    # if there only is one spin and else that the beta block is used.
    # Cocc = (view(orbcoeff, :, 1:n_elec[1], 1),
    #         view(orbcoeff, :, 1:n_elec[2], n_spin))

    # Usually the number of spin blocks in orbcoeff and the number
    # of densities to be computed matches,
    # unless only a single Fock matrix block was diagonalised,
    # but two different occupancies are present in alpha and beta
    # in this case we still need to compute two densities
    n_densities = n_spin
    if !problem.restricted || n_elec[1] != n_elec[2]
        n_densities = 2
    end

    # TODO Later we can use a more clever object for this
    density = Array{eltype(orbcoeff)}(undef, n_bas, n_bas, n_densities)
    for s in 1:n_densities
        os = min(s, n_spin)
        density[:, :, s] = (
            orbcoeff[:, :, os] * Diagonal(occupancies[s]) * orbcoeff[:, :, os]'
        )
    end
    return density
end
