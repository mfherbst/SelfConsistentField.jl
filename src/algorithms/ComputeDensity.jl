"""
Compute the new density from the new orbital coefficients
"""
struct ComputeDensity <: Algorithm end

function iterate(::ComputeDensity, rp::SubReport)
    lg = Logger(rp)
    log!(lg, "Computing density", :debug, :computedensity)

    density = compute_density(rp.problem, rp.state.orbcoeff)
    log!(lg, "New density", density, :debug, :computedensity)

    state = FockIterState(rp.state.fock, rp.state.error_pulay, rp.state.energies, rp.state.orbcoeff, rp.state.orben, density)
    newalg = ComputeDensity()
    return newalg, new_subreport(newalg, state, lg, rp)
end

# Separate function to allow density to be computed in guesses
function compute_density(problem::ScfProblem, orbcoeff::Array)
    n_elec = problem.n_elec
    n_bas, n_orb, n_spin = size(orbcoeff)
    @assert n_spin == 1 || n_spin == 2 "Only one or two spin components allowed"
    @assert n_elec[1] < n_orb && n_elec[2] < n_orb begin
        "n_elec of alpha and beta should be less than n_orb"
    end

    # Alpha and beta occupied coefficients
    # Notice that for the beta coefficients n_spin is used.
    # This makes sure that the same block is used for alpha and beta
    # if there only is one spin and else that the beta block is used.
    Cocc = (view(orbcoeff, :, 1:n_elec[1], 1),
            view(orbcoeff, :, 1:n_elec[2], n_spin))

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
        density[:, :, s] = Cocc[s] * Cocc[s]'
    end
    return density
end
