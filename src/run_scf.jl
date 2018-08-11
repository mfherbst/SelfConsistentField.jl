using Printf

# TODO Find a way to get the print statements out
function run_scf(problem::ScfProblem, guess_density::AbstractArray;
                 max_iter=100, damping_max_error_norm=0.25,
                 kwargs...)
    # Setup accelerators and SCF-global objects
    damping = FixedDamping(problem; kwargs...)
    scfconv = nothing
    switched_to_diis = false

    n_iter = 0 # Number of iterations
    n_applies = NaN # Number of applies of the fock matrix, which was required

    # Build initial iterate
    fock, error_pulay, energies = compute_fock_matrix(problem, guess_density; kwargs...)
    iterate = FockIterState(fock, error_pulay, energies, nothing, nothing)

    @printf("%5s %14s %14s %14s %15s %12s\n",
            "iter", "e1e", "e2e", "etot", "scf_error", "n_applies")
    for i in 1:max_iter
        n_iter = i

        # Extrapolate next iterate
        iterate = compute_next_iterate(damping, iterate)

        # and perform a step to progress
        newiterate = roothan_step(problem, iterate; kwargs...)

        # Compute convergence state
        scfconv = check_convergence(iterate, newiterate; kwargs...)

        # Print current status
        energies = newiterate.energies
        @printf(" %4d %14.8f %14.8f %14.8f %16.9g %12d\n",
                i, energies["1e"], energies["2e"], energies["total"], scfconv.error_norm,
                n_applies)

        # If converged end iteration
        if scfconv.is_converged break end

        if scfconv.error_norm < 0.25 && !switched_to_diis
            println("**** Switching on DIIS ****")
            damping = cDIIS(problem; kwargs...)
            switched_to_diis = true
        end

        iterate = newiterate
    end

    # Return results
    return Dict(
        "n_iter"=>n_iter,
        "n_applies"=>n_applies,
        "problem"=>problem,
        "orben"=>iterate.orben,
        "orbcoeff"=>iterate.orbcoeff,
        "density"=>compute_density(problem,iterate.orbcoeff),
        "converged"=> scfconv.is_converged,
        "fock"=>iterate.fock,
        "energies"=>iterate.energies,
        "error_norm"=>scfconv.error_norm,
        "energy_change"=>scfconv.energy_change,
    )
end
