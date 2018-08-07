using Printf

# TODO Find a way to get the print statements out
function run_scf(problem::ScfProblem, guess_density::AbstractArray;
                 max_iter=100, damping_max_error_norm=1e-2,
                 damping_max_energy_total_change=0.0025,
                 kwargs...)
    # Setup accelerators and SCF-global objects
    damping = FixedDamping(problem; kwargs...)
    scfconv = nothing

    # Build initial iterate
    fock, error, energies = compute_fock(problem, guess_density; kwargs...)
    iterate = FockIterState(fock, error, energies, nothing, nothing)

    @printf("%5s %14s %14s %14s %15s %12s\n",
            "iter", "e1e", "e2e", "etot", "scf_error", "n_applies")
    for i in 1:max_iter
        # Extrapolate next iterate
        iterate = compute_next_iterate(damping, iterate)

        # and perform a step to progress
        newiterate = roothan_step(problem, iterate; kwargs...)
        n_applies = NaN # Number of applies of the fock matrix, which was required

        # Compute convergence state
        scfconv = check_convergence(iterate, newiterate; kwargs...)

        # Print current status
        energies = newiterate.energies
        @printf(" %4d %14.8f %14.8f %14.8f %15.9g %12d\n",
                i, energies["energy_1e"], energies["energy_2e"],
                energies["energy_total"], scfconv.error_norm,
                n_applies)

        # If converged end iteration
        if scfconv.is_converged break end

        remove_damping = (scfconv.error_norm < damping_max_error_norm
                          || abs(scfconv.energy_change["energy_total"])
                          < damping_max_energy_total_change)
        if remove_damping && damping != nothing
            println(repeat(" ", 18), "**** Removing any damping ****")
            damping = nothing
        end
        iterate = newiterate
    end

    # Return results
    return Dict(
        "orben"=>iterate.orben,
        "orbcoeff"=>iterate.orbcoeff,
        "density"=>compute_density(problem,iterate.orbcoeff),
        "is_converged"=> scfconv.is_converged,
        "fock"=>iterate.fock,
        "energies"=>iterate.energies,
        "error_norm"=>scfconv.error_norm,
        "energy_change"=>scfconv.energy_change,
    )
end
