# TODO Output as in molsturm
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
    iterate = compute_next_iterate(damping, iterate)

    println("iter      etot         echange     error norm")
    for i in 1:max_iter
        newiterate = roothan_step(problem, iterate; kwargs...)

        # Compute convergence state
        scfconv = check_convergence(iterate, newiterate; kwargs...)

        # Print current status
        print("$i   ", newiterate.energies["energy_total"], "   ")
        print(scfconv.energy_change["energy_total"], "   ")
        println(scfconv.error_norm)

        # If converged end iteration
        if scfconv.is_converged break end

        remove_damping = (scfconv.error_norm < damping_max_error_norm
                          || abs(scfconv.energy_change["energy_total"])
                          < damping_max_energy_total_change)
        if remove_damping && damping != nothing
            println("  ... removing any damping ... ")
            iterate = postprocess_iterate(damping, iterate)
            damping = nothing
        end

        # Else extrapolate next iterate
        iterate = compute_next_iterate(damping, newiterate)
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
