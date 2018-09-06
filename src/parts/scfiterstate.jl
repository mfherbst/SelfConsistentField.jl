Base.IteratorSize(::ScfIterState) = Base.SizeUnknown()

# Iteration functions for algorithms
Base.iterate(report::ScfIterState) = (report, report.history[end])
function Base.iterate(report::ScfIterState, laststepstate::StepState)

    laststepstate.is_converged && return nothing

    iterresult = iterate(report.algorithm, laststepstate)
    
    # If the underlying algorithm is done, we are done :-)
    if iterresult == nothing

        log!(report.history[end], "Algorithm is done but convergence was not reached.", :info)
        return nothing
    end

    # store results
    (algorithm, stepstate) = iterresult
    push!(report.history, stepstate)

    # update iterate reference in report
    report.algorithm = stepstate.algorithm
    report.iterate = stepstate.iterate
    report.error_norm = stepstate.error_norm
    report.energy_change = stepstate.energy_change
    report.is_converged = stepstate.is_converged

    # log results
    log!(stepstate, @sprintf(" %4d %14.8f %14.8f %14.8f %16.9g %12d", length(report.history) - 1, report.iterate.energies["1e"], report.iterate.energies["2e"], report.iterate.energies["total"], report.error_norm, NaN), :info)

    return (report, stepstate)
end

function Base.convert(::Type{T}, rp::ScfIterState) where {T <: AbstractDict}
    Dict(
         "n_iter"=>length(rp.history) - 1,
         "n_applies"=>NaN,
         "problem"=>rp.problem,
         "orben"=>rp.iterate.orben,
         "orbcoeff"=>rp.iterate.orbcoeff,
         "density"=>compute_density(rp.problem,rp.iterate.orbcoeff),
         "converged"=> rp.is_converged,
         "fock"=>rp.iterate.fock,
         "energies"=>rp.iterate.energies,
         "error_norm"=>rp.error_norm,
         "energy_change"=>rp.energy_change
    )
end
