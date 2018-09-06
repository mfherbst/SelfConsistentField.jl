Base.IteratorSize(::Report) = Base.SizeUnknown()

# Iteration functions for algorithms
Base.iterate(report::Report) = (report, report.history[end])
function Base.iterate(report::Report, lastsubreport::StepState)

    !ismissing(lastsubreport.convergence) && lastsubreport.convergence.is_converged && return nothing

    iterresult = iterate(report.algorithm, lastsubreport)
    
    # If the underlying algorithm is done, we are done :-)
    if iterresult == nothing

        log!(report.history[end], "Algorithm is done but convergence was not reached.", :info)
        return nothing
    end

    # store results
    (algorithm, subreport) = iterresult
    push!(report.history, subreport)

    # log results
    log!(subreport, @sprintf(" %4d %14.8f %14.8f %14.8f %16.9g %12d", length(report.history) - 1, report.state.energies["1e"], report.state.energies["2e"], report.state.energies["total"], !ismissing(report.convergence) ? report.convergence.error_norm : NaN, NaN), :info)

    # update state reference in report
    report.algorithm = subreport.algorithm
    report.state = subreport.state
    report.convergence = subreport.convergence
    return (report, subreport)
end

function Base.convert(::Type{T}, rp::Report) where {T <: AbstractDict}
    Dict(
         "n_iter"=>length(rp.history) - 1,
         "n_applies"=>NaN,
         "problem"=>rp.problem,
         "orben"=>rp.state.orben,
         "orbcoeff"=>rp.state.orbcoeff,
         "density"=>compute_density(rp.problem,rp.state.orbcoeff),
         "converged"=> rp.convergence.is_converged,
         "fock"=>rp.state.fock,
         "energies"=>rp.state.energies,
         "error_norm"=>rp.convergence.error_norm,
         "energy_change"=>rp.convergence.energy_change
    )
end
