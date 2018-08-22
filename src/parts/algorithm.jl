# Algorithms need to be initialized
Base.iterate(::T) where {T <: Algorithm} = throw("You need to initialize the algorithm first")

function (new_algorithm_type::Type{T})(old_algorithm::Algorithm, options...) where {T<:Algorithm}
    ChainedAlgorithm(old_algorithm, new_algorithm_type(options...))
end

function iterate(::Algorithm, ::SubReport)
    error("Please overload the function iterate")
end

function initialize_if_neccessary(alg::Algorithm, problem::ScfProblem, state::ScfIterState, softdefaults::Defaults)
    applicable(initialize, alg, problem, state, softdefaults) ? initialize(alg, problem, state, softdefaults) : state
end

# Iteration functions for algorithms
Base.iterate(report::Report) = (report, report.history[end])
function Base.iterate(report::Report, lastsubreport::SubReport)
    iterresult = iterate(report.algorithm, lastsubreport)
    
    # If the underlying algorithm is done, we are done :-)
    if iterresult == nothing
        log!(report.history[end], "Algorithm is done", :debug)
        return nothing
    end

    # store results
    (algorithm, subreport) = iterresult
    push!(report.history, subreport)

    # calculate convergence information
    report.convergence = check_convergence(report.state, subreport.state)

    # log results
    log!(subreport, @sprintf(" %4d %14.8f %14.8f %14.8f %16.9g %12d\n", length(report.history), report.state.energies["1e"], report.state.energies["2e"], report.state.energies["total"], report.convergence.error_norm, NaN), :info)

    # update state reference in report
    report.state = subreport.state
    return (report, subreport)
end

