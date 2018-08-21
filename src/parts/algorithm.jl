# Algorithms need to be initialized
Base.iterate(::T) where {T <: Algorithm} = throw("You need to initialize the algorithm first")

function (new_algorithm_type::Type{T})(old_algorithm::Algorithm, options...) where {T<:Algorithm}
    ChainedAlgorithm(old_algorithm, new_algorithm_type(options...))
end

function iterate(::Algorithm, ::SubReport)
    error("Please overload function appy")
end

function initialize_if_neccessary(alg::Algorithm, problem::ScfProblem, state::ScfIterState, softdefaults::Defaults)
    applicable(initialize, alg, problem, state, softdefaults) ? initialize(alg, problem, state, softdefaults) : state
end

# Iteration functions for algorithms
Base.iterate(report::Report) = report.state
function Base.iterate(report::Report, ::ScfIterState)
    subreport = SubReport(missing, report.problem, missing, Vector{ReportMessage}(), report.history[end], report.loglevel, report)
    res = iterate(report.algorithm, subreport)
    
    # If the underlying algorithm is done, we are done :-)
    if res == nothing
        println("tatitata")
        log!(subreport, "Algorithm is done", :debug)
        return nothing
    end

    # calculate convergence information
    report.convergence = check_convergence(report.iterate, subreport.iterate)

    # store results
    report.state = subreport.state
    push!(report.history, res[2])

    # log results
    log!(subreport, @sprintf(" %4d %14.8f %14.8f %14.8f %16.9g %12d\n", length(history), report.convergence.energies["1e"], report.convergence.energies["2e"], "etot", "scf_error", "n_applies"), :info)

    return report, res[2]
end

