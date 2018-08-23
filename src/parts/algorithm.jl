# Algorithms need to be initialized
Base.iterate(::T) where {T <: Algorithm} = throw("You need to initialize the algorithm first")

"""
    In case we have an algorithm with an execution condition in Tuple form it
    needs to be converted into an algorithm before application
"""
function (new_algorithm_type::Type{T})(args...) where {T<:Algorithm}
    conv_to_alg(i) = args[i] isa Tuple{Algorithm, Function} ?
                        ConditionalExec(args[i][1], args[i][2]) :
                        args[i]

    new_algorithm_type(ntuple(conv_to_alg, length(args))...)
end

#function (new_algorithm_type::Type{T})(old_algorithm::Algorithm, options...) where {T<:Algorithm}
#    ScfPipeline(old_algorithm, new_algorithm_type(options...))
#end

function iterate(::Algorithm, ::SubReport)
    error("Please overload the function iterate")
end

function initialize_if_neccessary(alg::Algorithm, problem::ScfProblem, state::ScfIterState, params::Parameters)
    applicable(initialize, alg, problem, state, params) ? initialize(alg, problem, state, params) : state
end

# Iteration functions for algorithms
Base.iterate(report::Report) = (report, report.history[end])
function Base.iterate(report::Report, lastsubreport::SubReport)
    iterresult = iterate(report.algorithm, lastsubreport)
    
    # If the underlying algorithm is done, we are done :-)
    if iterresult == nothing
        log!(report.history[end], "Algorithm is done", :info)
        return nothing
    end

    # store results
    (algorithm, subreport) = iterresult
    push!(report.history, subreport)

    # log results
    log!(subreport, @sprintf(" %4d %14.8f %14.8f %14.8f %16.9g %12d", length(report.history) - 1, report.state.energies["1e"], report.state.energies["2e"], report.state.energies["total"], !ismissing(report.convergence) ? report.convergence.error_norm : NaN, NaN), :info)

    # update state reference in report
    report.state = subreport.state
    report.convergence = subreport.convergence
    return (report, subreport)
end

