"""
    Algorithms need to be initialized before usage. Initialization returns a
    Report which can be iterated.
"""
Base.iterate(::T) where {T <: Algorithm}         = error("Algorithms need to be initialized first.")
iterate(::T, ::SubReport) where {T <: Algorithm} = error("iterate(::Algorithm, ::Subreport) needs" *
                                                         "to be implemented by each algorithm.")

"""
    Short hand notation for conditional execution of algorithms needs to be
    expanded.

    In case we have an algorithm with an execution condition in Tuple form it
    needs to be converted into ConditionalExec before application.
"""
function (new_algorithm_type::Type{T})(args...) where {T<:Algorithm}
    conv_to_alg(i) = args[i] isa Tuple{Algorithm, Function} ?
                        ConditionalExec(args[i][1], args[i][2]) :
                        args[i]

    newargs = ntuple(conv_to_alg, length(args))
    args == newargs && throw(MethodError(new_algorithm_type, args))
    new_algorithm_type(newargs...)
end

"""
    Initialization is not neccessary for all Algorithms. Checks if algorithm
    needs to be initialized and executes the right method if it needs to.
"""
function initialize_if_neccessary(alg::Algorithm, problem::ScfProblem, state::ScfIterState, params::Parameters)
    applicable(initialize, alg, problem, state, params) ? initialize(alg, problem, state, params) : state
end
