mutable struct ConditionalExec <: Algorithm
    algorithm::Algorithm
    condition::Function
end

ConditionalExec(uninit::UninitialisedAlgorithm, condition::Function; params...) = invoke(ConditionalExec, Tuple{Vararg{Any,N} where N}, uninit, condition; params...)

ConditionalExec(problem::ScfProblem, state::ScfIterState, lg::Logger, algorithm::Algorithm, condition::Function; params...) = ConditionalExec(algorithm, condition)

function iterate(ce::ConditionalExec, subreport::SubReport)
    lg = Logger(subreport)

    if ce.condition(subreport)
        log!(lg, "applying algorithm", typeof(ce.algorithm))
        res = iterate(ce.algorithm, subreport)
        isempty(res) && return res

        resalg, newrp = res
        newalg = ConditionalExec(resalg, ce.condition)
        return newalg, new_subreport(newalg, lg, newrp)
    end
    return ce, subreport
end
