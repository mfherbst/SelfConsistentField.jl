mutable struct ConditionalExec <: Algorithm
    algorithm::Algorithm
    condition::Function
end

ConditionalExec(uninit::UninitialisedAlgorithm, condition::Function; params...) = invoke(ConditionalExec, Tuple{Vararg{Any,N} where N}, uninit, condition; params...)

ConditionalExec(problem::ScfProblem, state::ScfIterState, lg::Logger, algorithm::Algorithm, condition::Function; params...) = ConditionalExec(algorithm, condition)

function notify(ce::ConditionalExec, subreport::StepState)
    if applicable(notify, ce.algorithm, subreport)
        algorithm, subreport = notify(ce.algorithm, subreport)
        newce = ConditionalExec(algorithm, ce.condition)
        return newce, subreport
    else
        return ce, subreport
    end
end

function iterate(ce::ConditionalExec, subreport::StepState)
    lg = Logger(subreport)

    if ce.condition(subreport)
        log!(lg, "applying algorithm", typeof(ce.algorithm))
        res = iterate(ce.algorithm, subreport)
        isempty(res) && return res

        resalg, newrp = res
        newalg = ConditionalExec(resalg, ce.condition)
        return newalg, StepState(newalg, lg, newrp)
    else
        if applicable(notify, ce.algorithm, subreport)
            algorithm, subreport = notify(ce.algorithm, subreport)
            newce = ConditionalExec(algorithm, ce.condition)
            return newce, subreport
        end
    end
    return ce, subreport
end
