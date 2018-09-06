mutable struct ConditionalExec <: Algorithm
    algorithm::Algorithm
    condition::Function
end

ConditionalExec(uninit::UninitialisedAlgorithm, condition::Function; params...) = invoke(ConditionalExec, Tuple{Vararg{Any,N} where N}, uninit, condition; params...)

ConditionalExec(problem::ScfProblem, iterate::Iterate, lg::Logger, algorithm::Algorithm, condition::Function; params...) = ConditionalExec(algorithm, condition)

function notify(ce::ConditionalExec, stepstate::StepState)
    if applicable(notify, ce.algorithm, stepstate)
        algorithm, stepstate = notify(ce.algorithm, stepstate)
        newce = ConditionalExec(algorithm, ce.condition)
        return newce, stepstate
    else
        return ce, stepstate
    end
end

function iterate(ce::ConditionalExec, stepstate::StepState)
    lg = Logger(stepstate)

    if ce.condition(stepstate)
        log!(lg, "applying algorithm", typeof(ce.algorithm))
        res = iterate(ce.algorithm, stepstate)
        isempty(res) && return res

        resalg, newrp = res
        newalg = ConditionalExec(resalg, ce.condition)
        return newalg, StepState(newalg, lg, newrp)
    else
        if applicable(notify, ce.algorithm, stepstate)
            algorithm, stepstate = notify(ce.algorithm, stepstate)
            newce = ConditionalExec(algorithm, ce.condition)
            return newce, stepstate
        end
    end
    return ce, stepstate
end
