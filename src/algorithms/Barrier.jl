mutable struct Barrier <: Algorithm
    algorithm::Algorithm
    future_algorithm::Algorithm
    condition::Function
    changed::Bool
end

function Barrier(::ScfProblem, iterate::Iterate, lg::Logger, alg1::Algorithm, alg2::Algorithm, switchcondition::Function; params...)
    Barrier(alg1, alg2, switchcondition, false)
end

function notify(barrier::Barrier, stepstate::StepState)
    algorithm = barrier.algorithm
    future_algorithm = barrier.future_algorithm

    if applicable(notify, barrier.algorithm, stepstate)
        algorithm, stepstate = notify(barrier.algorithm, stepstate)
    end

    if !barrier.changed
        if applicable(notify, barrier.future_algorithm, stepstate)
            future_algorithm, stepstate = notify(barrier.future_algorithm, stepstate)
        end
    end
    newbr = Barrier(algorithm, future_algorithm, barrier.condition, barrier.changed)
    return newbr, StepState(newbr, stepstate)
end

function iterate(barrier::Barrier, stepstate::StepState)
    changed = barrier.changed
    algorithm = barrier.algorithm
    future_algorithm = barrier.future_algorithm
    lg = Logger(stepstate)

    if !changed
        if barrier.condition(stepstate)
            changed = true
            algorithm = barrier.future_algorithm
            log!(lg, "Switching to algorithm", typeof(barrier.future_algorithm), :info, :changealgorithm)
        else
            if applicable(notify, barrier.future_algorithm, stepstate)
                future_algorithm, stepstate = notify(barrier.future_algorithm, stepstate)
            end
        end
    end

    res = iterate(algorithm, stepstate)
    res == nothing && return res

    resalgorithm, newstepstate = res
    newalg = Barrier(resalgorithm, future_algorithm, barrier.condition, changed)
    return newalg, StepState(newalg, lg, newstepstate)
end
