mutable struct Barrier <: Algorithm
    algorithm::Algorithm
    future_algorithm::Algorithm
    condition::Function
    changed::Bool
end

function Barrier(::ScfProblem, state::ScfIterState, lg::Logger, alg1::Algorithm, alg2::Algorithm, switchcondition::Function; params...)
    Barrier(alg1, alg2, switchcondition, false)
end

function notify(barrier::Barrier, subreport::SubReport)
    algorithm = barrier.algorithm
    future_algorithm = barrier.future_algorithm

    if applicable(notify, barrier.algorithm, subreport)
        algorithm, subreport = notify(barrier.algorithm, subreport)
    end

    if !barrier.changed
        if applicable(notify, barrier.future_algorithm, subreport)
            future_algorithm, subreport = notify(barrier.future_algorithm, subreport)
        end
    end
    newbr = Barrier(algorithm, future_algorithm, barrier.condition, barrier.changed)
    return newbr, SubReport(newbr, subreport)
end

function iterate(barrier::Barrier, subreport::SubReport)
    changed = barrier.changed
    algorithm = barrier.algorithm
    future_algorithm = barrier.future_algorithm
    lg = Logger(subreport)

    if !changed
        if barrier.condition(subreport)
            changed = true
            algorithm = barrier.future_algorithm
            log!(lg, "Switching to algorithm", typeof(barrier.future_algorithm), :info, :changealgorithm)
        else
            if applicable(notify, barrier.future_algorithm, subreport)
                future_algorithm, subreport = notify(barrier.future_algorithm, subreport)
            end
        end
    end

    res = iterate(algorithm, subreport)
    res == nothing && return res

    resalgorithm, newsubreport = res
    newalg = Barrier(resalgorithm, future_algorithm, barrier.condition, changed)
    return newalg, SubReport(newalg, lg, newsubreport)
end
