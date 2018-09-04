mutable struct Barrier <: Algorithm
    algorithm::Algorithm
    future_algorithm::Algorithm
    condition::Function
    changed::Bool
end

function Barrier(::ScfProblem, state::ScfIterState, lg::Logger, alg1::Algorithm, alg2::Algorithm, switchcondition::Function; params...)
    Barrier(alg1, alg2, switchcondition, false)
end

function iterate(barrier::Barrier, subreport::SubReport)
    changed = barrier.changed
    algorithm = barrier.algorithm
    lg = Logger(subreport)

    if !changed && barrier.condition(subreport)
            changed = true
            algorithm = barrier.future_algorithm
            log!(lg, "Switching to algorithm", typeof(barrier.future_algorithm), :info, :changealgorithm)
    end

    res = iterate(algorithm, subreport)
    res == nothing && return res

    resalgorithm, newsubreport = res
    newalg = Barrier(resalgorithm, barrier.future_algorithm, barrier.condition, changed)
    return newalg, new_subreport(newalg, lg, newsubreport)
end
