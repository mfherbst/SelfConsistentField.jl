mutable struct Barrier <: Algorithm
    algorithm1::Algorithm
    algorithm2::Algorithm
    changecondition::Function
    changed::Bool
    Barrier(alg1::Algorithm, alg2::Algorithm, condition::Function) = new(alg1, alg2, condition, false)
end

function setup(barrier::Barrier, problem::ScfProblem, state::ScfIterState, params::Parameters)
    alg1 = setup_if_neccessary(barrier.algorithm1, problem, state, params)
    alg2 = setup_if_neccessary(barrier.algorithm2, problem, state, params)
    Barrier(alg1, alg2, barrier.changecondition, barrier.changed)
end

function iterate(barrier::Barrier, subreport::SubReport)
    if !barrier.changed
        if barrier.changecondition(subreport)
            barrier.changed = true
            log!(subreport, "Switching to algorithm", typeof(barrier.algorithm2), :info, :changealgorithm)
        end
    end
    if !barrier.changed
        barrier.algorithm1, rp = iterate(barrier.algorithm1, subreport)
    else
        barrier.algorithm2, rp = iterate(barrier.algorithm2, subreport)
    end
    return barrier, rp
end
