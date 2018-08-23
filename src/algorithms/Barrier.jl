mutable struct Barrier <: Algorithm
    algorithm1::Algorithm
    algorithm2::Algorithm
    changecondition::Function
    changed::Bool
    Barrier(alg1::Algorithm, alg2::Algorithm, condition::Function) = new(alg1, alg2, condition, false)
end

function initialize(ca::Barrier, problem::ScfProblem, state::ScfIterState, params::Parameters)
    initialize_if_neccessary(ca.algorithm1, problem, state, params)
    initialize_if_neccessary(ca.algorithm2, problem, state, params)
end

function iterate(ca::Barrier, subreport::SubReport)
    if !ca.changed
        if ca.changecondition(subreport)
            ca.changed = true
            log!(subreport, "Switching to algorithm", typeof(ca.algorithm2), :info, :changealgorithm)
        end
    end
    ca.changed ? iterate(ca.algorithm2, subreport) : iterate(ca.algorithm1, subreport)
end
