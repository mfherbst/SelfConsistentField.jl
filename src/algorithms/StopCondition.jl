mutable struct StopCondition <: Algorithm
    algorithm::Algorithm
    condition::Function
end

function initialize(sc::StopCondition, problem::ScfProblem, state::ScfIterState, params::Parameters)
    initialize_if_neccessary(sc.algorithm, problem, state, params)
end

function iterate(sc::StopCondition, subreport::SubReport)
    sc.condition(subreport) ? iterate(sc.algorithm, subreport) : nothing
end