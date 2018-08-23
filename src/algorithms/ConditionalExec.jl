mutable struct ConditionalExec <: Algorithm
    algorithm::Algorithm
    condition::Function
end

function initialize(ce::ConditionalExec, problem::ScfProblem, state::ScfIterState, params::Parameters)
    initialize_if_neccessary(ce.algorithm, problem, state, params)
end

function iterate(ce::ConditionalExec, subreport::SubReport)
    ce.condition(subreport) ? iterate(ce.algorithm, subreport) : subreport
end
