mutable struct ConditionalExec <: Algorithm
    algorithm::Algorithm
    condition::Function
end

function setup(ce::ConditionalExec, problem::ScfProblem, state::ScfIterState, params::Parameters)
    newalg = setup_if_neccessary(ce.algorithm, problem, state, params)
    ConditionalExec(newalg, ce.condition)
end

function iterate(ce::ConditionalExec, subreport::SubReport)
    ce.algorithm, rp = ce.condition(subreport) ? iterate(ce.algorithm, subreport) : (ce.algorithm, subreport)
    return ce, rp
end
