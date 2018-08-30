abstract type Algorithm end

const LogLevel = Dict{Symbol, Set}
const Parameters = Dict{Symbol, Any}

mutable struct SubReport
    algorithm::Union{Missing, Algorithm}
    problem::ScfProblem
    state::Union{Missing, ScfIterState}
    convergence::Union{Missing, ScfConvergence}
    messages::Vector{ReportMessage}
    source::Union{Nothing, SubReport}
    loglevel::LogLevel
end

mutable struct Report
    problem::ScfProblem
    state::ScfIterState
    convergence::Union{Missing, ScfConvergence}
    algorithm::Algorithm
    history::Vector{SubReport}
    loglevel::LogLevel
end

struct InitReport
    loglevel::LogLevel
    messages::Vector{ReportMessage}
end
InitReport(loglevel::LogLevel) = InitReport(loglevel, Vector{ReportMessage}())
InitReport(initrp::InitReport) = InitReport(initrp.loglevel, Vector{ReportMessage}())
