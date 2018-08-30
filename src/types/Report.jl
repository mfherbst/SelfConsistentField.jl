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

struct Logger
    loglevel::LogLevel
    messages::Vector{ReportMessage}
end
Logger(loglevel::LogLevel) = Logger(loglevel, Vector{ReportMessage}())
Logger(rp::SubReport) = Logger(rp.loglevel, Vector{ReportMessage}())
Logger(logger::Logger) = Logger(logger.loglevel, Vector{ReportMessage}())
