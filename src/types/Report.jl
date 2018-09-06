abstract type Algorithm end

const LogLevel = Dict{Symbol, Set}
const Parameters = Dict{Symbol, Any}

mutable struct LogMessage
    msg::String
    data::Any
end
LogMessage(msg::String) = LogMessage(msg, nothing)

mutable struct SubReport
    algorithm::Union{Missing, Algorithm}
    problem::ScfProblem
    state::Union{Missing, ScfIterState}
    convergence::Union{Missing, ScfConvergence}
    messages::Vector{LogMessage}
    previous::Union{Nothing, SubReport}
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
    messages::Vector{LogMessage}
end
Logger(loglevel::LogLevel) = Logger(loglevel, Vector{LogMessage}())
Logger(rp::SubReport) = Logger(rp.loglevel, Vector{LogMessage}())
Logger(logger::Logger) = Logger(logger.loglevel, Vector{LogMessage}())

function SubReport(algorithm::Algorithm, state::ScfIterState, logger::Logger, rp::SubReport)
    SubReport(algorithm, rp.problem, state, rp.convergence, logger.messages, rp, logger.loglevel)
end

function SubReport(algorithm::Algorithm, logger::Logger, rp::SubReport)
    SubReport(algorithm, rp.state, logger, rp)
end

function SubReport(algorithm::Algorithm, rp::SubReport)
    SubReport(algorithm, Logger(rp), rp)
end
