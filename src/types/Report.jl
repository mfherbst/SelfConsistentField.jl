abstract type Algorithm end

const LogLevel = Dict{Symbol, Set}
const Parameters = Dict{Symbol, Any}

mutable struct LogMessage
    msg::String
    data::Any
end
LogMessage(msg::String) = LogMessage(msg, nothing)

mutable struct StepState
    algorithm::Union{Missing, Algorithm}
    problem::ScfProblem
    iterate::Union{Missing, Iterate}
    convergence::Union{Missing, ScfConvergence}
    messages::Vector{LogMessage}
    previous::Union{Nothing, StepState}
    loglevel::LogLevel
end

mutable struct Report
    problem::ScfProblem
    iterate::Iterate
    convergence::Union{Missing, ScfConvergence}
    algorithm::Algorithm
    history::Vector{StepState}
    loglevel::LogLevel
end

struct Logger
    loglevel::LogLevel
    messages::Vector{LogMessage}
end
Logger(loglevel::LogLevel) = Logger(loglevel, Vector{LogMessage}())
Logger(rp::StepState) = Logger(rp.loglevel, Vector{LogMessage}())
Logger(logger::Logger) = Logger(logger.loglevel, Vector{LogMessage}())

function StepState(algorithm::Algorithm, iterate::Iterate, logger::Logger, rp::StepState)
    StepState(algorithm, rp.problem, iterate, rp.convergence, logger.messages, rp, logger.loglevel)
end

function StepState(algorithm::Algorithm, logger::Logger, rp::StepState)
    StepState(algorithm, rp.iterate, logger, rp)
end

function StepState(algorithm::Algorithm, rp::StepState)
    StepState(algorithm, Logger(rp), rp)
end
