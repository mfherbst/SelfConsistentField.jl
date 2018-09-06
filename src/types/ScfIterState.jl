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
    error_norm::Float64          # Frobenius error norm of the Pulay error
    energy_change::Dict{String, Float64}   # Change in the energy terms
    is_converged::Bool           # Is SCF converged
    messages::Vector{LogMessage}
    loglevel::LogLevel
    previous::Union{Nothing, StepState}
end

mutable struct ScfIterState
    problem::ScfProblem
    iterate::Iterate
    error_norm::Float64          # Frobenius error norm of the Pulay error
    energy_change::Dict{String, Float64}   # Change in the energy terms
    is_converged::Bool           # Is SCF converged
    algorithm::Algorithm
    loglevel::LogLevel
    history::Vector{StepState}
end

struct Logger
    loglevel::LogLevel
    messages::Vector{LogMessage}
end
Logger(loglevel::LogLevel) = Logger(loglevel, Vector{LogMessage}())
Logger(rp::StepState) = Logger(rp.loglevel, Vector{LogMessage}())
Logger(logger::Logger) = Logger(logger.loglevel, Vector{LogMessage}())

function StepState(algorithm::Algorithm, iterate::Iterate, logger::Logger, rp::StepState)
    StepState(algorithm, rp.problem, iterate, rp.error_norm, rp.energy_change, rp.is_converged, logger.messages, logger.loglevel, rp)
end

function StepState(algorithm::Algorithm, logger::Logger, rp::StepState)
    StepState(algorithm, rp.iterate, logger, rp)
end

function StepState(algorithm::Algorithm, rp::StepState)
    StepState(algorithm, Logger(rp), rp)
end
