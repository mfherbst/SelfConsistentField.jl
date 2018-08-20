mutable struct Report
    problem::ScfProblem
    state::ScfIterState
    algorithm::Algorithm
    history::Vector{SubReport}
    loglevel::Dict{Symbol, Set}
end

mutable struct SubReport
    algorithm::Union{Missing, Algorithm}
    problem::ScfProblem
    state::Union{Missing, ScfIterState}
    messages::Vector{ReportMessage}
    datasource::Union{Nothing, SubReport}
    loglevel::Dict{Symbol, Set}
    report::Report # There are cases in which algorithms need information about
                   # past events and the whole Algorithm setup
end

abstract type Algorithm end

# Algorithms need to be initialized
Base.iterate(::T) = throw("You need to initialize the algorithm first") where T <: Algorithm

mutable struct ChainedAlgorithm <: Algorithm
    algorithms::Vector{Algorithm}
    reports::Vector{SubReport}
    done::Bool
end

function Base.iterate(chainedalg::ChainedAlgorithm, subreport::SubReport)
    chainedalg.done && return nothing

    # track if all algorithms are done
    done = true

    # Copy the subreport for the new iteration
    subsubreport = copy(subreport)
    for algorithm in chainedalg.algorithms
        subsubreport = SubReport(missing, report.problem, missing, algorithm.toactivate(datasource), Vector{ReportEntry}(), subsubreport, subreport.loglevel, report)
        res = iterate(algorithm, subsubreport)

        # if the iteration is not done, reset done variable
        if res != nothing
            done = false
            push!(chainedalg.reports, subsubreport)
        end
    end

    subreport.state = subsubreport.state
    (chainedalg, subreport)
end

function Base.(new_algorithm_type::Type{T})(old_algorithm::Algorithm, options...) where {T<:Algorithm}
    ChainedAlgorithm(old_algorithm, new_algorithm_type(options...))
end

struct Initializer <: Algorithm end

function initialize(algoritm::Algorithm, problem::ScfProblem, guess::AbstractGuess = RandomGuess();
                    softdefaults::Dict{Symbol, Any} = Dict{Symbol, Any}(),
                    harddefaults::Dict{Symbol, Any} = Dict{Symbol, Any}(),
                    loglevel::Dict{Symbol, Set})
    # As long as we do not have a complete report yet we cannot log messages easily
    # Create a fake subreport and use the usual log! function. Copy the messages over later.
    logmessages = Vector{ReportEntry}()
    fakesubreport = SubReport(Initializer, problem, missing, false, logmessages, nothing, Report(problem, ScfIterState([], [], Dict{String, Float64}(), nothing, nothing, nothing), [], [[]], loglevel))
    log!(fakesubreport, "Starting Initialization", :info)


    # Run initial algorithm configuration and start with the guess.

    # The first algorithm needs to be callable with Union{Guess, ScfIterState},
    # the rest may theoretically only accept ScfIterState.
    log!(fakesubreport, "begin initializing algorithms", :debug, :subalginit)
    state::Union{Guess, ScfIterState} = guess
    for subalg in algorithm
        log!(fakesubreport, "initializing " * String(repr("text/plain", Stuff)), :debug, :subalginit)
        state = initialize(subalg, problem, state, softdefaults)
    end
    log!(fakesubreport, "done initializing algorithms", :debug, :subalginit)

    # Hard defaults override the values produced by automatic algorithm
    # configuration
    apply_hard_defaults(algorithm, harddefaults)

    # Construct report and add already existing log messages.
    log!(fakesubreport, "setting up initial report", :debug,)
    report = Report(problem, state, algorithm, Vector{Vector{SubReport}}(), loglevel)
    subreport = SubReport(Initializer, problem, missing, true, logmessages, nothing, report)
    push!(report.history, logmessages)

    # return new report
    report
end

# Iteration functions for algorithms
Base.iterate(report::Report) = report.state
function Base.iterate(report::Report, ::ScfIterState)
    subreport = SubReport(missing, report.problem, missing, Vector{ReportMessage}(), report.history[end], report.loglevel, report)
    res = iterate(report.algorithm, subreport)
    
    # If the underlying algorithm is done, we are done :-)
    res == nothing && return nothing
    
    push!(report.history, res[2])
    return report, res[2]
end

# Logging
mutable struct ReportMessage
    msg::String
    data::Any
    location::Symbol
end

ReportMessage(msg::String, data:Any) = ReportMessage(msg, data, :memory)
ReportMessage(msg::String) = ReportMessage(msg, nothing, :memory)

function log!(subreport::SubReport, msg::String, data::Any, level::Symbol...)
    if any(map(x -> x ∈ subreport.report.loglevel[:report], level))
        push!(subreport.messages, ReportMessage(msg, data, :memory))
    end
    if any(map(x -> x ∈ subreport.report.loglevel[:stdout], level))
        print(level[end], ': ', msg, ' ')
        show(stdout, MIME("text/plain"), data)
        println()
    end
    data
end


function log!(subreport::SubReport, msg::String, level::Symbol...)
    log!(subreport, msg, nothing, level)
end

function logger(subreport::SubReport)
    args -> log!(subreport, args)
end
