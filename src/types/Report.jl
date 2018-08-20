mutable struct Report
    problem::ScfProblem
    state::ScfIterState
    algorithm::Algorithm
    history::Vector{Vector{SubReport}}
    loglevel::Dict{Symbol, Set}
end

mutable struct SubReport
    subalgorithm::SubAlgorithm
    problem::ScfProblem
    state::Union{Missing, ScfIterState}
    active::Bool
    messages::Vector{ReportEntry}
    datasource::Union{Nothing, SubReport}
    report::Report # There are cases in which algorithms need information about
                   # past events and the current main algorithm
end

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

# An Algorithm is a list of (SubAlgorithm, Function) pairs for which
# apply(SubAlgorithm, state::SubReport) is called in order of occurence if
# Function applied to the previous SubReport returns true.
const Algorithm = Vector{Tuple{alg::SubAlgorithm, toactivate::Function}}
abstract type SubAlgorithm end

function Base.(new_subalgorithm_type::Type{T})(old_subalgorithm::SubAlgorithm, options...) where {T<:SubAlgorithm}
    [old_subalgorithm, new_subalgorithm_type(options...)]
end

function Base.(new_subalgorithm_type::SubAlgorithm)(old_algorithm::Algorithm, options...)
    pushfirst!(old_algorithm, new_subalgorithm_type(options...))
end

function Base.(subalgorithm::SubAlgorithm)(subreport::SubReport)
    push!(subreport.report.history[end], apply(subalgorithm, subreport))
end

struct Initializer <: SubAlgorithm end
notify(::Initializer, ::SubReport) = false
apply(::Initializer, ::SubReport) = throw(InvalidStateException)

function initialize(algoritm::Algorithm, problem::ScfProblem, guess::AbstractGuess = RandomGuess();
                    softdefaults::Dict{Symbol, Any} = Dict{Symbol, Any}(), harddefaults::Dict{Symbol, Any} = Dict{Symbol, Any}(),
                    loglevel::Dict{Symbol, Set})
    # As long as we do not have a complete report yet we cannot log messages easily
    # Create a fake subreport and use the usual log! function. Copy the messages over later.
    logmessages = Vector{ReportEntry}()
    fakesubreport = SubReport(Initializer, problem, missing, false, logmessages, nothing, Report(problem, ScfIterState([], [], Dict{String, Float64}(), nothing, nothing, nothing), [], [[]], loglevel))
    log!(fakesubreport, "Starting Initialization", :info)


    # Run initial subalgorithm configuration and start with the guess.

    # The first algorithm needs to be callable with Union{Guess, ScfIterState},
    # the rest may theoretically only accept ScfIterState.
    log!(fakesubreport, "begin initializing subalgorithms", :debug, :subalginit)
    state::Union{Guess, ScfIterState} = guess
    for subalg in algorithm
        log!(fakesubreport, "initializing " * String(repr("text/plain", Stuff)), :debug, :subalginit)
        state = initialize(subalg, problem, state, softdefaults)
    end
    log!(fakesubreport, "done initializing subalgorithms", :debug, :subalginit)

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
Base.iterate(report::Report) = report
function Base.iterate(report::Report, ::ScfIterState)
    subreport = report.history[end][end]
    # if no subalgorithm was active in the last iteration, return nothing to abort iteration
    if all(map(sub -> !sub.active, report.history[end]))
        log!(subreport, "breaking iteration since no subalgorithms were active in the last one", :info, :breakmessage)
        log!(subreport, "breaking state", report, :debug, :toplevelbreakstate)
        return nothing
    else
        push!(report.history, Vector{SubReport}())
        for subalgorithm in report.algorithm
            subreport = SubReport(subalgorithm.alg, report.problem, missing, subalgorithm.toactivate(datasource), Vector{ReportEntry}(), subreport, report)
            if applicable(notify, subalgorithm.alg, subreport)
                notify(subalgorithm.alg, subreport)
            end
            if subreport.active
                push!(apply(subreport.subalgorithm, subreport), report.history[end])
            end
        end
    end
end
