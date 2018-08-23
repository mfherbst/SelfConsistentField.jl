abstract type Algorithm end

const LogLevel = Dict{Symbol, Set}
const Parameters = Dict{Symbol, Any}

mutable struct ReportMessage
    msg::String
    data::Any
end
ReportMessage(msg::String) = ReportMessage(msg, nothing)

# There is a recursive dependency between a Report and Subreport. We break this
# dependency by creating a constructor that automatically sets T as SubReport

# In Julia 1.1 type redifinition should work. If you live in the future, please implement
#    mutable struct SubReport end
# and remove the current workaround
mutable struct Report{T}
    problem::ScfProblem
    state::ScfIterState
    convergence::Union{Missing, ScfConvergence}
    algorithm::Algorithm
    history::Vector{T}
    loglevel::LogLevel
end

mutable struct SubReport
    algorithm::Union{Missing, Algorithm}
    problem::ScfProblem
    state::Union{Missing, ScfIterState}
    convergence::Union{Missing, ScfConvergence}
    messages::Vector{ReportMessage}
    source::Union{Nothing, SubReport}
    loglevel::LogLevel
    report::Report # There are cases in which algorithms need information about
                   # past events and the whole Algorithm setup
end
# TODO remove reference to report in SubReport

Report(p,s,c,a,h,l) = Report{SubReport}(p,s,c,a,h,l)
