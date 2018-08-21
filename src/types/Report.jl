abstract type Algorithm end

mutable struct ReportMessage
    msg::String
    data::Any
    location::Symbol
end
ReportMessage(msg::String, data::Any) = ReportMessage(msg, data, :memory)
ReportMessage(msg::String) = ReportMessage(msg, nothing, :memory)

# There is a recursive dependency between a Report and Subreport. We break this
# dependency by creating a constructor that automatically sets T as SubReport

# In Julia 1.1 type redifinition should work. If you live in the future, please implement
#    mutable struct SubReport end
# and remove the current workaround
mutable struct Report{T}
    problem::ScfProblem
    state::ScfIterState
    algorithm::Algorithm
    history::Vector{T}
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

Report(p,s,a,h,l) = Report{SubReport}(p,s,a,h,l)
