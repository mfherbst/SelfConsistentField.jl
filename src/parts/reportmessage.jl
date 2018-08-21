# Logging
function log!(subreport::SubReport, msg::String, data::Any, level::Symbol...)
    if any(map(x -> x ∈ subreport.report.loglevel[:report], level))
        push!(subreport.messages, ReportMessage(msg, data, :memory))
    end
    if any(map(x -> x ∈ subreport.report.loglevel[:stdout], level))
        print(level[end], ": ", msg, " ")
        show(stdout, MIME("text/plain"), data)
        println()
    end
    data
end

function log!(subreport::SubReport, msg::String, level::Symbol...)
    log!(subreport, msg, nothing, level)
end

function logger(subreport::SubReport)
    args -> log!(subreport, args...)
end

# During initialization no functional subreport is available yet. This function
# returns a logging function with a fake subreport which supports everything
# the logging function needs.
function logger(logmessages::Vector{ReportMessage}, loglevel::Dict{Symbol, Set})
    #fakesubreport = SubReport(Initializer(), ScfProblem(), missing, false, logmessages, nothing, Report(problem, FockIterState([], [], Dict{String, Float64}(), nothing, nothing, nothing), [], [[]], loglevel))
    # TODOneeds to be implemented somehow
    (args...) -> println(args)
end
