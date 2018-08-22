# Logging
function log!(rp::SubReport, msg::String, data::Any, level::Symbol...)
    if haskey(rp.loglevel, :report)
        if !isempty(level ∩ rp.loglevel[:report])
            push!(rp.messages, ReportMessage(msg, data, :memory))
        end
    end
    if haskey(rp.loglevel, :stdout)
        if !isempty(level ∩ rp.loglevel[:stdout])
            print(level[end], ": ", msg, " ")
            data != nothing && show(stdout, MIME("text/plain"), data)
            println()
        end
    end
    data
end

function log!(rp::SubReport, msg::String, level::Symbol...)
    log!(rp, msg, nothing, level...)
end

function logger(rp::SubReport)
    args -> log!(rp, args...)
end

# During initialization no functional subreport is available yet. This function
# returns a logging function with a fake subreport which supports everything
# the logging function needs.
function logger(logmessages::Vector{ReportMessage}, loglevel::Dict{Symbol, Set})
    #fakesubreport = SubReport(Initializer(), ScfProblem(), missing, false, logmessages, nothing, Report(problem, FockIterState([], [], Dict{String, Float64}(), nothing, nothing, nothing), [], [[]], loglevel))
    # TODOneeds to be implemented somehow
    (args...) -> println(args)
end
