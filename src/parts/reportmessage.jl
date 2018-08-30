# Logging
function log!(rp::Union{Logger, SubReport}, msg::String, data::Any, level::Symbol...)
    if haskey(rp.loglevel, :report)
        if !isempty(level ∩ rp.loglevel[:report])
            push!(rp.messages, ReportMessage(msg, data))
        end
    end
    if haskey(rp.loglevel, :stdout)
        if !isempty(level ∩ rp.loglevel[:stdout])
            if data != nothing
                print(level[end], ": ", msg, " ")
                show(stdout, MIME("text/plain"), data)
                println()
            else
                println(level[end], ": ", msg)
            end
        end
    end
    data
end

function log!(rp::Union{Logger, SubReport}, msg::String, level::Symbol...)
    log!(rp, msg, nothing, level...)
end
