# Logging
function log!(rp::Union{Logger, StepState}, msg::String, data::Any, level::Symbol...)
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

function log!(rp::Union{Logger, StepState}, msg::String, level::Symbol...)
    log!(rp, msg, nothing, level...)
end
