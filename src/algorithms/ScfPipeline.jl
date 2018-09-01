mutable struct ScfPipeline <: Algorithm
    algorithms::Vector{Algorithm}
    subreports::Vector{SubReport}
end

function ScfPipeline(::ScfProblem, ::ScfIterState, lg::Logger, algorithms::Algorithm...)
    log!(lg, "constructing ScfPipeline", :debug, :scfpipeline, :setup)
    ScfPipeline(collect(algorithms), Vector{SubReport}())
end

function copy(sp::ScfPipeline)
    ScfPipeline(map(copy, sp.algorithms), Base.copy(sp.subreports))
end

function iterate(scfpipeline::ScfPipeline, rp::SubReport)
    lg = Logger(rp)
    newpipe = ScfPipeline(Vector{Algorithm}(), Vector{SubReport}())

    if !ismissing(rp.convergence)
        rp.convergence.is_converged && return nothing
    end

    # Copy the subreport for the new iteration
    subsubreport = rp
    for algorithm in scfpipeline.algorithms
        log!(lg, "Applying Algorithm ", typeof(algorithm), :debug)
        res = iterate(algorithm, subsubreport)

        # if the iteration is not done, reset done variable
        if res != nothing
            newalgorithm, subsubreport = res
            push!(newpipe.algorithms, newalgorithm)
            push!(newpipe.subreports, subsubreport)

            if !ismissing(subsubreport.convergence) ? subsubreport.convergence.is_converged : false
                log!(subsubreport, "Convergence reached", :debug)
                break
            end
        else
            return nothing
        end
    end

    newpipe, new_subreport(newpipe, subsubreport.state, lg, rp)
end
