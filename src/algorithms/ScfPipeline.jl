mutable struct ScfPipeline <: Algorithm
    algorithms::Vector{Algorithm}
    subreports::Vector{SubReport}
end

function ScfPipeline(problem::ScfProblem, state::ScfIterState, lg::Logger, algorithms::Algorithm...)
    ScfPipeline(collect(algorithms), Vector{SubReport}())
end

function copy(sp::ScfPipeline)
    ScfPipeline(map(copy, sp.algorithms), copy(sp.subreports))
end

function iterate(scfpipeline::ScfPipeline, rp::SubReport)
    lg = Logger(rp)
    pipe = copy(scfpipeline)

    if !ismissing(rp.convergence)
        rp.convergence.is_converged && return nothing
    end

    # Copy the subreport for the new iteration
    subsubreport = rp
    for algorithm in pipe.algorithms
        log!(lg, "Applying Algorithm ", typeof(algorithm), :debug)
        res = iterate(algorithm, subsubreport)

        # if the iteration is not done, reset done variable
        if res != nothing
            _, subsubreport = res
            push!(pipe.subreports, subsubreport)

            if !ismissing(subsubreport.convergence) ? subsubreport.convergence.is_converged : false
                log!(subsubreport, "Convergence reached", :debug)
                break
            end
        else
            return nothing
        end
    end

    pipe, new_subreport(pipe, subsubreport, lg, rp)
end
