mutable struct ScfPipeline <: Algorithm
    algorithms::Vector{Algorithm}
    subreports::Vector{StepState}
end

function ScfPipeline(::ScfProblem, ::Iterate, lg::Logger, algorithms::Algorithm...)
    log!(lg, "constructing ScfPipeline", :debug, :scfpipeline, :setup)
    ScfPipeline(collect(algorithms), Vector{StepState}())
end

ScfPipeline(uninit1::UninitialisedAlgorithm, uninit2::UninitialisedAlgorithm; params...) = invoke(ScfPipeline, Tuple{Vararg{Any,N} where N}, uninit1, uninit2; params...)

function notify(sp::ScfPipeline, subreport::StepState)
    new_algorithms = Vector{Algorithm}()
    new_subreports = Vector{StepState}()
    rp = subreport
    for algorithm in sp.algorithms
        if applicable(notify, ce.algorithm, subreport)
            alg, rp = notify(algorithm, rp)
            push!(new_algorithms, alg)
            push!(new_subreports, rp)
        end
    end
    newsp = ScfPipeline(new_algorithms, new_subreports)
    return newsp, StepState(newsp, rp)
end

function iterate(scfpipeline::ScfPipeline, rp::StepState)
    lg = Logger(rp)
    newpipe = ScfPipeline(Vector{Algorithm}(), Vector{StepState}())

    !ismissing(rp.convergence) && rp.convergence.is_converged && return nothing

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

            if !ismissing(subsubreport.convergence) && subsubreport.convergence.is_converged
                log!(subsubreport, "Convergence reached", :debug, :scfpipeline)
                break
            end
        else
            return nothing
        end
    end

    newpipe, StepState(newpipe, lg, subsubreport)
end
