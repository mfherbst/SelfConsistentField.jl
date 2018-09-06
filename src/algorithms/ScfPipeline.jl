mutable struct ScfPipeline <: Algorithm
    algorithms::Vector{Algorithm}
    stepstates::Vector{StepState}
end

function ScfPipeline(::ScfProblem, ::Iterate, lg::Logger, algorithms::Algorithm...)
    log!(lg, "constructing ScfPipeline", :debug, :scfpipeline, :setup)
    ScfPipeline(collect(algorithms), Vector{StepState}())
end

ScfPipeline(uninit1::UninitialisedAlgorithm, uninit2::UninitialisedAlgorithm; params...) = invoke(ScfPipeline, Tuple{Vararg{Any,N} where N}, uninit1, uninit2; params...)

function notify(sp::ScfPipeline, stepstate::StepState)
    new_algorithms = Vector{Algorithm}()
    new_stepstates = Vector{StepState}()
    rp = stepstate
    for algorithm in sp.algorithms
        if applicable(notify, ce.algorithm, stepstate)
            alg, rp = notify(algorithm, rp)
            push!(new_algorithms, alg)
            push!(new_stepstates, rp)
        end
    end
    newsp = ScfPipeline(new_algorithms, new_stepstates)
    return newsp, StepState(newsp, rp)
end

function iterate(scfpipeline::ScfPipeline, rp::StepState)
    lg = Logger(rp)
    newpipe = ScfPipeline(Vector{Algorithm}(), Vector{StepState}())

    rp.is_converged && return nothing

    # Copy the stepstate for the new iteration
    substepstate = rp
    for algorithm in scfpipeline.algorithms
        log!(lg, "Applying Algorithm ", typeof(algorithm), :debug)
        res = iterate(algorithm, substepstate)

        # if the iteration is not done, reset done variable
        if res != nothing
            newalgorithm, substepstate = res
            push!(newpipe.algorithms, newalgorithm)
            push!(newpipe.stepstates, substepstate)

            if substepstate.is_converged
                log!(substepstate, "Convergence reached", :debug, :scfpipeline)
                break
            end
        else
            return nothing
        end
    end

    newpipe, StepState(newpipe, lg, substepstate)
end
