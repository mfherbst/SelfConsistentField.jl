mutable struct ScfPipeline <: Algorithm
    algorithms::Vector{Algorithm}
    reports::Vector{SubReport}
    done::Bool
    function ScfPipeline(algorithms::Algorithm...)
        new(collect(algorithms), Vector{SubReport}(), false)
    end
end

function initialize(sp::ScfPipeline, problem::ScfProblem, state::ScfIterState, params::Parameters)
    for algorithm in sp.algorithms
        state = initialize_if_neccessary(algorithm, problem, state, params)
    end
    return state
end

function iterate(scfpipeline::ScfPipeline, subreport::SubReport)
    scfpipeline.done && return nothing

    # track if all algorithms are done
    done = true

    # Copy the subreport for the new iteration
    subsubreport = subreport
    for algorithm in scfpipeline.algorithms
        log!(subsubreport, "Applying Algorithm ", typeof(algorithm), :debug)
        res = iterate(algorithm, subsubreport)

        # if the iteration is not done, reset done variable
        if res != nothing
            done = false
            _, subsubreport = res
            push!(scfpipeline.reports, subsubreport)
        else
            return nothing
        end
    end

    scfpipeline, subsubreport
end
