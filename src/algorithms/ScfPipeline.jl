mutable struct ScfPipeline <: Algorithm
    algorithms::Vector{Algorithm}
    subreports::Vector{SubReport}
    function ScfPipeline(algorithms::Algorithm...)
        new(collect(algorithms), Vector{SubReport}())
    end
end

function setup(sp::ScfPipeline, problem::ScfProblem, state::ScfIterState, params::Parameters)
    algorithms_setup = Vector{Algorithm}()
    for algorithm in sp.algorithms
        push!(algorithms_setup, setup_if_neccessary(algorithm, problem, state, params))
    end
    return ScfPipeline(algorithms_setup...)
end

function iterate(scfpipeline::ScfPipeline, subreport::SubReport)
    if !ismissing(subreport.convergence)
        subreport.convergence.is_converged && return nothing
    end

    # Copy the subreport for the new iteration
    subsubreport = subreport
    for algorithm in scfpipeline.algorithms
        log!(subsubreport, "Applying Algorithm ", typeof(algorithm), :debug)
        res = iterate(algorithm, subsubreport)

        # if the iteration is not done, reset done variable
        if res != nothing
            _, subsubreport = res
            push!(scfpipeline.subreports, subsubreport)

            if !ismissing(subreport.convergence) ? subreport.convergence.is_converged : false
                log!(subsubreport, "Convergence reached", :debug)
                break
            end
        else
            return nothing
        end
    end

    scfpipeline, subsubreport
end
