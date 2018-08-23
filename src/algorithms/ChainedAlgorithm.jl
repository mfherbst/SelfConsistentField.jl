mutable struct ChainedAlgorithm <: Algorithm
    algorithms::Vector{Algorithm}
    reports::Vector{SubReport}
    done::Bool
    function ChainedAlgorithm(algorithms::Algorithm...)
        new(collect(algorithms), Vector{SubReport}(), false)
    end
end

function initialize(ca::ChainedAlgorithm, problem::ScfProblem, state::ScfIterState, params::Parameters)
    for algorithm in ca.algorithms
        state = initialize_if_neccessary(algorithm, problem, state, params)
    end
    return state
end

function iterate(chainedalg::ChainedAlgorithm, subreport::SubReport)
    chainedalg.done && return nothing

    # track if all algorithms are done
    done = true

    # Copy the subreport for the new iteration
    subsubreport = subreport
    for algorithm in chainedalg.algorithms
        subsubreport = iterate(algorithm, subsubreport)

        # if the iteration is not done, reset done variable
        if subsubreport != nothing
            done = false
            push!(chainedalg.reports, subsubreport)
        end
    end

    (chainedalg, subsubreport)
end
