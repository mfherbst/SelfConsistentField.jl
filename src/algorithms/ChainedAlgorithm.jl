mutable struct ChainedAlgorithm <: Algorithm
    algorithms::Vector{Algorithm}
    reports::Vector{SubReport}
    done::Bool
end

function initialize(ca::ChainedAlgorithm, problem::ScfProblem, state::ScfIterState, softdefaults::Defaults)
    subalgstate = copy(state)
    for algorithm in ca.algorithms
        subalgstate = initialize_if_neccessary(algorithm, problem, subalgstate, softdefaults)
    end
end

function Base.iterate(chainedalg::ChainedAlgorithm, subreport::SubReport)
    chainedalg.done && return nothing

    # track if all algorithms are done
    done = true

    # Copy the subreport for the new iteration
    subsubreport = copy(subreport)
    for algorithm in chainedalg.algorithms
        subsubreport = SubReport(missing, report.problem, missing, algorithm.toactivate(datasource), Vector{ReportEntry}(), subsubreport, subreport.loglevel, report)
        res = iterate(algorithm, subsubreport)

        # if the iteration is not done, reset done variable
        if res != nothing
            done = false
            push!(chainedalg.reports, subsubreport)
        end
    end

    subreport.state = subsubreport.state
    (chainedalg, subreport)
end
