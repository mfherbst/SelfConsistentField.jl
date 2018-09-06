struct ConvergenceCheck <: Algorithm
    convergencecondition::Function
    olditerate::Union{Missing, Iterate}
end

function ConvergenceCheck(::ScfProblem, iterate::Iterate, lg::Logger, convergencecondition::Function; store_initial_matrix = true, params...)
    ConvergenceCheck(convergencecondition, store_initial_matrix ? iterate : missing)
end

copy(cc::ConvergenceCheck) = ConvergenceCheck(cc.convergencecondition, cc.olditerate)

function ConvergenceCheck(problem::ScfProblem, iterate::Iterate, lg::Logger;
                          max_error_norm::Number = 5e-7, max_energy_total_change::Number = 1.25e-07,
                          max_energy_1e_change::Number = 5e-5, params...)

    function convergencecondition(error_norm::Float64, energy_change::Dict{String, Float64})
        !(error_norm >= max_error_norm ||
          norm(energy_change["total"]) >= max_energy_total_change ||
          norm(energy_change["1e"]) >= max_energy_1e_change)
    end
    ConvergenceCheck(problem, iterate, lg, convergencecondition)
end

function iterate(cc::ConvergenceCheck, rp::StepState)
    newcc = ConvergenceCheck(cc.convergencecondition, rp.iterate)
    lg = Logger(rp)

    rp.is_converged && return nothing

    if ismissing(cc.olditerate)
        return newcc, StepState(newcc, rp)
    else
        error_norm = compute_error_norm(rp.problem, rp.iterate)
        energy_change = Dict{String, Float64}()
        for (key, oval) in cc.olditerate.energies
            nval = rp.iterate.energies[key]
            val = nval - oval
            energy_change[key] = val
        end
        is_converged = cc.convergencecondition(error_norm, energy_change)

        return newcc, StepState(newcc, rp.problem, rp.iterate, error_norm, energy_change, is_converged, lg.messages, rp.loglevel, rp)
    end
end
