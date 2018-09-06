struct ConvergenceCheck <: Algorithm
    convergencecondition::Function
    olditerstate::Union{Missing, ScfIterState}
end

function ConvergenceCheck(::ScfProblem, state::ScfIterState, lg::Logger, convergencecondition::Function; store_initial_matrix = true, params...)
    ConvergenceCheck(convergencecondition, store_initial_matrix ? state : missing)
end

copy(cc::ConvergenceCheck) = ConvergenceCheck(cc.convergencecondition, cc.olditerstate)

function ConvergenceCheck(problem::ScfProblem, state::ScfIterState, lg::Logger;
                          max_error_norm::Number = 5e-7, max_energy_total_change::Number = 1.25e-07,
                          max_energy_1e_change::Number = 5e-5, params...)

    function convergencecondition(convstate::ScfConvergence)
        !(convstate.error_norm >= max_error_norm ||
          norm(convstate.energy_change["total"]) >= max_energy_total_change ||
          norm(convstate.energy_change["1e"]) >= max_energy_1e_change)
    end
    ConvergenceCheck(problem, state, lg, convergencecondition)
end

function iterate(cc::ConvergenceCheck, rp::StepState)
    newcc = ConvergenceCheck(cc.convergencecondition, rp.state)
    lg = Logger(rp)

    !ismissing(rp.convergence) && rp.convergence.is_converged && return nothing

    if ismissing(cc.olditerstate)
        return newcc, StepState(newcc, rp)
    else
        error_norm = compute_error_norm(rp.problem, rp.state)
        energy_change = Dict{String, Float64}()
        for (key, oval) in cc.olditerstate.energies
            nval = rp.state.energies[key]
            val = nval - oval
            energy_change[key] = val
        end

        updated_conv_data = ScfConvergence(error_norm, energy_change, false)
        is_converged = cc.convergencecondition(updated_conv_data)
        convergence = ScfConvergence(updated_conv_data.error_norm, updated_conv_data.energy_change, is_converged)
        log!(lg, "new convergence", convergence, :debug, :convergencecheck)

        return newcc, StepState(newcc, rp.problem, rp.state, convergence, lg.messages, rp, rp.loglevel)
    end
end
