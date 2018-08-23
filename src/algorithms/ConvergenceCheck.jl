mutable struct ConvergenceCheck <: Algorithm
    convergencecondition::Function
    olditerstate::Union{Missing, ScfIterState}
    ConvergenceCheck(convergencecondition::Function) = new(convergencecondition, missing)
end

function ConvergenceCheck(; max_error_norm::Number = 5e-7, max_energy_total_change::Number = 1.25e-07, max_energy_1e_change::Number = 5e-5)
    function convergencecondition(convstate::ScfConvergence)
        (convstate.error_norm >= max_error_norm) |
        (convstate.energy_change["total"] >= max_energy_total_change) |
        (convstate.energy_change["1e"] >= max_energy_1e_change) ?
        false : true
    end
    ConvergenceCheck(convergencecondition)
end

function iterate(cc::ConvergenceCheck, subreport::SubReport)
    rp = new_subreport(subreport)
    rp.state = subreport.state

    if !ismissing(cc.olditerstate)
        error_norm = compute_error_norm(rp.problem, rp.state)
        energy_change = Dict{String, Float64}()
        for (key, oval) in cc.olditerstate.energies
            nval = rp.state.energies[key]
            val = nval - oval
            energy_change[key] = val
        end

        rp.convergence = ScfConvergence(error_norm, energy_change, false)
        rp.convergence.is_converged = cc.convergencecondition(rp.convergence)
    end

    # Store new iterstate for next iteration
    cc.olditerstate = rp.state

    if !ismissing(rp.convergence) ? rp.convergence.is_converged : false
        subreport.convergence = rp.convergence
        return nothing
    end
    return cc, rp
end
