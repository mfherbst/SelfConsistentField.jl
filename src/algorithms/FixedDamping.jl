#
# Simple damping
#
struct FixedDamping <: Algorithm
    damping::Union{Missing, Float64}                   # Damping coefficient

    prev_state::Union{Nothing,ScfIterState}  # Previous state
end

function FixedDamping(problem::ScfProblem, iterstate::ScfIterState, lg::Logger; damping=0.4, params...)
    FixedDamping(damping, iterstate)
end

function iterate(fd::FixedDamping, rp::StepState)
    lg = Logger(rp)
    log!(lg, "running FixedDamping", :debug)


    prev_mat = fd.prev_state.fock
    cur_mat = rp.state.fock
    ret = (1 - fd.damping) * prev_mat + fd.damping * cur_mat
    state = update_fock_matrix(rp.state, ret)

    newdamping = FixedDamping(fd.damping, state)
    return newdamping, StepState(newdamping, state, lg, rp)
end
