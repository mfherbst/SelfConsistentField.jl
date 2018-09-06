#
# Simple damping
#
struct FixedDamping <: Algorithm
    damping::Union{Missing, Float64}                   # Damping coefficient

    prev_iterate::Union{Nothing,ScfIterState}  # Previous iterate
end

function FixedDamping(problem::ScfProblem, iterstate::ScfIterState, lg::Logger; damping=0.4, params...)
    FixedDamping(damping, iterstate)
end

function iterate(fd::FixedDamping, rp::StepState)
    lg = Logger(rp)
    log!(lg, "running FixedDamping", :debug)


    prev_mat = get_iterate_matrix(fd.prev_iterate)
    cur_mat = get_iterate_matrix(rp.state)
    ret = (1 - fd.damping) * prev_mat + fd.damping * cur_mat
    iterate = update_iterate_matrix(rp.state, ret)

    newdamping = FixedDamping(fd.damping, iterate)
    return newdamping, StepState(newdamping, iterate, lg, rp)
end
