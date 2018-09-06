#
# Simple damping
#
struct FixedDamping <: Algorithm
    damping::Union{Missing, Float64}                   # Damping coefficient

    prev_iterate::Union{Nothing,Iterate}  # Previous iterate
end

function FixedDamping(problem::ScfProblem, iteriterate::Iterate, lg::Logger; damping=0.4, params...)
    FixedDamping(damping, iteriterate)
end

function iterate(fd::FixedDamping, rp::StepState)
    lg = Logger(rp)
    log!(lg, "running FixedDamping", :debug)


    prev_mat = fd.prev_iterate.fock
    cur_mat = rp.iterate.fock
    ret = (1 - fd.damping) * prev_mat + fd.damping * cur_mat
    iterate = update_fock_matrix(rp.iterate, ret)

    newdamping = FixedDamping(fd.damping, iterate)
    return newdamping, StepState(newdamping, iterate, lg, rp)
end
