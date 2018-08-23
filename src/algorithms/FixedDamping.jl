#
# Simple damping
#
mutable struct FixedDamping <: Algorithm
    damping::Union{Missing, Float64}                   # Damping coefficient

    prev_iterate::Union{Nothing,ScfIterState}  # Previous iterate
    function FixedDamping(; damping_value = missing)
        new(damping_value, nothing)
    end
end

function initialize(fd::FixedDamping, problem::ScfProblem, iterstate::ScfIterState, params::Parameters)
    # TODO needs to be implemented the same way in all other algorithms
    fd.damping = coalesce(fd.damping, get(params, :damping_value, missing), 0.4)
    return iterstate
end

function iterate(fd::FixedDamping, subreport::SubReport)
    rp = new_subreport(subreport)
    rp.state = compute_next_iterate(fd, rp.source.state)
    log!(rp, "inside FixedDamping", :debug)
    return fd, rp
end

function compute_next_iterate(acc::FixedDamping, iterate::ScfIterState)
    if acc.prev_iterate != nothing
        prev_mat = get_iterate_matrix(acc.prev_iterate)
        cur_mat = get_iterate_matrix(iterate)
        damp = acc.damping
        ret = (1 - damp) * prev_mat + damp * cur_mat
        iterate = update_iterate_matrix(iterate, ret)
        acc.prev_iterate = iterate
        return iterate
    else
        acc.prev_iterate = iterate
        return iterate
    end
end
