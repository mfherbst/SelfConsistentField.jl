#
# Simple damping
#
mutable struct FixedDamping <: Accelerator
    damping::Float64                   # Damping coefficient

    prev_iterate::Union{Nothing,ScfIterState}  # Previous iterate
    function FixedDamping(problem::ScfProblem; damping_value=0.4, kwargs...)
        new(damping_value, nothing)
    end
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
