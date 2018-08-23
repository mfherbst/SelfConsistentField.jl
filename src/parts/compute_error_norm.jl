function compute_error_norm(::ScfProblem, state::ScfIterState)
    norm(reshape(state.error_pulay, :))
end
