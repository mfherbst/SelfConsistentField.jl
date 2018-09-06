function compute_error_norm(::ScfProblem, iterate::Iterate)
    norm(reshape(iterate.error_pulay, :))
end
