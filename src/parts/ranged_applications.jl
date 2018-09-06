function before_errnorm(threshold::Number)
    function enable_when(rp::StepState)
        !ismissing(rp.convergence) ? rp.convergence.error_norm > threshold : true
    end
end

function after_errnorm(threshold::Number)
    function enable_when(rp::StepState)
        !ismissing(rp.convergence) ? rp.convergence.error_norm < threshold : false
    end
end

function between_errnorm(start::Number, stop::Number)
    function enable_when(rp::StepState)
        !ismissing(rp.convergence) ? start > rp.convergence.error_norm > stop : false
    end
end

