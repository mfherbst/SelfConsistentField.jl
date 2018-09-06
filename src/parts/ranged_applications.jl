function before_errnorm(threshold::Number)
    function enable_when(rp::StepState)
        rp.error_norm > threshold
    end
end

function after_errnorm(threshold::Number)
    function enable_when(rp::StepState)
        rp.error_norm < threshold
    end
end

function between_errnorm(start::Number, stop::Number)
    function enable_when(rp::StepState)
        start > rp.error_norm > stop
    end
end

