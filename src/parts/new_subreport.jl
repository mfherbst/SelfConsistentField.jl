function new_subreport(algorithm::Algorithm, rp::SubReport)
    SubReport(algorithm, rp.problem, missing, rp.convergence, Vector{ReportMessage}(), rp, rp.loglevel)
end

function new_subreport(rp::SubReport)
    new_subreport(copy(rp.algorithm), rp)
end
