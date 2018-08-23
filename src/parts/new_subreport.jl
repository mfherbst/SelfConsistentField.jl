function new_subreport(rp::SubReport)
    SubReport(rp.algorithm, rp.problem, missing, rp.convergence, Vector{ReportMessage}(), rp, rp.loglevel, rp.report)
end
