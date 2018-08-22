function new_subreport(rp::SubReport)
    SubReport(rp.algorithm, rp.problem, missing, Vector{ReportMessage}(), rp, rp.loglevel, rp.report)
end
