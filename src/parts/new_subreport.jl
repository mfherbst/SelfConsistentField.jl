function new_subreport(algorithm::Algorithm, state::ScfIterState, logger::Logger, rp::SubReport)
    SubReport(algorithm, rp.problem, state, rp.convergence, logger.messages, rp, logger.loglevel)
end

function new_subreport(algorithm::Algorithm, logger::Logger, rp::SubReport)
    new_subreport(algorithm, rp.state, logger, rp)
end

function new_subreport(algorithm::Algorithm, rp::SubReport)
    new_subreport(algorithm, Logger(rp), rp)
end
