struct Initializer <: Algorithm end

function apply_hard_defaults(algorithm::Algorithm, defaults::Dict{Symbol, Any}, pre_logger::Function)
    pre_logger("TODO: apply_hard_defaults need to be implemented!!!", :warn)
end

function initialize(algorithm::Algorithm, problem::ScfProblem, guess_density::AbstractArray;
                    softdefaults::Defaults = Defaults(),
                    harddefaults::Defaults = Defaults(),
                    loglevel::LogLevel)

    # As long as we do not have a complete report yet we cannot log messages easily
    # use the pre_logger infrastructure to create a logging function
    logmessages = Vector{ReportMessage}()
    pre_logger = logger(logmessages, loglevel)
    pre_logger("Starting Initialization", :info)
    
    # Build initial iteration state
    pre_logger("Building initial iterate", :debug, :subalginit)
    fock, error_pulay, energies = compute_fock_matrix(problem, guess_density)
    pre_initstate = FockIterState(fock, error_pulay, energies, nothing, nothing, guess_density)

    # Run initial algorithm configuration
    pre_logger("initializing " * String(repr("text/plain", typeof(algorithm))), :debug, :subalginit)
    initstate = initialize_if_neccessary(algorithm, problem, pre_initstate, softdefaults)

    # Hard defaults override the values produced by automatic algorithm
    # configuration
    apply_hard_defaults(algorithm, harddefaults, pre_logger)

    # Construct report and add already existing log messages.
    pre_logger("setting up initial report", :debug,)
    report = Report(problem, initstate, algorithm, Vector{SubReport}(), loglevel)
    subreport = SubReport(Initializer(), problem, initstate, logmessages, nothing, loglevel, report)
    push!(report.history, subreport)

    # log a fancy header
    log!(subreport, @sprintf("%5s %14s %14s %14s %15s %12s\n", "iter", "e1e", "e2e", "etot", "scf_error", "n_applies"), :info)

    # return new report
    report
end

