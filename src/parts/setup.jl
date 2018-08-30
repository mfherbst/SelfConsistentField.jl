
struct UninitialisedAlgorithm
    algorithmtype::DataType
    args::Tuple
    params::Dict{Symbol, Any}
end

"""
    Create UninitializedAlgorithm given an Algorithm without full initalisation
    information. This is the case if no problem and state was given to the
    constructor. UninitialisedAlgorithm instances are promoted to Algorithm
    using the initalise function called by

    setup(::UninitialisedAlgorithm, ::ScfProblem, ::ScfIterState)
"""
function (new_algorithm_type::Type{T})(args...; params...) where {T<:Algorithm}
    function conv_to_alg(i)
        args[i] isa Tuple{Union{UninitialisedAlgorithm,Algorithm}, Function} ?
            ConditionalExec(args[i][1], args[i][2]; params...) : args[i]
    end

    # expand Tuple notation in arguments for conditional execution
    newargs = ntuple(conv_to_alg, length(args))
    #args == newargs && throw(MethodError(new_algorithm_type, args))
    args == newargs ? UninitialisedAlgorithm(new_algorithm_type, args, params) :
                      new_algorithm_type(args...; params...)
end

function initialise(uninit::UninitialisedAlgorithm, problem::ScfProblem, initstate::ScfIterState, lg::Logger; global_params...)
    for arg in uninit.args
        if arg isa UninitialisedAlgorithm
            initlg = Logger(rp)
            algorithm = initialise(uninit, problem, initstate, initlg; global_params..., params...)
            log!(lg, "Logger", initlg, :debug, :initreport)
        end
    end
    # Not all algorithms have to be initialized
    if applicable(uninit.algorithmtype, problem, initstate, lg, uninit.args...)
        uninit.algorithmtype(problem, initstate, lg, uninit.args...; global_params..., uninit.params...)
    end
end

struct Setup <: Algorithm end

function setup(uninit::UninitialisedAlgorithm, problem::ScfProblem, initstate::ScfIterState;
                    params::Parameters = Parameters(), loglevel::LogLevel)

    lg = Logger(loglevel)
    log!(lg, "Starting Initialization", :debug)
    
    # Run initial algorithm configuration
    initlg = Logger(loglevel)
    algorithm = initialise(uninit, problem, initstate, initlg; params...)
    log!(lg, "Logger", initlg, :debug, :initreport)

    # Construct report and add already existing log messages.
    log!(lg, "setting up initial report", :debug, :firstreportsetup)
    report = Report(problem, initstate, missing, algorithm, Vector{SubReport}(), loglevel)

    # log a fancy header
    log!(lg, @sprintf("%5s %14s %14s %14s %15s %12s", "iter", "e1e", "e2e", "etot", "scf_error", "n_applies"), :info)

    subreport = SubReport(Setup(), problem, initstate, missing, lg.messages, nothing, loglevel)
    push!(report.history, subreport)

    # return new report
    report
end


"""
    Builds initial iteration state using given guess density.
"""
function build_initial_state(problem::ScfProblem, guess_density::AbstractArray)
    fock, error_pulay, energies = compute_fock_matrix(problem, guess_density)
    FockIterState(fock, error_pulay, energies, nothing, nothing, guess_density)
end

"""
    Alternative setup function using a guess density instead of an ScfProblem
"""
function setup(uninit::UninitialisedAlgorithm, problem::ScfProblem, guess_density::AbstractArray;
                    params::Parameters = Parameters(), loglevel::LogLevel)
    setup(uninit, problem, build_initial_state(problem, guess_density), params = params, loglevel = loglevel)
end
