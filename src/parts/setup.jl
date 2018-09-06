
struct UninitialisedAlgorithm
    algorithmtype::DataType
    args::Tuple
    params::Dict{Symbol, Any}
end

"""
    Create UninitializedAlgorithm given an Algorithm without full initalisation
    information. This is the case if no problem and iterate was given to the
    constructor. UninitialisedAlgorithm instances are promoted to Algorithm
    using the initalise function called by

    setup(::UninitialisedAlgorithm, ::ScfProblem, ::Iterate)
"""
function (new_algorithm_type::Type{T})(args...; params...) where {T<:Algorithm}
    function conv_to_alg(i)
        args[i] isa Tuple{Union{UninitialisedAlgorithm,Algorithm}, Function} ?
            ConditionalExec(args[i][1], args[i][2]; params...) : args[i]
    end

    # expand Tuple notation in arguments for conditional execution
    newargs = ntuple(conv_to_alg, length(args))
    args == newargs ? UninitialisedAlgorithm(new_algorithm_type, args, params) :
                      new_algorithm_type(newargs...; params...)
end

function initialise(uninit::UninitialisedAlgorithm, problem::ScfProblem, inititerate::Iterate, lg::Logger; global_params...)
    function init_uninitialised(i)
        !(uninit.args[i] isa UninitialisedAlgorithm) ? uninit.args[i] : begin
                initlg = Logger(lg)
                algorithm = initialise(uninit.args[i], problem, inititerate, initlg; global_params..., uninit.args[i].params...)
                log!(lg, "Logger", initlg, :debug, :initreport)
                algorithm
            end
    end
    initialised_args = ntuple(init_uninitialised, length(uninit.args))
    log!(lg, "algorithm type to be expanded", uninit.algorithmtype, :debug)
    uninit.algorithmtype(problem, inititerate, lg, initialised_args...; global_params..., uninit.params...)
end

struct Setup <: Algorithm end

function setup(uninit::UninitialisedAlgorithm, problem::ScfProblem, inititerate::Iterate;
                    params::Parameters = Parameters(), loglevel::LogLevel)

    lg = Logger(loglevel)
    log!(lg, "Starting Initialization", :debug)
    
    # Run initial algorithm configuration
    initlg = Logger(loglevel)
    algorithm = initialise(uninit, problem, inititerate, initlg; params...)
    if algorithm isa UninitialisedAlgorithm
        error("Some algorithm could not be initialised using its constructor.")
    end
    log!(lg, "Logger", initlg, :debug, :initreport)

    # Construct report and add already existing log messages.
    log!(lg, "setting up initial report", :debug, :firstreportsetup)
    report = ScfIterState(problem, inititerate, missing, algorithm, Vector{StepState}(), loglevel)

    # log a fancy header
    log!(lg, @sprintf("%5s %14s %14s %14s %15s %12s", "iter", "e1e", "e2e", "etot", "scf_error", "n_applies"), :info)

    stepstate = StepState(Setup(), problem, inititerate, missing, lg.messages, nothing, loglevel)
    push!(report.history, stepstate)

    # return new report
    report
end


"""
    Builds initial iteration iterate using given guess density.
"""
function build_initial_iterate(problem::ScfProblem, guess_density::AbstractArray)
    fock, error_pulay, energies = compute_fock_matrix(problem, guess_density)
    FockIterState(fock, error_pulay, energies, nothing, nothing, guess_density)
end

"""
    Alternative setup function using a guess density instead of an ScfProblem
"""
function setup(uninit::UninitialisedAlgorithm, problem::ScfProblem, guess_density::AbstractArray;
                    params::Parameters = Parameters(), loglevel::LogLevel)
    setup(uninit, problem, build_initial_iterate(problem, guess_density), params = params, loglevel = loglevel)
end
