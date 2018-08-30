#!/usr/bin/env julia

using TensorOperations
using SelfConsistentField
include("print_results.jl")

function parse_args()
    have_error = false
    restricted = nothing
    ofile = nothing

    next = 1
    for i in next:length(ARGS)
        next = i
        if ARGS[i] == "--unrestricted"
            restricted = false
        elseif ARGS[i] == "--restricted"
            restricted = true
        elseif ARGS[i] in ["-h", "--help"]
            println("$PROGRAM_FILE [--unrestricted|--restricted] "
                    * "<integral_file> [<outfile>]")
        else
            break
        end
    end

    if next > length(ARGS)
        error("No input hdf5 integral file given")
    end

    # Now parse integral file
    if ! isfile(ARGS[next])
        error("Not a valid hdf5 integral file: $(ARGS[next])")
    end
    intfile = ARGS[next]
    next += 1

    # If it exists, parse output file
    if next <= length(ARGS)
        ofile = ARGS[next]
        next += 1
    end

    if next <= length(ARGS)
        error("Excess commandline arguments: $(ARGS[next:end])")
    end
    return restricted, intfile, ofile
end

function print_info(problem)
    scftype = "ROHF"
    if problem.restricted
        if is_closed_shell(problem)
            scftype = "RHF"
        end
    else
        scftype = "UHF"
    end

    println("Problem information:")
    println("   № α electrons:     ", problem.n_elec[1])
    println("   № β electrons:     ", problem.n_elec[2])
    println("   № basis functions: ", size(problem.h_core)[1])
    println("   SCF type:          $scftype")
end

function compute_termwise_hf_energies!(res, problem)
    n_spin = size(res["density"], 3)
    energies = res["energies"]
    @assert n_spin == 1 || n_spin == 2

    # Get view to α and β density
    Dα = view(res["density"], :, :, 1)
    Dβ = view(res["density"], :, :, n_spin)

    for (label, term) in problem.terms_0e
        energies[label] = term
    end
    for (label, term) in problem.terms_1e
        energies[label] = tr(Dα * term) + tr(Dβ * term)
    end

    @assert length(problem.terms_2e) == 1
    JKBuilder = problem.terms_2e["coulomb+exchange"]
    @assert isa(JKBuilder, JKBuilderFromTensor)

    eri = JKBuilder.eri
    @tensor begin
        energy_J = 1/2 * scalar(Dα[μ,ν] * eri[α,β,μ,ν] * Dα[α,β] +
                                Dα[μ,ν] * eri[α,β,μ,ν] * Dβ[α,β] +
                                Dβ[μ,ν] * eri[α,β,μ,ν] * Dα[α,β] +
                                Dβ[μ,ν] * eri[α,β,μ,ν] * Dβ[α,β])
        energy_K = -1/2 * scalar(Dα[μ,ν] * eri[μ,β,α,ν] * Dα[α,β] +
                                 Dβ[μ,ν] * eri[μ,β,α,ν] * Dβ[α,β])
    end
    energies["coulomb"] = energy_J
    energies["exchange"] = energy_K
end

function hartree_fock(intfile; restricted=nothing, ofile=nothing)
    println("Reading file $intfile")
    system, integrals = load_integral_hdf5(intfile)
    println()

    problem = assemble_hf_problem(system, integrals; restricted=restricted)
    print_info(problem)
    println()

    #guess_method = "hcore"
    guess_method = "hcore"
    if startswith(integrals.discretisation.basis_type, "sturmian")
        guess_method = "random"
    end

    harddefaults = Dict(
        :max_error_norm=>5e-7,
        :max_energy_total_change=>1.25e-07,
    )

    guess_density = compute_guess(problem, system.coords, system.atom_numbers;
                                  method = guess_method)
    if ! problem.restricted && is_closed_shell(problem)
        break_spin_symmetry(guess_density)
    end

    #ecdiis = FallbackMechanism(EDIIS(), cDIIS(); n_fallback_iterations = 5, fallback_predicate = nrtuff)
    roothan = ScfPipeline(
                          ComputeOrbitals(),
                          ComputeDensity(),
                          ComputeFock()
                         )

    algorithm = ScfPipeline(
        roothan,
        cDIIS(),
    #    (
    #     Barrier(
    #             EDIIS(),
    #             (cDIIS(), before_errnorm(10e-5)),
    #             after_errnorm(10e-2)
    #            ),
    #     before_errnorm(10e-7)
    #    ),
    #    (
    #     FixedDamping(),
    #     between_errnorm(10e-2, 10e-7)
    #    ),
        ConvergenceCheck(; max_error_norm = 10e-10)
    )

    #acceleration = ediis, cediis 
    #
    #ScfAlgorithm(
    #    UpdateCoefficients(),
    #    UpdateFock(),
    #    (acceleration, should_apply_accerator),
    #    damping,
    #    ConvergenceCheck(is_converged)
    #)

    loglevel = Dict{Symbol, Set}(:stdout => Set([:info, :warn]))
    solver = setup(algorithm, problem, guess_density; :loglevel => loglevel)
    collect(solver)
    res = convert(Dict, solver)

    #res = run_scf(problem, guess_density; params...)

    if ! res["converged"]
        error("SCF failed to converge")
    end
    compute_termwise_hf_energies!(res, problem)

    print_results(problem, res)
    if ofile != nothing
        dump_molsturm_hdf5(integrals, res, ofile)
    end
end

function main()
    restricted, intfile, ofile = parse_args()
    hartree_fock(intfile; restricted=restricted, ofile=ofile)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
