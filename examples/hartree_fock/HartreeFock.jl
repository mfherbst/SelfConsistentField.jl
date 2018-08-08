using SelfConsistentField
include("../../test/util/read_hdf5.jl")
include("print_results.jl")

function parse_args()
    error = false
    restricted = nothing
    if length(ARGS) == 1
        intfile = ARGS[1]
    elseif length(ARGS) == 2
        if ARGS[1] == "--unrestricted"
            restricted = false
        elseif ARGS[1] == "--restricted"
            restricted = true
        else
            error = true
        end
        intfile = ARGS[2]
    else
        error = true
    end
    if error
        println("$PROGRAM_FILE [--unrestricted|--restricted] <integral_file>")
        exit(1)
    end
    return restricted, intfile
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

function main()
    restricted, intfile = parse_args()

    println("Reading file $intfile")
    system, integrals = read_hdf5(intfile)
    println()

    problem = assemble_hf_problem(system, integrals; restricted=restricted)
    print_info(problem)
    println()

    # Compute HCore guess
    guess_density = compute_guess_hcore(problem, system.coords, system.atom_numbers)

    params = Dict(
        :max_error_norm=>5e-7,
        :max_energy_total_change=>1.25e-07,
    )
    res = run_scf(problem, guess_density; params...)
    print_results(problem, integrals, res)
end

main()
