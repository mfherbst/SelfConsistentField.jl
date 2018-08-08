using SelfConsistentField
include("JKBuilder.jl")
include("read_hdf5.jl")
include("print_results.jl")

if length(ARGS) != 2 || ! (ARGS[1] in ["--restricted", "--unrestricted"])
    println("$PROGRAM_FILE [--unrestricted|--restricted] <integral_file>")
    exit(1)
end

restricted = true
if ARGS[1] == "--unrestricted"
    restricted = false
end
intfile = ARGS[2]

println("Reading file $intfile")
system, integrals = read_hdf5(intfile)
println()

problem = assemble_hf_problem(system, integrals; restricted=restricted)

scftype = "ROHF"
if problem.restricted
    if is_closed_shell(problem)
        scftype = "RHF"
    end
else
    scftype = "UHF"
end

println("Problem information:")
println("   # α electrons:     ", system.n_elec[1])
println("   # β electrons:     ", system.n_elec[2])
println("   # basis functions: ", size(integrals.overlap_bb)[1])
println("   SCF type:          $scftype")
println()

# Compute HCore guess
guess_density = compute_guess_hcore(problem, system.coords, system.atom_numbers)

params = Dict(
    :max_error_norm=>5e-7,
    :max_energy_total_change=>1.25e-07,
)
res = run_scf(problem, guess_density; params...)

if res["is_converged"]
    println("SCF converged.")
    println()
    print_energies(problem, integrals, res)
    println()
    print_mo_occupation(problem, res)
end
