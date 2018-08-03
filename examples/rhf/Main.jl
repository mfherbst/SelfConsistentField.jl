using SelfConsistentField
include("RestrictedClosedProblem.jl")
include("../common/Setup.jl")

if length(ARGS) != 1
	println("Provide integral file on first arg")
	exit(1)
end
intfile = ARGS[1]
system, integrals = read_hdf5(intfile)

println()
println("Running restricted SCF for $intfile")
println()

nelec = system.nelecs[1] + system.nelecs[2]
println("Number of electrons:  $nelec")
nbas = size(integrals.electron_repulsion_bbbb)[1]
nocc = div(nelec, 2)
println("Nbas: $nbas    nocc: $nocc")

ene_nuc_rep = compute_nuclear_repulsion(system)
println("Nuclear repulsion energy: $ene_nuc_rep")

problem = RestrictedClosedProblem(
	ene_nuc_rep,
	integrals.kinetic_bb + integrals.nuclear_attraction_bb,
	integrals.electron_repulsion_bbbb,
	integrals.overlap_bb, nocc, nbas)
guess_density = compute_guess_hcore(problem, system.coords, system.atom_numbers)
params = Dict(
	:max_error_norm=>5e-7,
	:max_energy_total_change=>1.25e-07,
	:max_iter=>45,
)
res = run_scf(problem, guess_density; params...)

if res["is_converged"]
	println("SCF converged.")
	println()
	print_energies(problem, integrals, res)
	println()
	print_mo_occupation(problem, res)
end
