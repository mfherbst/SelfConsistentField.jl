using SelfConsistentField
include("UnrestrictedProblem.jl")
include("../common/Setup.jl")

if length(ARGS) != 1
	println("Provide integral file on first arg")
	exit(1)
end
intfile = ARGS[1]
system, integrals = read_hdf5(intfile)

println()
println("Running unrestricted SCF for $intfile")
println()

n_bas = size(integrals.electron_repulsion_bbbb)[1]
n_occ_a = system.nelecs[1]
n_occ_b = system.nelecs[2]
n_occ = (n_occ_a, n_occ_b)
println("Number of electrons:  $n_occ_a alpha and $n_occ_b beta")
println("Nbas: $n_bas    nocca: $n_occ_a    noccb: $n_occ_b")

ene_nuc_rep = compute_nuclear_repulsion(system)
println("Nuclear repulsion energy: $ene_nuc_rep")

problem = UnrestrictedProblem(
	ene_nuc_rep,
	integrals.kinetic_bb + integrals.nuclear_attraction_bb,
	integrals.electron_repulsion_bbbb,
	integrals.overlap_bb, n_occ, n_bas)

# Compute HCore guess density and duplicate it on the spin components
guess_density = compute_guess_hcore(problem, system.coords, system.atom_numbers)
guess_density = cat(3, guess_density, guess_density)

# Break the symmetry in alpha and beta densities by adding a random diagonal
guess_diagonal = 1e-3 * randn(n_bas)
guess_density[:, :, 1] += diagm(guess_diagonal)
guess_density[:, :, 2] -= diagm(guess_diagonal)

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
