using TensorOperations
# using LinearAlgebra  v0.7.0

struct RestrictedClosedProblem <: ScfProblem
	# This struct defines the SCF problem to be solved

	# Nuclear repulsion energy
	energy_nuc_rep::Float64

	# Core Hamiltonian
	Hcore::AbstractArray

	# Electron repulsion integrals
	eri::AbstractArray

	# Overlap matrix
	overlap::AbstractArray

	# Number of occupied orbitals
	n_occ::Tuple{Int, Int}

	# Number of orbitals to be computed
	n_orbs::Int
end

"""
Update a Fock iterate using the provided orbital energies
and orbital coefficients
"""
function SelfConsistentField.compute_fock(
	problem::RestrictedClosedProblem, density::AbstractArray;
	kwargs...
)
	# TODO Split into compute fock and compute error

	# Aliases
	n_bas, _, n_spin = size(density)
	eri = problem.eri
	Hcore = problem.Hcore
	S = problem.overlap

	# Extract the alpha density
	assert(size(density) == (n_bas, n_bas, 1))
	Da = view(density, :,:,1)

	# Compute J + K matrix
	@tensor JKa[μ,ν] := eri[μ,ν,κ,λ] * 2 * Da[κ,λ] - eri[κ,ν,μ,λ] * Da[κ,λ]

	# Because we restrict ourselves to the alpha block, the acutal
	# energy is be twice the trace. For the 2e energy this cancels the factor 1/2
	energy_two_elec = trace(JKa * Da) # * 2 / 2
	energy_one_elec = 2 * trace(Hcore * Da)

	# Build alpha Fock matrix and error matrix
	focka = Hcore + JKa
	error = 2 * (S * Da * focka - focka * Da * S)
	# Factor 2 is used in the error, since both Fock matrix and density
	# are actually present identially in the beta block as well

	# Reshape both into expected shape
	error = reshape(error, (n_bas, n_bas, 1))
	fock = reshape(focka, (n_bas, n_bas, 1))

	total = energy_one_elec + energy_two_elec + problem.energy_nuc_rep
	energies = Dict("energy_1e"=> energy_one_elec,
		        "energy_2e"=> energy_two_elec,
			"energy_total"=> total)

	return fock, error, energies
end
