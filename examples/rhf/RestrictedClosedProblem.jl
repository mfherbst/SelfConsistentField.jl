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
	n_occ::Int

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
	eri = problem.eri
	Hcore = problem.Hcore
	S = problem.overlap
	D = density

	@tensoropt JK[μ,ν] := 2 * eri[α,β,μ,ν] * D[α,β] - eri[μ,β,α,ν] * D[α,β]
	energy_two_elec = trace(JK * D)
	energy_one_elec = 2 * trace(Hcore * D)

	fock = Hcore + JK
	error = S * D * fock - fock * D * S

	total = energy_one_elec + energy_two_elec + problem.energy_nuc_rep
	energies = Dict("energy_1e"=> energy_one_elec,
		        "energy_2e"=> energy_two_elec,
			"energy_total"=> total)

	return fock, error, energies
end

function SelfConsistentField.get_n_spin(problem::RestrictedClosedProblem)
	return 1
end
