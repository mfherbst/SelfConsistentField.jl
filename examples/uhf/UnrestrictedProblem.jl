using TensorOperations
using LinearAlgebra: tr

struct UnrestrictedProblem <: ScfProblem
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
    n_occ::Tuple{Int,Int}

    # Number of orbitals to be computed
    n_orbs::Int
end

"""
Update a Fock iterate using the provided orbital energies
and orbital coefficients
"""
function SelfConsistentField.compute_fock(
    problem::UnrestrictedProblem, density::AbstractArray;
    kwargs...
)
    # TODO Split into compute fock and compute error

    # Aliases
    n_bas, _, n_spin = size(density)
    eri = problem.eri
    Hcore = problem.Hcore
    S = problem.overlap

    # Extract the alpha and beta density
    @assert size(density) == (n_bas, n_bas, 2)
    Da = view(density, :,:,1)
    Db = view(density, :,:,2)

    @assert size(Da) == (n_bas, n_bas)
    @assert size(Db) == (n_bas, n_bas)
    @assert size(eri) == (n_bas, n_bas, n_bas, n_bas)

    # Compute J, Ka and Kb matrices
    @tensor begin
        J[μ,ν] := eri[α,β,μ,ν] * Da[α,β] + eri[α,β,μ,ν] * Db[α,β]
        JKa[μ,ν] := J[μ,ν] - eri[μ,β,α,ν] * Da[α,β]
        JKb[μ,ν] := J[μ,ν] - eri[μ,β,α,ν] * Db[α,β]
    end
    J = 0  # Free memory in J

    energy_two_elec = 1/2 * (tr(JKa * Da) + tr(JKb * Db))
    energy_one_elec = tr(Hcore * Da) + tr(Hcore * Db)

    # Build alpha fock matrix and error matrix
    fock = Array{eltype(density)}(n_bas, n_bas, n_spin)
    fock[:, :, 1] = Hcore + JKa
    fock[:, :, 2] = Hcore + JKb

    error = Array{eltype(density)}(n_bas, n_bas, n_spin)
    for s in 1:n_spin
        error[:, :, s] = (S * view(density, :, :, s) * view(fock, :, :, s)
                          - view(fock, :, :, s) * view(density, :, :, s) * S)
    end

    total = energy_one_elec + energy_two_elec + problem.energy_nuc_rep
    energies = Dict("1e"=> energy_one_elec,
                    "2e"=> energy_two_elec,
                    "total"=> total)

    return fock, error, energies
end
