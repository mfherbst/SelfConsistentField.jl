using TensorOperations
using LinearAlgebra: tr, I

struct RestrictedOpenProblem <: ScfProblem
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
Function to build Roothaan's effective Fock  matrix, reference
http://www-theor.ch.cam.ac.uk/people/ross/thesis/node15.html
and pyscf ROHF procedure

fock_ab      Tuple of alpha, beta Fock matrix
density_ab   Tuple of alpha, beta density matrix
overlap      Overlap matrix
"""
function build_roothaan_effective_fock(fock_ab, density_ab, overlap)
    Fa, Fb = fock_ab
    Da, Db = density_ab
    S = overlap
    n_bas, _ = size(overlap)

    # The effective ROHF Fock matrix in MOs is constructed as
    #
    #          |  closed   open    virtual
    #  ----------------------------------------
    #  closed  |    Fc      Fb       Fc
    #  open    |    Fb      Fc       Fa
    #  virtual |    Fc      Fa       Fc
    #  ----------------------------------------
    #
    # where
    Fc = (Fa + Fb) / 2

    # Projectors for closed-shell, open-shell and virtual spaces
    Pc = Db * S
    Po = (Da - Db) * S
    Pv = I - Da * S

    # Now build the effective Fock matrix. The idea is to
    # build the lower triangle only and than use the hermitian
    # conjugate to build the full thing. Implies that a factor
    # 1/2 is needed along the diagonal

    Feff = (
          Pc' * Fc * Pc / 2
        + Po' * Fb * Pc     + Po' * Fc * Po / 2
        + Pv' * Fc * Pc     + Pv' * Fa * Po     + Pv' * Fc * Pv / 2
    )
    return Feff + Feff'
end


"""
Update a Fock iterate using the provided orbital energies
and orbital coefficients
"""
function SelfConsistentField.compute_fock(
   problem::RestrictedOpenProblem, density::AbstractArray;
   kwargs...
)
    # TODO Split into compute fock and compute error
    # The compute error function should check the number of spin
    # components in the passed density matrix. For ROHF one would
    # just pass Dtot, whereas for RHF just Da and for UHF (Da, Db) as a tuple.

    # Aliases
    n_bas, _, _ = size(density)
    n_docc = min(problem.n_occ...)
    n_socc = max(problem.n_occ...)
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

    # Build alpha and beta fock matrix and transform to MO
    focka = Hcore + JKa
    fockb = Hcore + JKb

    #
    # Up to here exactly as in UHF

    # Build Fock and error matrix (using the total density)
    fock = build_roothaan_effective_fock((focka, fockb), (Da, Db), S)

    Dtot = Da + Db
    error = S * Dtot * fock - fock * Dtot * S

    # Resize fock matrix to 1 spin component as required
    fock = reshape(fock, n_bas, n_bas, 1)

    # From here it is exactly as in RHF or UHF
    #
    total = energy_one_elec + energy_two_elec + problem.energy_nuc_rep
    energies = Dict("1e"=> energy_one_elec,
                    "2e"=> energy_two_elec,
                    "total"=> total)

    return fock, error, energies
end
