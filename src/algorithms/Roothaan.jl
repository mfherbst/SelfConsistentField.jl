struct Roothaan <: Algorithm end

"""
Perform one step of the Roothaan algorithm if the Fock matrix is iterated
"""
function Base.iterate(::Roothaan, rp::SubReport)
    orben, orbcoeff = compute_orbitals(rp.problem, rp.state.fock)
    density = compute_density(rp.problem, orbcoeff)
    fock, error_pulay, energies = compute_fock_matrix(rp.problem, density)
    return FockIterState(fock, error_pulay, energies, orbcoeff, orben, density)
end
