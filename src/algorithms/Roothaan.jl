struct Roothaan <: Algorithm end

"""
Perform one step of the Roothaan algorithm if the Fock matrix is iterated
"""
function iterate(roothan::Roothaan, rp::SubReport)
    newrp = new_subreport(rp)
    orben, orbcoeff = compute_orbitals(rp.problem, rp.state.fock)
    density = compute_density(rp.problem, orbcoeff)
    fock, error_pulay, energies = compute_fock_matrix(rp.problem, density)
    newrp.state = FockIterState(fock, error_pulay, energies, orbcoeff, orben, density)
    return roothan, newrp
end
