"""
Break the spin symmetry in the α and the β block of the
passed density.
"""
function break_spin_symmetry(density)
    n_bas, _, n_spin = size(density)
    if n_spin == 1
        return density  # Only single spin component,
        #                 no breaking can be done
    end

    # Ideally one would only throw away the density matrix
    # elements which connect different atoms, but this is
    # a little sophisticated for us for now. Instead we will
    # make the β density tridiagonal
    density[:,:,2] = Tridiagonal(view(density, :, :, 2))
    return density
end
