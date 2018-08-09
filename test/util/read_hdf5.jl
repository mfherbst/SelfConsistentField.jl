using HDF5

"""
Read an hdf5 file
"""
function read_hdf5(file)
    file = h5open(file)

    eri = read(file, "electron_repulsion_bbbb")
    T = read(file,"kinetic_bb")
    V = read(file,"nuclear_attraction_bb")
    S = read(file,"overlap_bb")
    JK = JKBuilderFromTensor(eri)

    nelec = convert(Tuple{Int64, Int64}, (read(file,"system/nelec")...,))
    atnums = read(file, "system/atom_numbers")
    coords = read(file, "system/coords")

    close(file)

    # coordinates need to be transformed from row-major
    # to column-major
    coords = Array{Float64}(adjoint(coords))

    # Build named tuples and return
    system = (coords=coords, atom_numbers=atnums, n_elec=nelec)
    integrals = (kinetic_bb=T, nuclear_attraction_bb=V,
                 overlap_bb=S, electron_repulsion_bbbb=eri,
                 jk_builder_bb=JK)
    return system, integrals
end

