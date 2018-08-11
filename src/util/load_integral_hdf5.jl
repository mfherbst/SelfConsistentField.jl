using HDF5

"""
Read an hdf5 file of integral data
"""
function load_integral_hdf5(hdf5_file)
    file = h5open(hdf5_file)

    eri = read(file, "electron_repulsion_bbbb")
    T = read(file,"kinetic_bb")
    V = read(file,"nuclear_attraction_bb")
    S = read(file,"overlap_bb")
    JK = JKBuilderFromTensor(eri)

    nelec = convert(Tuple{Int64, Int64}, (read(file,"system/nelec")...,))
    atnums = read(file, "system/atom_numbers")
    coords = read(file, "system/coords")

    discretisation=nothing
    basis_type = read(file, "discretisation/basis_type")
    if startswith(basis_type, "sturmian")
        discretisation = (
            basis_type=read(file, "discretisation/basis_type"),
            nlm_basis=read(file, "discretisation/nlm_basis"),
            k_exp=read(file, "discretisation/k_exp"),
            n_max=read(file, "discretisation/n_max"),
            l_max=read(file, "discretisation/l_max"),
            m_max=read(file, "discretisation/m_max"),
        )
    elseif basis_type == "gaussian"
        discretisation = (
            basis_type=read(file, "discretisation/basis_type"),
            basis_set_name=read(file, "discretisation/basis_set_name"),
        )
    end
    close(file)

    # coordinates need to be transformed from row-major
    # to column-major
    coords = Array{Float64}(adjoint(coords))

    # Build named tuples and return
    system = (coords=coords, atom_numbers=atnums, n_elec=nelec)
    integrals = (kinetic_bb=T, nuclear_attraction_bb=V,
                 overlap_bb=S, electron_repulsion_bbbb=eri,
                 jk_builder_bb=JK, discretisation=discretisation)
    return system, integrals
end

