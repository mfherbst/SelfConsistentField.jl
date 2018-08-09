#!/usr/bin/env python3
import molsturm
import h5py
import sys
import numpy as np

if len(sys.argv) != 2:
    raise SystemExit("Need hdf5 integrals file as first arg.")
ifile = sys.argv[1]

with h5py.File(ifile, "r") as h5f:
    nelec = list(h5f["system/nelec"])
    atom_numbers = np.array(h5f["system/atom_numbers"])
    coords = np.array(h5f["system/coords"])

    system = molsturm.System(
        atoms=atom_numbers,
        coords=coords,
        electrons=nelec,
    )
    params = molsturm.ScfParameters()
    params.system = system
    params["scf/print_iterations"] = True

    basis_type = str(h5f["discretisation/type"][0])
    if basis_type == "gaussian":
        bas = str(h5f["discretisation/basis_set_name"][0])
        params.basis = molsturm.construct_basis("gaussian", params.system,
                                                basis_set_name=bas)
    elif basis_type == "sturmian/atomic":
        k_exp = float(h5f["discretisation/k_exp"].value)
        n_max = int(h5f["discretisation/n_max"].value)
        l_max = int(h5f["discretisation/l_max"].value)
        m_max = int(h5f["discretisation/m_max"].value)
        params.basis = molsturm.construct_basis("sturmian/atomic", params.system,
                                                k_exp=k_exp, n_max=n_max, l_max=l_max,
                                                m_max=m_max)


res = molsturm.self_consistent_field(params)
molsturm.print_convergence_summary(res)
molsturm.print_energies(res)
molsturm.print_mo_occupation(res)
