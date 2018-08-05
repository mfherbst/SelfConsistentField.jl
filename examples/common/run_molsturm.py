#!/usr/bin/env python3
import molsturm
import h5py
import sys
import numpy as np

if len(sys.argv) != 3:
    raise SystemExit("Need hdf5 integrals file as first arg and basis set as second")
ifile = sys.argv[1]
bas = sys.argv[2]

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
params.basis = molsturm.construct_basis("gaussian", params.system,
                                        basis_set_name=bas)
params["scf/print_iterations"] = True

res = molsturm.self_consistent_field(params)
molsturm.print_convergence_summary(res)
molsturm.print_energies(res)
molsturm.print_mo_occupation(res)
