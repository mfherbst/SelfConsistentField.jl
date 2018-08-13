#!/usr/bin/env python3
import pyscf
import h5py
import sys
import numpy as np

if len(sys.argv) != 3:
    raise SystemExit("Need args: hdf5 integrals file, restricted")
ifile = sys.argv[1]

if sys.argv[2] == "true":
    restricted = True
elif sys.argv[2] == "false":
    restricted = False
else:
    raise SystemExit("restricted == true|false")

with h5py.File(ifile, "r") as h5f:
    nelec = list(h5f["system/nelec"])
    atom_numbers = np.array(h5f["system/atom_numbers"])
    coords = np.array(h5f["system/coords"])
    basis_set_name = str(h5f["discretisation/basis_set_name"].value)

mol = pyscf.gto.Mole()
mol.unit = "Bohr"
mol.atom = [(str(atom_numbers[i]), coords[i]) for i in range(len(atom_numbers))]
mol.nelec = nelec
mol.basis = basis_set_name
mol.verbose = 4
mol.build()

if mol.nelectron == 1 or mol.spin == 0:
    mf = pyscf.scf.RHF(mol)
    if not restricted:
        raise SystemExit("Closed-shell systems are always restricded")
else:
    if restricted:
        mf = pyscf.scf.ROHF(mol)
    else:
        mf = pyscf.scf.UHF(mol)
mf.kernel()
