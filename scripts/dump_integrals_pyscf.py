#!/usr/bin/env python3
from collections import namedtuple
import pyscf
import pyscf.gto
import h5py
import numpy as np

mol = pyscf.gto.Mole()
mol.unit = "Bohr"
mol.atom = [
    ("O", [0, 0, 0]),
    ("H", [0, 0, 1.795239827225189]),
    ("H", [1.693194615993441, 0, -0.599043184453037]),
]
mol.basis = "cc-pvqz"
mol.spin = 0  # singlet
mol.build()
nao = mol.nao_nr()

Integral = namedtuple("Integral", ["name", "shape"])
int_map = {
    "int2e":      Integral("electron_repulsion_bbbb", 4 * (nao, )),
    "int1e_nuc":  Integral("nuclear_attraction_bb", 2 * (nao, )),
    "int1e_kin":  Integral("kinetic_bb", 2 * (nao, )),
    "int1e_ovlp": Integral("overlap_bb", 2 * (nao, )),
}

with h5py.File("integrals.hdf5", "w") as h5f:
    for key, integral in int_map.items():
        h5f.create_dataset(
            integral.name,
            data=mol.intor(key).reshape(integral.shape),
            compression="gzip"
        )

    sys = h5f.create_group("system")
    sys.create_dataset("nelec", data=np.array(mol.nelec))
    sys.create_dataset("atom_numbers", data=mol.atom_charges())
    sys.create_dataset("coords", data=np.array(mol.atom_coords()))
