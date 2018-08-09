#!/usr/bin/env python3
from collections import namedtuple
import pyscf
import pyscf.gto
import h5py
import numpy as np


def dump_integrals_gaussian(atoms, coords, electrons,
                            basis_set_name, ifile=None):
    mol = pyscf.gto.Mole()
    mol.unit = "Bohr"
    mol.atom = mol.atom = [(str(atoms[i]), coords[i])
                           for i in range(len(atoms))]
    mol.nelec = electrons
    mol.basis = basis_set_name
    mol.build()
    nao = mol.nao_nr()

    Integral = namedtuple("Integral", ["name", "shape"])
    int_map = {
        "int2e":      Integral("electron_repulsion_bbbb", 4 * (nao, )),
        "int1e_nuc":  Integral("nuclear_attraction_bb", 2 * (nao, )),
        "int1e_kin":  Integral("kinetic_bb", 2 * (nao, )),
        "int1e_ovlp": Integral("overlap_bb", 2 * (nao, )),
    }

    if ifile is None:
        ifile = "integrals_{}_{}.hdf5".format("".join(atoms), basis_set_name)

    with h5py.File(ifile, "w") as h5f:
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

        basis = h5f.create_group("discretisation")
        basis.create_dataset("type",
                             data=np.array(["gaussian"],
                                           dtype=h5py.special_dtype(vlen=str)))
        basis.create_dataset("basis_set_name",
                             data=np.array([basis_set_name],
                                           dtype=h5py.special_dtype(vlen=str)))


def main():
    atoms = ["O", "H", "H"]
    coords = [[0, 0, 0], [0, 0, 1.795239827225189],
              [1.693194615993441, 0, -0.599043184453037]]
    electrons = (5, 5)
    basis_set_name = "cc-pvdz"

    dump_integrals_gaussian(atoms, coords, electrons, basis_set_name)


if __name__ == "__main__":
    main()
