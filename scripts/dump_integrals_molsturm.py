#!/usr/bin/env python3
import molsturm
import molsturm.integrals
import h5py
import numpy as np


def write_integrals_to_hdf5(params, h5f):
    system = params.system

    integrals = ["electron_repulsion_bbbb", "nuclear_attraction_bb",
                 "kinetic_bb", "overlap_bb"]
    for key in integrals:
        h5f.create_dataset(
            key, data=getattr(molsturm.integrals, key)(params),
            compression="gzip"
        )

    sys = h5f.create_group("system")
    sys.create_dataset("nelec", data=np.array((system.n_alpha, system.n_beta)))
    sys.create_dataset("atom_numbers", data=system.atom_numbers)
    sys.create_dataset("coords", data=system.coords)


def dump_integrals_gaussian(atoms, coords, electrons,
                            basis_set_name, ifile=None):
    system = molsturm.System(atoms, coords, electrons)

    params = molsturm.ScfParameters()
    params.system = system
    params.basis = molsturm.construct_basis("gaussian", params.system,
                                            basis_set_name=basis_set_name)

    if ifile is None:
        latoms = [a.lower() for a in atoms]
        ifile = "integrals_{}_{}.hdf5".format("".join(latoms), basis_set_name)

    with h5py.File(ifile, "w") as h5f:
        write_integrals_to_hdf5(params, h5f)

        g_discr = h5f.create_group("discretisation")
        g_discr.create_dataset("basis_type", data="gaussian",
                               dtype=h5py.special_dtype(vlen=str))
        g_discr.create_dataset("basis_set_name",
                               dtype=h5py.special_dtype(vlen=str),
                               data=basis_set_name)
        g_discr.create_dataset("has_real_harmonics", data=1, dtype=np.uint8)


def dump_integrals_sturmian(atoms, coords, electrons,
                            k_exp, n_max, l_max, m_max,
                            ifile=None):
    system = molsturm.System(atoms, coords, electrons)

    params = molsturm.ScfParameters()
    params.system = system
    basis = molsturm.construct_basis("sturmian/atomic", params.system,
                                     k_exp=k_exp, n_max=n_max,
                                     l_max=l_max, m_max=m_max)
    basis.backend = "cs_reference"
    params.basis = basis

    if ifile is None:
        latoms = [a.lower() for a in atoms]
        ifile = ("integrals_{}_{:02d}{:02d}{:02d}_{:.4f}.hdf5"
                 "".format("".join(latoms), n_max, l_max, m_max, k_exp))

    with h5py.File(ifile, "w") as h5f:
        write_integrals_to_hdf5(params, h5f)

        g_discr = h5f.create_group("discretisation")
        g_discr.create_dataset("basis_type", data="sturmian/atomic",
                               dtype=h5py.special_dtype(vlen=str))
        g_discr.create_dataset("nlm_basis", data=basis.functions)
        g_discr.create_dataset("k_exp", data=k_exp)
        g_discr.create_dataset("n_max", data=n_max)
        g_discr.create_dataset("l_max", data=l_max)
        g_discr.create_dataset("m_max", data=m_max)
        g_discr.create_dataset("has_real_harmonics",
                               data=basis.has_real_harmonics, dtype=np.uint8)


def main():
    atoms = ["O"]
    coords = [[0, 0, 0]]
    electrons = (5, 3)

    k_exp = 3.638
    n_max = 5
    l_max = 3
    m_max = 3

    dump_integrals_sturmian(atoms, coords, electrons,
                            k_exp, n_max, l_max, m_max)

    dump_integrals_gaussian(atoms, coords, electrons,
                            basis_set_name="3-21g")


if __name__ == "__main__":
    main()
