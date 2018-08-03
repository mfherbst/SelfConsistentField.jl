#!/usr/bin/env python3
import molsturm
import molsturm.integrals
import h5py
import numpy as np

system = molsturm.System(
    atoms=["O", "H", "H"],
    coords=[[0, 0, 0],
            [0, 0, 1.795239827225189],
            [1.693194615993441, 0, -0.599043184453037]]
)
#system = molsturm.System("Be")
#system = molsturm.System("Ar")
params = molsturm.ScfParameters()
params.system = system
params.basis = molsturm.construct_basis("gaussian", params.system,
                                        basis_set_name="cc-pvdz")
integrals = ["electron_repulsion_bbbb", "nuclear_attraction_bb",
             "kinetic_bb", "overlap_bb"]

with h5py.File("integrals_h2o_dz.hdf5", "w") as h5f:
    for key in integrals:
        h5f.create_dataset(
            key, data=getattr(molsturm.integrals, key)(params),
            compression="gzip"
        )

    sys = h5f.create_group("system")
    sys.create_dataset("nelec", data=np.array((system.n_alpha, system.n_beta)))
    sys.create_dataset("atom_numbers", data=system.atom_numbers)
    sys.create_dataset("coords", data=system.coords)
