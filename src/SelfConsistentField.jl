module SelfConsistentField

using LinearAlgebra
using TensorOperations
using DataStructures
using Optim

# Types
export ScfProblem, TwoElectronBuilder, JKBuilderFromTensor

# Functions
export run_scf, compute_nuclear_repulsion, assemble_hf_problem
export compute_pulay_error, is_closed_shell

# Utilities, which help in setting up an SCF problem or storing
# the result
export compute_guess, dump_molsturm_hdf5, load_integral_hdf5
export break_spin_symmetry

include("types/ScfIterState.jl")
include("types/ScfProblem.jl")
include("types/Accelerator.jl")
include("types/DiisState.jl")

include("parts/compute_density.jl")
include("parts/compute_orbitals.jl")
include("parts/compute_pulay_error.jl")
include("parts/compute_fock_matrix.jl")
include("parts/check_convergence.jl")
include("parts/assemble_hf_problem.jl")
include("parts/misc.jl")

include("algorithms/JKBuilderFromTensor.jl")
include("algorithms/Roothaan.jl")
include("algorithms/FixedDamping.jl")
include("algorithms/cDIIS.jl")
include("algorithms/EDIIS.jl")

include("util/guess.jl")
include("util/dump_molsturm_hdf5.jl")
include("util/load_integral_hdf5.jl")
include("util/break_spin_symmetry.jl")

include("parts/diis.jl")

include("run_scf.jl")

end # module
