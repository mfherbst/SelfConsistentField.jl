module SelfConsistentField

using LinearAlgebra
using TensorOperations
using DataStructures
using Optim
using Printf

# Types
export ScfProblem, TwoElectronBuilder, JKBuilderFromTensor

# Algorithms
export ScfPipeline, ConditionalExec, ConvergenceCheck, Barrier, Roothaan, cDIIS, EDIIS, FixedDamping, initialize

# Functions
export run_scf, compute_nuclear_repulsion, assemble_hf_problem
export compute_pulay_error, is_closed_shell
export after_errnorm, before_errnorm, between_errnorm

# Utilities, which help in setting up an SCF problem or storing
# the result
export compute_guess, dump_molsturm_hdf5, load_integral_hdf5
export break_spin_symmetry

include("types/ScfIterState.jl")
include("types/ScfProblem.jl")
include("types/ScfConvergence.jl")
include("types/Report.jl")
include("types/DiisState.jl")

include("parts/compute_density.jl")
include("parts/compute_orbitals.jl")
include("parts/compute_pulay_error.jl")
include("parts/compute_fock_matrix.jl")
include("parts/check_convergence.jl")
include("parts/assemble_hf_problem.jl")
include("parts/misc.jl")
include("parts/algorithm.jl")
include("parts/reportmessage.jl")
include("parts/initialize.jl")
include("parts/new_subreport.jl")
include("parts/compute_error_norm.jl")

include("algorithms/JKBuilderFromTensor.jl")
include("algorithms/Roothaan.jl")
include("algorithms/FixedDamping.jl")
include("algorithms/cDIIS.jl")
include("algorithms/EDIIS.jl")
include("algorithms/ScfPipeline.jl")
include("algorithms/ConditionalExec.jl")
include("algorithms/ConvergenceCheck.jl")
include("algorithms/Barrier.jl")

include("util/guess.jl")
include("util/dump_molsturm_hdf5.jl")
include("util/load_integral_hdf5.jl")
include("util/break_spin_symmetry.jl")

include("parts/diis.jl")

include("run_scf.jl")

end # module
