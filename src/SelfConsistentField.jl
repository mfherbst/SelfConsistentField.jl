module SelfConsistentField

using LinearAlgebra
using TensorOperations
using DataStructures
using Optim
using Printf

# Types
export ScfProblem, TwoElectronBuilder, JKBuilderFromTensor

# Algorithms
export ScfPipeline, ConditionalExec, ConvergenceCheck,
       Barrier, ComputeOrbitals, ComputeDensity, ComputeFock,
       cDIIS, EDIIS, FixedDamping, setup

# Functions
export compute_nuclear_repulsion, assemble_hf_problem
export compute_pulay_error, is_closed_shell
export after_errnorm, before_errnorm, between_errnorm

# Utilities, which help in setting up an SCF problem or storing
# the result
export compute_guess, dump_molsturm_hdf5, load_integral_hdf5
export break_spin_symmetry

include("types/ScfIterState.jl")
include("types/ScfProblem.jl")
include("types/ScfConvergence.jl")
include("types/ReportMessage.jl")
include("types/Report.jl")
include("types/DiisState.jl")

include("parts/compute_pulay_error.jl")
include("parts/check_convergence.jl")
include("parts/assemble_hf_problem.jl")
include("parts/misc.jl")
include("parts/algorithm.jl")
include("parts/report.jl")
include("parts/ranged_applications.jl")
include("parts/reportmessage.jl")
include("parts/setup.jl")
include("parts/new_subreport.jl")
include("parts/compute_error_norm.jl")

include("algorithms/ComputeDensity.jl")
include("algorithms/ComputeOrbitals.jl")
include("algorithms/ComputeFock.jl")
include("algorithms/JKBuilderFromTensor.jl")
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

end # module
