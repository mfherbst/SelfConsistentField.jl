module SelfConsistentField

using LinearAlgebra
using TensorOperations
using DataStructures

export run_scf, ScfProblem, compute_guess_hcore
export assemble_hf_problem, compute_nuclear_repulsion
export compute_pulay_error, TwoElectronBuilder, is_closed_shell
export JKBuilderFromTensor

include("types/ScfIterState.jl")
include("types/ScfProblem.jl")
include("types/Accelerator.jl")

include("auxiliary/compute_density.jl")
include("auxiliary/compute_orbitals.jl")
include("auxiliary/compute_pulay_error.jl")
include("auxiliary/compute_fock_matrix.jl")
include("auxiliary/check_convergence.jl")
include("auxiliary/assemble_hf_problem.jl")
include("auxiliary/misc.jl")

include("guess/hcore.jl")

include("algorithms/JKBuilderFromTensor.jl")
include("algorithms/Roothaan.jl")
include("algorithms/FixedDamping.jl")
include("algorithms/cDIIS.jl")

include("run_scf.jl")

end # module
