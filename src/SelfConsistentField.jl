module SelfConsistentField

using LinearAlgebra: norm, Hermitian, eigen
using TensorOperations

export run_scf, ScfProblem, compute_guess_hcore

include("types/ScfIterState.jl")
include("types/ScfProblem.jl")
include("types/Accelerator.jl")

include("auxiliary/compute_density.jl")
include("auxiliary/compute_orbitals.jl")
include("auxiliary/check_convergence.jl")

include("guess/hcore.jl")

include("algorithms/Roothaan.jl")
include("algorithms/FixedDamping.jl")

include("run_scf.jl")

end # module
