"""
Datatype for an SCF problem, regardless wether restricted or unrestricted
"""
abstract type ScfProblem end

function compute_fock(problem::ScfProblem, density::AbstractArray; kwargs...)
	throw("Overload compute_fock for your ScfProblem")
end
