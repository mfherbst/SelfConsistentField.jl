"""
Datatype for an SCF problem, regardless wether restricted or unrestricted
"""
abstract type ScfProblem end

function get_n_spin(problem::ScfProblem)
	throw("Overload get_n_spin for your ScfProblem")
end

function compute_fock(problem::ScfProblem, density::AbstractArray; kwargs...)
	throw("Overload compute_fock for your ScfProblem")
end
