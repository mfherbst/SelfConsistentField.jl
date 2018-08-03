"""
Type describing the SCF quantity to be iterated by the SCF map
(e.g. Fock matrix)
"""
abstract type ScfIterState end

"""
Struct of the state if the Fock matrix is iterated
"""
struct FockIterState <: ScfIterState
	# Fock or Kohn-Sham matrix
	fock::AbstractArray

	# Matrix describing the error in fock
	error::AbstractArray

	# Named tuple or dict of the computed energy terms
	energies::Dict{String, Float64}

	# Orbital coefficients and energies
	orbcoeff::Union{Void,AbstractArray}
	orben::Union{Void,AbstractArray}
end

function get_iterate_matrix(iterate::FockIterState)
	return iterate.fock
end

function update_iterate_matrix(iterate::FockIterState, fock::AbstractArray)
	return FockIterState(fock, iterate.error, iterate.energies,
			     iterate.orbcoeff, iterate.orben)
end

