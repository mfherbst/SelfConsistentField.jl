"""
Type describing the SCF quantity to be iterated by the SCF map
(e.g. Fock matrix)
"""
abstract type ScfIterState end

"""
Struct of the state if the Fock matrix is iterated
"""
mutable struct FockIterState <: ScfIterState
    # Fock or Kohn-Sham matrix
    fock::AbstractArray

    # Matrix describing the Pulay error in fock
    error_pulay::AbstractArray

    # Named tuple or dict of the computed energy terms
    energies::Dict{String, Float64}

    # Orbital coefficients and energies
    orbcoeff::Union{Nothing,AbstractArray}
    orben::Union{Nothing,AbstractArray}

    # Density
    density::Union{Nothing,AbstractArray}
end

function update_fock_matrix(iterate::FockIterState, fock::AbstractArray)
    return FockIterState(fock, iterate.error_pulay, iterate.energies,
                         iterate.orbcoeff, iterate.orben, iterate.density)
end

