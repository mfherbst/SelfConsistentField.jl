"""Type for the iterate objects helping to build two-electron terms"""
abstract type TwoElectronBuilder end

function add_2e_term!(out::AbstractArray,
                      builder::TwoElectronBuilder,
                      density::AbstractArray,
                      total_density::AbstractArray)
    error("Please overload add_2e_term! for your TwoElectronBuilder type")
end

"""
Type for defining a problem for an SCF
"""
struct ScfProblem
    # -----------------------
    # System information

    # Overlap matrix
    overlap::AbstractMatrix

    # Number of occupied orbitals in alpha and beta
    n_elec::Tuple{Int,Int}

    # Number of orbitals to be computed (same for alpha and beta)
    n_orb::Int

    # Should the SCF be restricted or not
    restricted::Bool

    # -----------------------
    # Operator terms

    # Zero-electron terms
    terms_0e::Dict{String,Number}

    # One-electron terms
    terms_1e::Dict{String,AbstractMatrix}

    # Coulomb and exchang term
    terms_2e::Dict{String,TwoElectronBuilder}

    # -----------------------
    # Caches

    # Core Hamiltonian (i.e. sum of all terms_1e)
    h_core::AbstractMatrix  # TODO Should be explicitly made a stored array

    # TODO Probably some factorisation of the overlap
    #      matrix is sensible, since it is used on the RHS a lot.

    # -----------------------
    # Functions

    # computes next fock, error_pulay and energies
    compute_fock_matrix::Function

    """
    Setup an ScfProblem taking some system-specific data as well as the terms
    of the problem matrix.

    overlap        Overlap matrix
    n_elec          Tuple of alpha and beta occupied electrons.
    n_orb          Number of orbitals to be computed (same for alpha and beta)
    restricted     Should a restricted (true) or an unrestricted (false) SCF be run
    terms_0e       The energies of the non-electron dependent terms (e.g. nuclear attraction)
    terms_1e       The matrices representing the density-independent one-electron terms
                   (e.g. nuclear attraction and electronic kinetic energy)
    terms_2e       The two-electron terms, which need to be updated for each SCF iteration
                   with a new density.
    """
    function ScfProblem(overlap::AbstractMatrix, n_elec::Tuple{Int,Int}, n_orb::Int,
                        restricted::Bool,
                        terms_0e::Dict{String,Number},
                        terms_1e::Dict{String,AbstractMatrix},
                        terms_2e::Dict{String,TwoElectronBuilder})
        @assert length(n_elec) == 2 "n_elec needs to have exactly two elements"

        n_bas, _ = size(overlap)
        @assert size(overlap) == (n_bas, n_bas)  "Overlap must be a quadratic matrix"

        h_core = zeros(eltype(overlap), n_bas, n_bas)
        for (label, term) in terms_1e
            @assert size(term) == (n_bas, n_bas) begin
                "Term '$label ' must have a shape of $n_bas x $n_bas"
            end
            h_core .+= term
        end

        new(overlap, n_elec, n_orb, restricted, terms_0e, terms_1e, terms_2e, h_core)
    end
end
