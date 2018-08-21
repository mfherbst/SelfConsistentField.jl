"""
Compute nuclear repulsion energy from a system.

The system is expected to contain the members
     coords (n_atoms × 3 array)     atom positions
     atom_numbers  (n_atoms array)  atomic charges
"""
function compute_nuclear_repulsion(system)
    sum = 0
    for (A, ZA) in enumerate(system.atom_numbers)
        for B in 1:A-1
            dist = norm(system.coords[A, :] - system.coords[B, :])
            ZB = system.atom_numbers[B]
            sum += ZA * ZB / dist
        end
    end
    return sum
end

"""
Build a standard Hartree-Fock problem
Compute nuclear repulsion energy from a system.

The system object is expected to contain the members
     coords (n_atoms × 3 array)     atom positions
     atom_numbers  (n_atoms array)  atomic charges
     n_elec                         number of electrons as tuple α, β

The integrals object is expected to contain the members
     overlap_bb       Overlap matrix object
     kinetic_bb       Kinetic energy matrix object
     nuclear_attraction_bb     Nuclear attraction potential matrix object
     jk_builder_bb    TwoElectronBuilder for Coulomb+Exchange matrix (as a sum)

Extra terms can be provided as kwargs to the function call
"""
function assemble_hf_problem(system, integrals; n_orb=nothing, restricted=nothing,
                             kwargs...)
    overlap = integrals.overlap_bb

    if restricted == nothing
        # true if same number of electrons, else false
        restricted = is_closed_shell(system)
    end

    if n_orb == nothing
        n_orb = size(overlap)[1]
    end

    terms_0e = Dict{String, Number}(
        "nuclear_repulsion" => compute_nuclear_repulsion(system),
    )
    terms_1e = Dict{String, AbstractMatrix}(
        "kinetic" => integrals.kinetic_bb,
        "nuclear_attraction" => integrals.nuclear_attraction_bb,
    )
    terms_2e = Dict{String, TwoElectronBuilder}(
        "coulomb+exchange" => integrals.jk_builder_bb,
    )

    # Deal with extra terms
    for (symbol, term) in kwargs
        # Convert kwargs symbol to label for the term
        label = String(symbol)
        if isa(term, TwoElectronBuilder)
            terms_2e[label] = term
        elseif isa(term, AbstractMatrix)
            terms_1e[label] = term
        elseif isa(term, Number)
            terms_0e[label] = term
        else
            @assert false "the term labelled $label does not have a recognised type."
        end
    end

    return ScfProblem(overlap, system.n_elec, n_orb, restricted,
                      terms_0e, terms_1e, terms_2e, compute_fock_matrix)
end

