struct ComputeFock <: Algorithm end

function iterate(::ComputeFock, rp::SubReport)
    lg = Logger(rp)
    log!(lg, "Calculating new fock matrix", :debug, :computefock)

    fock, error_pulay, energies = compute_fock_matrix(rp.problem, rp.state.density)
    log!(lg, "new fock matrix", fock, :debug, :computefock)
    log!(lg, "new error_pulay", error_pulay, :debug, :computefock)
    log!(lg, "new energies", energies, :debug, :computefock)

    state = FockIterState(fock, error_pulay, energies, rp.state.orbcoeff, rp.state.orben, rp.state.density)
    return ComputeFock(), new_subreport(ComputeFock(), state, lg, rp)
end

"""
Build Fock matrix or Kohn-Sham matrix along with the associated energies
and errors for the density passed as second argument.

For restricted open-shell and unrestricted an alpha and a beta density are expected,
for restricted closed-shell HF a single density is expected.
"""
function compute_fock_matrix(problem::ScfProblem, density::AbstractArray; kwargs...)
    # A lot of the n_spin stuff could go if a sensible data structure
    # was used for the density, which could implicitly take care of the
    # alpha and beta density

    overlap = problem.overlap
    _, _, n_spin = size(density)
    h_core = problem.h_core
    n_bas = size(h_core)[1]

    @assert size(density) == (n_bas, n_bas, n_spin)
    @assert n_spin == 1 || n_spin == 2

    # Compute total density
    total_density=0
    if n_spin == 1
        total_density = 2 * view(density, :, :, 1)
    elseif n_spin == 2
        total_density = view(density, :, :, 1) + view(density, :, :, 2)
    end

    # Compute 2 electron terms
    G = zeros(n_bas, n_bas, n_spin)
    for (label, builder) in problem.terms_2e
        # Build the term and add it to G
        add_2e_term!(G, builder, density, total_density)
    end
    @assert size(G) == (n_bas, n_bas, n_spin)

    # Compute energies
    energies = Dict(
        "0e" => sum(v for (k,v) in problem.terms_0e),
        "1e" => tr(h_core * total_density)
    )

    # Form tr(G*density) for alpha and beta separately
    if n_spin == 1
        # We only track the alpha component and implicitly have alpha = beta,
        # thus we need to multiply by 2 to make up for this.
        # Additionally the two-electron energy gets a factor of 1/2 to avoid
        # double counting, this these two factors just cancel.
        energies["2e"] = tr(view(G, :,:,1) * view(density, :,:,1)) # * 2 / 2
    elseif n_spin == 2
        energies["2e"] = 1/2 * (tr(view(G, :,:,1) * view(density,:,:,1))
                                + tr(view(G, :,:,2) * view(density,:,:,2)))
    end
    energies["total"] = energies["0e"] + energies["1e"] + energies["2e"]

    # Build fock matrix (notice that broadcasting will take care of the
    # deviating shapes in case different spin components are present)
    fock = G .+ h_core

    error_pulay = 0
    if n_spin == 2 && problem.restricted
        # Build effective fock matrix for restricted open-shell
        fock = compute_roothaan_effective_fock(fock, density, overlap)
        error_pulay = compute_pulay_error(fock, total_density, overlap)
        fock = reshape(fock, n_bas, n_bas, 1)
    else
        error_pulay = compute_pulay_error(fock, density, overlap)
    end

    return fock, error_pulay, energies
end

"""
Function to build Roothaan's effective Fock  matrix, reference
http://www-theor.ch.cam.ac.uk/people/ross/thesis/node15.html
and pyscf ROHF procedure

fock         alpha, beta Fock matrix
density      alpha, beta density matrix
overlap      Overlap matrix
"""
function compute_roothaan_effective_fock(fock, density, overlap)
    n_bas, _, n_spin = size(density)

    @assert size(density) == (n_bas, n_bas, n_spin)
    @assert size(fock) == (n_bas, n_bas, n_spin)
    @assert size(overlap) == (n_bas, n_bas)

    Da = view(density, :, :, 1)
    Db = view(density, :, :, 2)
    Fa = view(fock, :, :, 1)
    Fb = view(fock, :, :, 2)
    S = overlap

    # The effective ROHF Fock matrix in MOs is constructed as
    #
    #          |  closed   open    virtual
    #  ----------------------------------------
    #  closed  |    Fc      Fb       Fc
    #  open    |    Fb      Fc       Fa
    #  virtual |    Fc      Fa       Fc
    #  ----------------------------------------
    #
    # where
    Fc = (Fa + Fb) / 2

    # Projectors for closed-shell, open-shell and virtual spaces
    Pc = Db * S
    Po = (Da - Db) * S
    Pv = I - Da * S

    # Now build the effective Fock matrix. The idea is to
    # build the lower triangle only and than use the hermitian
    # conjugate to build the full thing. Implies that a factor
    # 1/2 is needed along the diagonal

    Feff = (
          Pc' * Fc * Pc / 2
        + Po' * Fb * Pc     + Po' * Fc * Po / 2
        + Pv' * Fc * Pc     + Pv' * Fa * Po     + Pv' * Fc * Pv / 2
    )
    return Feff + Feff'
end
