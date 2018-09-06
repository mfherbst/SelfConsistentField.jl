using HDF5
using LinearAlgebra

"""
Take a matrix which has either one or two spin components
and make an α-β block diagonal matrix out of it
"""
function build_ab_matrix(mat)
    n_bas = size(mat, 1)

    Mα = Mβ = nothing
    if ndims(mat) == 2
        Mα = view(mat, :, :)
        Mβ = view(mat, :, :)
    elseif ndims(mat) == 3
        n_spin = size(mat, 3)
        @assert n_spin == 1 || n_spin == 2

        Mα = view(mat, :, :, 1)
        Mβ = view(mat, :, :, n_spin)
    end
    @assert Mα != nothing && Mβ != nothing

    ret = [ Mα 0*I; 0*I Mβ ]
    @assert size(ret) == (2 * n_bas, 2 * n_bas)
    return ret
end

"""
Take a matrix which has either one or two spin components
and concatenate the  α and β block after another
"""
function concat_spin(mat)
    Mα = Mβ = nothing
    if ndims(mat) == 2
        n_spin = size(mat, 2)
        @assert n_spin == 1 || n_spin == 2

        Mα = view(mat, :, 1)
        Mβ = view(mat, :, n_spin)
    elseif ndims(mat) == 3
        n_spin = size(mat, 3)
        @assert n_spin == 1 || n_spin == 2

        Mα = view(mat, :, :, 1)
        Mβ = view(mat, :, :, n_spin)
    end
    @assert Mα != nothing && Mβ != nothing
    [ Mα  Mβ ]
end

"""
Take a matrix in AOs and perform a transformation into the
molecular orbital basis
"""
function matrix_b2f(mat, orbcoeff)
    n_bas, n_orb, n_spin = size(orbcoeff)

    if ndims(mat) == 2
        if n_spin == 1
            mat = cat(mat; dims=3)
        elseif n_spin == 2
            mat = cat(mat, mat; dims=3)
        end
    end
    @assert size(mat) == (n_bas, n_bas, n_spin)

    mat_ff = Array{eltype(mat)}(undef, n_orb, n_orb, n_spin)
    for σ in 1:n_spin
        mat_ff[:,:,σ] = (view(orbcoeff, :, :, σ)' * view(mat, :, :, σ) *
                         view(orbcoeff, :, :, σ))
    end
    return build_ab_matrix(mat_ff)
end


"""
Dump the iterate returned from the run_scf
method, resembling the format of the molsturm HDF5
file as of molsturm 0.0.3 as much as possible
"""
function dump_molsturm_hdf5(integrals, scfres, file;
                            export_eri_bbbb=false, export_eri_ffff=false)
    n_bas = size(scfres["fock"], 1)
    problem = scfres["problem"]
    energies = scfres["energies"]
    energy_change = scfres["energy_change"]

    h5open(file, "w") do file
        attrs(file)["format_version"] = "0.1.0"
        attrs(file)["format_type"] = "hdf5"

        file["n_alpha"] = problem.n_elec[1]
        file["n_beta"] = problem.n_elec[2]
        file["n_bas"] = n_bas
        file["restricted"] = UInt8(problem.restricted)
        file["overlap_bb"] = build_ab_matrix(problem.overlap)

        file["n_orbs_alpha"] = problem.n_orb
        file["n_orbs_beta"] = problem.n_orb
        file["orben_f"] = concat_spin(scfres["orben"])
        file["fock_bb"] = build_ab_matrix(scfres["fock"])
        file["fock_ff"] = matrix_b2f(scfres["fock"], scfres["orbcoeff"])
        file["hcore_ff"] = matrix_b2f(problem.h_core, scfres["orbcoeff"])

        # orbcoeff is not a symmetric matrix, so the difference
        # in ordering between Julia and python comes into play
        # (Julia is column-major and python row-major),
        # such that we have to transpose.
        file["orbcoeff_bf"] = Array(adjoint(concat_spin(scfres["orbcoeff"])))

        file["n_mtx_applies"] = scfres["n_applies"]
        file["n_iter"] = scfres["n_iter"]

        file["energy_coulomb"] = energies["coulomb"]
        file["energy_exchange"] = energies["exchange"]
        file["energy_kinetic"] = energies["kinetic"]
        file["energy_nuclear_attraction"] = energies["nuclear_attraction"]
        file["energy_nuclear_repulsion"] = energies["nuclear_repulsion"]

        file["energy_ground_state"] = energies["total"]
        file["final_1e_energy_change"] = energy_change["1e"]
        file["final_error_norm"] = scfres["error_norm"]
        file["final_tot_energy_change"] = energy_change["total"]

        # TODO The following parameters are missing, that a user of
        #      the molsturm state could expect.
        #    input_parameters
        #    spin_squared
        #    eri_ffff
        #    overlap_ff
    end

    if export_eri_bbbb
        error("Exporting eri in basis functions not yet implemented.")
    end
    if export_eri_ffff
        error("Exporting eri in orbitals not yet implemented.")
    end

    sanitise_hdf5_hack(file)
end

"""
For some reason the hdf5 file created using dump_state_molsturm
is not compatible with molsturm per se. Probably it's actually an issue
on the python side. Until this is resolved we call python to explicitly
sanitise the file after it has been written.
"""
function sanitise_hdf5_hack(file)
    # The inner python code to execute

    code = """
    import h5py
    import numpy as np

    with h5py.File("$file") as f:
        # There is an issue with storing these attribute strings from julia,
        # so overwrite explicitly here.
        f.attrs["format_version"] = "0.1.0"
        f.attrs["format_type"] = "hdf5"

        # These keys are bool and need to be explicitly converted so
        for k in ["restricted"]:
            if k in f:
                f.create_dataset(k, data=f.pop(k), dtype=np.dtype("b1"))
    """

    # Create python Cmd object
    pycmd = `python3 -c "$code"`
    run(pycmd)
end

