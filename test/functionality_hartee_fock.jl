this_dir = dirname(PROGRAM_FILE)
include(joinpath(this_dir, "util", "load_integral_file.jl"))

function run_problem(intfile, restricted)
    intfile = joinpath(this_dir, "data", intfile)
    system, integrals = load_integral_file(intfile)
    problem = assemble_hf_problem(system, integrals; restricted=restricted)
    hcore_guess_density = compute_guess_hcore(problem, system.coords, system.atom_numbers)
    return run_scf(problem, hcore_guess_density)
end

rhf_be_321g = run_problem("integrals_be_3-21g.hdf5", true)
@test rhf_be_321g["energies"]["total"] ≈ -14.4868202421763 atol=1e-10

uhf_be_321g = run_problem("integrals_be_3-21g.hdf5", false)
@test uhf_be_321g["energies"]["total"] ≈ -14.4868202421763 atol=1e-10

rohf_c_321g = run_problem("integrals_c_3-21g.hdf5", true)
@test rohf_c_321g["energies"]["total"] ≈ -37.4803888099046 atol=1e-10

uhf_c_321g = run_problem("integrals_c_3-21g.hdf5", false)
@test uhf_c_321g["energies"]["total"] ≈ -37.4810698325847 atol=1e-10
