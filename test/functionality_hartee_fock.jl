this_dir = dirname(PROGRAM_FILE)

function run_problem(intfile, restricted)
    intfile = joinpath(this_dir, "data", intfile)
    system, integrals = load_integral_hdf5(intfile)
    problem = assemble_hf_problem(system, integrals; restricted=restricted)
    guess_density = compute_guess(problem, system.coords, system.atom_numbers;
                                  method="hcore")  # TODO better use random here

    if ! problem.restricted && is_closed_shell(problem)
        break_spin_symmetry(guess_density)
    end
    return run_scf(problem, guess_density)
end

begin
    rhf_be_321g = run_problem("integrals_be_3-21g.hdf5", true)
    @test rhf_be_321g["energies"]["total"] ≈ -14.4868202421763 rtol=1e-9
    @test rhf_be_321g["converged"]
end

begin
    uhf_be_321g = run_problem("integrals_be_3-21g.hdf5", false)
    @test uhf_be_321g["energies"]["total"] ≈ -14.4868202421763 rtol=1e-9
    @test uhf_be_321g["converged"]
end

begin
    rohf_c_321g = run_problem("integrals_c_3-21g.hdf5", true)
    @test rohf_c_321g["energies"]["total"] ≈ -37.4803888099046 rtol=1e-9
    @test rohf_c_321g["converged"]
end

begin
    uhf_c_321g = run_problem("integrals_c_3-21g.hdf5", false)
    @test uhf_c_321g["energies"]["total"] ≈ -37.4810698325847 rtol=1e-9
    @test uhf_c_321g["converged"]
end
