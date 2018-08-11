using LinearAlgebra: norm, tr
using Printf

function print_results(problem, res)
    if res["converged"]
        println("SCF converged.")
    else
        println("SCF failed to converge")
    end

    println()
    print_energies(res)
    println()
    print_mo_occupation(problem, res)
end

function print_energies(res)
    energies = res["energies"]

    println("Final energies")
    @printf("%20s = %15.10g\n", "coulomb", energies["coulomb"])
    @printf("%20s = %15.10g\n", "exchange", energies["exchange"])
    @printf("%20s = %15.10g\n", "kinetic", energies["kinetic"])
    @printf("%20s = %15.10g\n", "nuclear_attraction", energies["nuclear_attraction"])
    @printf("%20s = %15.10g\n", "nuclear_repulsion", energies["nuclear_repulsion"])
    println()

    @printf("%20s = %15.10g\n", "E_1e", energies["1e"])
    @printf("%20s = %15.10g\n", "E_2e", energies["2e"])
    @printf("%20s = %15.10g\n", "E electronic", energies["1e"] + energies["2e"])

    println()
    Epot = energies["nuclear_attraction"] + energies["0e"] + energies["2e"]
    @printf("%20s = %15.10g\n", "E_pot", Epot)
    @printf("%20s = %15.10g\n", "E_kin", energies["kinetic"])
    @printf("%20s = %15.10g\n", "virial ratio", -Epot / energies["kinetic"])

    println()
    @printf("%20s = %20.15g\n", "E_total", energies["total"])
end

function print_mo_occupation(problem, res)
    n_orb, n_spin = size(res["orben"])

    println("Orbital occupation")
    println("a                             b")
    for i in 1:n_orb
        aocc = bocc = " "
        if i <= problem.n_elec[1] aocc = "*" end
        if i <= problem.n_elec[2] bocc = "*" end
        if n_spin == 1
            ene = res["orben"][i]
            println("$aocc       $ene       $bocc")
        else
            enea = res["orben"][i, 1]
            eneb = res["orben"][i, 2]
            println("$aocc       $enea    |    $eneb       $bocc")
        end
    end
end
