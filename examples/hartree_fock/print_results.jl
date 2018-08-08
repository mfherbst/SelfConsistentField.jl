using LinearAlgebra: norm, tr

function print_energies(problem, integrals, res)
    _, _, n_spin = size(res["density"])

    # Alpha density
    Da = view(res["density"], :, :, 1)
    Db = if n_spin == 1 Da else view(res["density"], :, :, 2) end

    # Compute energies
    # TODO This should live somewhere else
    Ekin = tr(Da * integrals.kinetic_bb) + tr(Db * integrals.kinetic_bb)
    Enucattr = (tr(Da * integrals.nuclear_attraction_bb)
                + tr(Db * integrals.nuclear_attraction_bb))
    e0e = res["energies"]["0e"]
    e1e = res["energies"]["1e"]
    e2e = res["energies"]["2e"]

    println("Final energies")
    println("kinetic            ", Ekin)
    println("nuclear_attraction ", Enucattr)
    println("nuclear_repulsion  ", e0e)
    println()

    println("E_1e               ", e1e)
    println("E_2e               ", e2e)
    println("E electronic       ", e1e + e2e)

    println()
    Epot = Enucattr + e2e + e0e
    println("E_pot              ", Epot)
    println("E_kin              ", Ekin)
    println("virial ratio       ", -Epot / Ekin)

    println()
    println("E_total            ", res["energies"]["total"])
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
