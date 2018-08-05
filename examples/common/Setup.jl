using HDF5

struct Integrals
	kinetic_bb
	nuclear_attraction_bb
	overlap_bb
	electron_repulsion_bbbb
end

"""
Struct describing a system to be modelled
"""
struct System
	# Coordinates of the atoms, as array
	# n_atoms times 3
	coords::Array{Float64, 2}

	# Atom numbers
	atom_numbers::Array{Float64, 1}

	# Number of electrons
	nelecs::Array{Int, 1}  #Tuple: alpha, beta
end


"""
Compute nuclear repulsion energy
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
Read an hdf5 file
"""
function read_hdf5(file)
	file = h5open(file)

	eri = read(file, "electron_repulsion_bbbb")
	T = read(file,"kinetic_bb")
	V = read(file,"nuclear_attraction_bb")
	S = read(file,"overlap_bb")

	nelec = read(file,"system/nelec")
	atnums = read(file, "system/atom_numbers")
	coords = read(file, "system/coords")

	close(file)

	# coordinates need to be transformed from row-major
	# to column-major
	coords = coords'

	system = System(coords, atnums, nelec)
	integrals = Integrals(T, V, S, eri)
	return system, integrals
end

function print_energies(problem, integrals, res)
	_, _, n_spin = size(res["density"])

	# Alpha density
	Da = view(res["density"], :, :, 1)
	Db = if n_spin == 1 Da else view(res["density"], :, :, 2) end

	# Compute energies
	# TODO This should live somewhere else
	Ekin = trace(Da * integrals.kinetic_bb) + trace(Db * integrals.kinetic_bb)
	Enucattr = (trace(Da * integrals.nuclear_attraction_bb)
		    + trace(Db * integrals.nuclear_attraction_bb))

	println("Final energies")
	println("kinetic            ", Ekin)
	println("nuclear_attraction ", Enucattr)
	println("nuclear_repulsion  ", problem.energy_nuc_rep)
	println()

	e1e = res["energies"]["energy_1e"]
	e2e = res["energies"]["energy_2e"]
	println("E_1e               ", e1e)
	println("E_2e               ", e2e)
	println("E electronic       ", e1e + e2e)

	println()
	Epot = Enucattr + e2e + problem.energy_nuc_rep
	println("E_pot              ", Epot)
	println("E_kin              ", Ekin)
	println("virial ratio       ", -Epot / Ekin)

	println()
	println("E_total            ", res["energies"]["energy_total"])
end

function print_mo_occupation(problem, res)
	println("Orbital occupation")
	println("a                             b")
	for i in 1:length(res["orben"])
		aocc = bocc = " "
		if i <= problem.n_occ
			aocc = bocc = "*"
		end
		ene = res["orben"][i]
		println("$aocc     $ene       $bocc")
	end
end
