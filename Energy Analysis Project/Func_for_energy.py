from OUR_basic_setup import*
import math
import numpy as np


def residue_id(res):  # residue as ASP A32
    # res.id[1] is integer, we should use str to get the corresponding string
    return res.get_resname() + '.' + res.get_parent().id + str(res.id[1])


####################################
#### Getting Interface Residues ####
####################################

# Generate 2 sets for the residues of both chains

resA = set()
resE = set()
bond_pairs = set()
all_interface_atoms = set()


# Generate a list for the atoms
select = []

# Add all atoms of the structure to the list
for at in st.get_atoms():
    select.append(at)

# Prepare Neighbor Search
nbsearch = NeighborSearch(select)


# By a neighbor search, check if the atoms are closer than given distance angstrom

for at1, at2 in sorted(nbsearch.search_all(args.cutoff)):
    if at1.get_parent().get_parent() != at2.get_parent().get_parent():
        resA.add(at1.get_parent())
        resE.add(at2.get_parent())
        all_interface_atoms.add(at1)
        all_interface_atoms.add(at2)
        bond_pairs.add((at1.get_parent(),
                        at2.get_parent()))

# Generate a list for the atoms
ALL_interface_res = list(resA)+list(resE)


####################################
#### Getting Solvation Energies ####
####################################
# maybe we should write a function for getting solvation of two residues because we need that in mutation maybe?
#  or the script i did below will be able to calculte it? idk
# ΔG(SolvA-B) - ΔG(SolvA) - ΔG(SolvB)


def solv_energy_one_atom(atom, structure):
    solv_energy = 0
    atom_name = atom.id
    residue = atom.get_parent()
    res_id = residue.get_id()[1]
    chain = residue.get_parent().get_id()

    at = structure[0][chain][res_id][atom_name]
    if 'EXP_NACCESS' in at.xtra:
        asa = at.xtra['EXP_NACCESS']
        S_E = atom.xtra['vdw'].fsrf*float(asa)
        solv_energy = S_E
    return solv_energy


####################################
#### Electrostatic interaction  ####
####################################


def calc_elec_int2(at1, at2):
    dist = at2-at1
    E = 0
    if dist > 0:

        Mehler_dielectric = (
            (86.9525)/(1 - (7.7839*math.exp(-0.3153*dist))))-8.85525
        E = 332.16 * \
            ((at1.xtra['charge'] * at2.xtra['charge']) /
                (Mehler_dielectric*dist))
    return E


####################################
#### VanderWaal Energies        ####
####################################

def vdw_energy(at1, at2):
    vdw = 0
    r = at2-at1
    if r > 0:
        eps1 = at1.xtra['vdw'].eps
        sig1 = at1.xtra['vdw'].sig
        eps2 = at2.xtra['vdw'].eps
        sig2 = at2.xtra['vdw'].sig
        epsilon = np.sqrt(eps1 * eps2)
        sigma = np.sqrt(sig1 * sig2)
        vdw += 4*epsilon*((sigma/r)**12 - (sigma/r)**6)
    return vdw


####################################
#### Total Interaction energy   ####
####################################
# ΔG A-B = ΔG (electA-B) + ΔG (vdwA-B) + ΔG (SolvA-B) - ΔG (SolvA) - ΔG (SolvB)
def Total_interaction_Energies(resA, resE):
    vdw = 0
    electr_inter = 0
    sol_AB = 0
    sol_A = 0
    sol_B = 0

    # computing electrostatics
    for res1 in resA:
        for res2 in resE:
            for atom1 in res1.get_atoms():
                for atom2 in res2.get_atoms():
                    vdw += vdw_energy(atom1, atom2)
                    electr_inter += calc_elec_int2(atom1, atom2)

    # computing solvation energies
    for res in resA:
        for atom in res.get_atoms():
            sol_AB += solv_energy_one_atom(atom, st)

            sol_A += solv_energy_one_atom(atom, st_A)

    for res in resE:

        for atom in res.get_atoms():
            sol_AB += solv_energy_one_atom(atom, st)
            sol_B += solv_energy_one_atom(atom, st_B)
    delta_G = electr_inter + vdw + sol_AB - sol_A - sol_B
    print('electrostatic:', electr_inter)
    print('energy vdw:', vdw)
    print('solvation AB:', sol_AB)
    print('solvation A:', sol_A)
    print('solvation B:', sol_B)
    return delta_G


####################################################
### change in solvation energies due to mutations ##
####################################################
'''
Determine the effect of replacing each interface residue by Ala in the overall ΔG A-B, and make a
plot of the results obtained, highlighting those residues that are more relevant for the stability
of the interface. Discuss the results obtained in relation to the nature of the involved amino
acids(hint, as Ala side-chain is part of all others, except Gly, there no need to do the
      replacement, just take into account the different atoms)
'''


ALA = ['N', 'CA', 'C', 'O', 'CB']


def Energy_mutation_chainA(resA, resE, mut_index):
    vdw = 0
    electr_inter = 0
    solvation = {}
    sol_AB = 0
    sol_A = 0
    sol_B = 0

    # computing electrostatics
    res_tobe_mutated = resA[mut_index]
    for res1 in resA:
        for res2 in resE:

            for atom2 in res2.get_atoms():
                if atom2 not in solvation:
                    AB = solv_energy_one_atom(atom2, st)
                    B = solv_energy_one_atom(atom2, st_B)
                    solvation[atom2] = [AB, B]

                if res1 == res_tobe_mutated:
                    for atom1 in res1.get_atoms():
                        if atom1.id in ALA:
                            if atom1 not in solvation:
                                AB = solv_energy_one_atom(atom1, st)
                                A = solv_energy_one_atom(atom1, st_A)
                                solvation[atom1] = [AB, A]
                            vdw += vdw_energy(atom1, atom2)
                            electr_inter += calc_elec_int2(atom1, atom2)
                else:

                    for atom1 in res1.get_atoms():
                        if atom1 not in solvation:
                            AB = solv_energy_one_atom(atom1, st)
                            A = solv_energy_one_atom(atom1, st_A)
                            solvation[atom1] = [AB, A]
                        vdw += vdw_energy(atom1, atom2)
                        electr_inter += calc_elec_int2(atom1, atom2)
    for atom in solvation:

        sol_AB += solvation[atom][0]
        res = atom.get_parent()
        if res in resA:
            sol_A += solvation[atom][1]
        if res in resE:
            sol_B += solvation[atom][1]
    delta_G = electr_inter + vdw + sol_AB - sol_A - sol_B
    return delta_G


def Energy_mutation_chainB(resA, resE, mut_index):
    vdw = 0
    electr_inter = 0
    solvation = {}
    sol_AB = 0
    sol_A = 0
    sol_B = 0

    # computing electrostatics
    res_tobe_mutated = resE[mut_index]
    for res1 in resE:
        for res2 in resA:

            for atom2 in res2.get_atoms():
                if atom2 not in solvation:
                    AB = solv_energy_one_atom(atom2, st)
                    B = solv_energy_one_atom(atom2, st_A)
                    solvation[atom2] = [AB, B]

                if res1 == res_tobe_mutated:
                    for atom1 in res1.get_atoms():
                        if atom1.id in ALA:
                            if atom1 not in solvation:
                                AB = solv_energy_one_atom(atom1, st)
                                A = solv_energy_one_atom(atom1, st_B)
                                solvation[atom1] = [AB, A]
                            vdw += vdw_energy(atom1, atom2)
                            electr_inter += calc_elec_int2(atom1, atom2)
                else:

                    for atom1 in res1.get_atoms():
                        if atom1 not in solvation:
                            AB = solv_energy_one_atom(atom1, st)
                            A = solv_energy_one_atom(atom1, st_B)
                            solvation[atom1] = [AB, A]
                        vdw += vdw_energy(atom1, atom2)
                        electr_inter += calc_elec_int2(atom1, atom2)
    for atom in solvation:

        sol_AB += solvation[atom][0]
        res = atom.get_parent()
        if res in resA:
            sol_A += solvation[atom][1]
        if res in resE:
            sol_B += solvation[atom][1]
    delta_G = electr_inter + vdw + sol_AB - sol_A - sol_B
    return delta_G
