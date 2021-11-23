from OUR_basic_setup import*
import math
from math import np


def residue_id(res):  # residue as ASP A32
    # res.id[1] is integer, we should use str to get the corresponding string
    return res.get_resname() + '.' + res.get_parent().id + str(res.id[1])


####################################
#### Getting Interface Residues ####
####################################

# Generate 2 sets for the residues of both chains
resA = set()
resE = set()

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


####################################
#### Getting Solvation Energies ####
####################################
# maybe we should write a function for getting solvation of two residues because we need that in mutation maybe?
#  or the script i did below will be able to calculte it? idk


def solv_energy(residues_list):
    solv_energies = {}
    Total_solvE = 0
    for res in residues_list:
        delG_solv_residue = 0
        for atom in res.get_atoms():
            if 'EXP_NACCESS' in atom.xtra:
                S_E = atom.xtra['vdw'].fsrf * \
                    float(atom.xtra['EXP_NACCESS'])
                delG_solv_residue += S_E
                Total_solvE += S_E
        solv_energies[residue_id(res)] = delG_solv_residue
    return solv_energies, Total_solvE


solv_res_ChainA, T_A = solv_energy(resA)
solv_res_ChainB, T_B = solv_energy(resE)

print('\nTotal solavtion energy of chain A = ', T_A)
print('\nSolvation energy for each residue of Chain A :')
for ener in solv_res_ChainA:
    print(ener, ':', solv_res_ChainA[ener])
print('\nTotal solavtion energy of chain B = ', T_B)
print('\nSolvation energy for each residue of Chain B :')
for ener in solv_res_ChainB:
    print(ener, ':', solv_res_ChainB[ener])


print('\nSolvation energy of both chains : ', T_A+T_B)

print('\ndo we need Solvation energy of whole complex too? i.e not only of interface?')


#### Solvation ENergy of Complex ####

print('WHOLE COMPLEX NOT ONLY INTERRFACE')
ALL_RESIDUES = []
for res in st.get_residues():
    ALL_RESIDUES.append(res)
solv_all, T_all = solv_energy(ALL_RESIDUES)
for ener in solv_all:
    print(ener, ':', solv_all[ener])

print('TOTAL')
print(T_all)


####################################
#### Electrostatic interaction  ####
####################################

# THIS IS NOT CORRECT


# calculate total electrostaic .... do we need electrostatic of each atom too?
def electrostatics(atom1, atom2):
    # which r is used in this constant????
    #Mehler_dielectric = ((86.9525)/(1 - 7.7839*math.exp(-0.3153*r_ij)))-8.85525
    #E = 332.16*((qi*qj)/Mehler_dielectric*r_ij)

    return 'ENERGY'


def calc_elec_int2(at1, at2):
    dist = at2-at1
    Mehler_dielectric = ((86.9525)/(1 - 7.7839*math.exp(-0.3153*dist)))-8.85525
    E = 332.16*((at1.xtra['charge'] * at2.xtra['charge']) /
                (Mehler_dielectric*dist))
    return E

####################################
#### VanderWaal Energies        ####
####################################


# do something like


def vdw_energy(atom1, atom2):
    return 0


def distance(at1, at2):
    return at2-at1

# THIS IS WRONG


def vdw_parameters(epsilon_i, sigma_i, epsilon_j, sigma_j):
    Ai = 2*math.sqrt(epsilon_i) * pow(sigma_i, 6)
    Aj = 2*math.sqrt(epsilon_j) * pow(sigma_j, 6)
    Ci = 2*math.sqrt(epsilon_i) * pow(sigma_i, 3)
    Cj = 2*math.sqrt(epsilon_j) * pow(sigma_j, 3)
    return Ai, Aj, Ci, Cj


def vdw_energy(epsilon_i, sigma_i, epsilon_j, sigma_j, rij):
    Ai, Aj, Ci, Cj = vdw_parameters(epsilon_i, sigma_i, epsilon_j, sigma_j)
    t1 = (Ai*Aj)/pow(rij, 12)
    t2 = (Ci*Cj)/pow(rij, 6)
    vdw_E = t1-t2
    return vdw_E

# maybe the one below is right , but do you remeber in class he advised us to use
# formula having A and C? and he said we have to calculate A and C for pairs
# and we can use the same again and again whenever we get the same pair?


def vdw_energy(at1, at2):
    vdw = 0
    r = at2-at1
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


# dont pass energies as parameters call their functions instead
def delta_G_AB(electro_AB, vander_AB, solv_AB, solv_A, solv_B):
    # electrostatic and vander interaction between both chains will be electrostatic+vander of all atoms in interface of chain A and B
    return 0


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

# use function delta_G_AB in function mutation
# print delta_G_AB before mutation : delta_G_AB after mutation
# to calculate delta G mutation w ehave to chande one residue in residue list by alanine and calculate the total energy of chain and sum it with total energy of another chain
# so we will calculate delta G more than 100 time if w etake cutoff distance 8 for interface but its oky!!

# residdue list should contain all residues maybe


ALL_interface_res = resA+resE


ENERGIES_AFTER_MUTATION = {}
for i in range(len(ALL_interface_res)):

    Mutant = ALL_interface_res[i]
    ENERGIES_AFTER_MUTATION.append(Mutant)
    ALL_interface_res[i] = 'ALANINE'
    delta_G = delta_G_AB()
    ENERGIES_AFTER_MUTATION[Mutant] = delta_G


def Mutations(residue_list):
    '''
    f
    for residue in residue list

    '''

    return 0
