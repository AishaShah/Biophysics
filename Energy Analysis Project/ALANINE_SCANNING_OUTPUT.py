from tabulate import tabulate
import matplotlib.pyplot as plt
from Func_for_energy import*

####### TOTAL INTERACTION ENERGY ########
E = Total_interaction_Energies(resA, resE)

print('Total Interaction Energy')
print('ΔG A-B = ΔG (electA-B) + ΔG (vdwA-B) + ΔG (SolvA-B) - ΔG (SolvA) - ΔG (SolvB)')
print('ΔG A-B = ',  E)


Energy_list = {}


resA = list(resA)
resE = list(resE)
for i in range(len(resA)):

    Energy_After_Mutation = Energy_mutation_chainA(resA, resE, i)

    id = resA[i].get_id()[1]
    Energy_list[id] = [residue_id(resA[i]), Energy_After_Mutation]


for j in range(len(resE)):

    Energy_After_Mutation = Energy_mutation_chainB(resA, resE, j)

    id = resE[j].get_id()[1]
    Energy_list[id] = [residue_id(resE[j]), Energy_After_Mutation]

energies_for_plotting = []
for res in sorted(Energy_list):
    energies_for_plotting.append(Energy_list[res][1]-E)


# mutated residue , change in energy
data = []
for res in sorted(Energy_list):
    data.append([res, Energy_list[res][0], Energy_list[res]
                [1], Energy_list[res][1]-E])
print(tabulate(data, headers=["residue ID", "Residue Name", "ΔGAB", "ΔΔG"]))

plt.plot(energies_for_plotting, 'o')
plt.ylabel('Difference in Energies')
plt.show()


'''
41 TYR.A41 -75.83751585042594 7.990438536915278
455 LEU.E455 -77.80510576666208 6.022848620679142
456 PHE.E456 -75.28014470977638 8.547809677564842

486 PHE.E486 -68.33113359101328 15.496820796327938
489 TYR.E489 -71.65971714598007 12.168237241361155
505 TYR.E505 -69.69061262636986 14.137341760971367

500 THR.E500 -81.08583971869619 2.74211466864503
501 ASN.E501 -80.2515430162799 3.576411371061326
83 TYR.A83 -79.10323983938142 4.724714547959806
27 THR.A27 -79.70313934172754 4.124815045613687


42 GLN.A42 -87.5094920297493 -3.681537642408074
453 TYR.E453 -89.05260802719378 -5.224653639852562

'''
