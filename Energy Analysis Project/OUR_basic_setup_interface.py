#!/usr/bin/env python

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import sys
import os
from Bio.PDB.NeighborSearch import NeighborSearch

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from forcefield import VdwParamset

parser = argparse.ArgumentParser(
    prog='structure_setup',
    description='basic structure setup'
)

parser.add_argument(
    '--naccess',
    action='store',
    dest='naccess_bin',
    default=os.path.dirname(os.path.abspath(__file__)) +
    '/soft/NACCESS/naccess',
    help='Vdw parameters'
)
parser.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default=os.path.dirname(os.path.abspath(__file__)) + '/data/vdwprm',
    help='Vdw parameters'
)

parser.add_argument('pdb_file', help='Input PDB', type=open)
parser.add_argument('pdbqt_file', help='Input PDBQT', type=open)

args = parser.parse_args()

print("PDB File:", args.pdb_file.name)
print("PDBQT File:", args.pdbqt_file.name)
print("Parameters File:", args.vdwprm_file)

# loading VdW parameters
ff_params = VdwParamset(args.vdwprm_file)

parser = PDBParser(PERMISSIVE=1)
print('Parsing PDB', args.pdb_file.name)

# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file.name)

# We will use the xtra attribute in Bio.PDB.Atom to hold the new data
# Getting Charges and Atom type from PDBQT
print('Parsing PDBQT', args.pdbqt_file.name)
params = [{}]

# Fix aton numbers when they do not start in 1
i = 1
for at in st.get_atoms():
    at.serial_number = i
    i += 1

for line in args.pdbqt_file:
    line = line.rstrip()
    # Skip TER records from PDBQT
    if line.find('TER') != -1:
        continue
    params.append({'charge': line[69:76], 'type': line[77:].replace(' ', '')})

total_charge = 0.
for at in st.get_atoms():
    at.xtra['atom_type'] = params[at.serial_number]['type']
    at.xtra['charge'] = float(params[at.serial_number]['charge'])
    print(at, at.get_parent())

    at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]
    total_charge += at.xtra['charge']
print('Total Charge: {:8.2f}'.format(total_charge))

# Calculating surfaces
# Srf goes to .xtra['EXP_NACCESS'] field
srf = NACCESS_atomic(st[0], naccess_binary=args.naccess_bin)

# Simple Test for atom fields. Choose one atom from the PDB file. st[model][chain_id][res_number][atom_name]
'''
print(vars(st[0]['A'][42]['N']))
print('')
print(vars(st[0]['A'][42]['N'].xtra['vdw']))

# Getting this is ASA
print((st[0]['A'][42]['N'].xtra['EXP_NACCESS']))
# Getting fsrf that is sigma for solvation
print((st[0]['A'][42]['N'].xtra['vdw']).fsrf)
'''


def solv_energy(st, residues_list):
    solv_energies = {}
    for res in residues_list:
        delG_solv = 0
        for atom in res.get_atoms():
            # calculate solvation of each atom
            # get value of sigma and ASA
            if 'EXP_NACCESS' in atom.xtra:
                delG_solv += atom.xtra['vdw'].fsrf * \
                    float(atom.xtra['EXP_NACCESS'])
        solv_energies[res.id[1]] = delG_solv
    return solv_energies


'''
print('CHAIN A')
print(solv_energy(st, resA))
print('CHAIN E')
print(solv_energy(st, resE))
'''


def solv_energy(st):
    solv_energies = {}
    for res in st.get_residues():
        delG_solv = 0
        for atom in res.get_atoms():
            # calculate solvation of each atom
            # get value of sigma and ASA
            if 'EXP_NACCESS' in atom.xtra:
                delG_solv += atom.xtra['vdw'].fsrf * \
                    float(atom.xtra['EXP_NACCESS'])
        solv_energies[res.id[1]] = delG_solv
    return solv_energies


print(solv_energy(st))
