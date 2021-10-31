#!/usr/bin/env python
"""
Generate a list of all CA atoms of given residue type with coordinates 
Parameters: PDB file name, residue type. 
"""


import argparse
from Bio.PDB.PDBParser import PDBParser


def residue_id(res):
    return '{} {} {}'.format(res.get_resname(), res.get_parent().id, res.id[1])


def atom_id(at):
    return '{}.{}'.format(residue_id(at.get_parent()), str(at.id))


parser = argparse.ArgumentParser(
    description='list of all CA atoms of given residue name with coordinates')

parser.add_argument('pdb_file', help='file name')
parser.add_argument('residue', help='Residue name')
parser.add_argument('--chain', dest='chain',
                    help='Chain id (optional)', default='*')

args = parser.parse_args()
print('PDB:', args.pdb_file, 'Name of residue:', args.residue, '\n')

parser = PDBParser()

st = parser.get_structure('', args.pdb_file)
res = args.residue
chain = args.chain
selec = []


for at in st.get_atoms():
    if at.id == 'CA':
        if at.get_parent().get_resname() == res:
            if chain == '*' or chain == at.get_parent().get_parent().id:
                selec.append(at)
if selec != []:
    print('list of all CA atoms of all available {} residues with coordinates : '.format(res))
    for atom in selec:
        print(atom_id(atom), atom.get_coord())

else:
    print('Residue provided is not available')
