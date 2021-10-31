"""
Print distances between all atom pairs of two given residues
Parameters: PDB file name, Residue 1, Residue 2
"""

import sys
import argparse
from Bio.PDB import *


def residue_id(res):
    return '{} {} {}'.format(res.get_resname(), res.get_parent().id, res.id[1])


def atom_id(at):
    return '{}.{}'.format(residue_id(at.get_parent()), str(at.id))


parser = argparse.ArgumentParser(
    description='Print distances between all atom pairs of two given residues ')

parser.add_argument('pdb_file', help='file name in pdb format')
parser.add_argument(
    'residue_1', help='chain and id of residue 1 in format [Chain]Res_id')
parser.add_argument(
    'residue_2', help='chain and id of residue 2 in format [Chain]Res_id')
args = parser.parse_args()  # execute parser

chain1 = args.residue_1[0]
res1 = int(args.residue_1[1:])

chain2 = args.residue_2[0]
res2 = int(args.residue_2[1:])


print('PDB:', args.pdb_file, '\nChain_1', chain1, ',',
      'res 1:', res1, '\nChain_2', chain2, ',', 'res 2:', res2,)


parser = PDBParser()
st = parser.get_structure('', args.pdb_file)

Available_chains = [ch.id for ch in st[0]]
if (chain1 and chain2) not in Available_chains:
    print('Chains in this structure are : ', *Available_chains)
    print('Please provide the correct chains')
    sys.exit()

try:
    r1 = st[0][chain1][res1]
except KeyError:
    print('Residue Number {} not available in chain {}'.format(res1, chain1))
    sys.exit()
try:
    r2 = st[0][chain2][res2]
except KeyError:
    print('Residue Number {} not available in chain {}'.format(res2, chain2))
    sys.exit()

print('\nResidue {} of Chain {} is {} '.format(res1, chain1, r1.get_resname()))
print('Residue {} of Chain {} is {} \n'.format(res2, chain2, r2.get_resname()))

print('Distance between all atom pairs of given two residues is as follows: \n')
for at1 in r1.get_atoms():
    for at2 in r2.get_atoms():
        dist = at2-at1
        print('Atom 1 : ', atom_id(at1))
        print('Atom 2 : ', atom_id(at2))
        print('Distance : ', dist)
        print('\n')
