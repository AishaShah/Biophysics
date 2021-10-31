#!/usr/bin/env python
"""
Generate a list of all atoms for a given residue number
Parameters: PDB file name, Residue number (Including Chain if applicable)
"""
import sys
from Bio.PDB import *

import argparse
from Bio.PDB.NeighborSearch import NeighborSearch


def residue_id(res):
    return '{} {} {}'.format(res.get_resname(), res.get_parent().id, res.id[1])


def atom_id(at):
    return '{}.{}'.format(residue_id(at.get_parent()), str(at.id))


parser = argparse.ArgumentParser(
    description='list of all atoms for a given residue number, chain(optional),model number(optional)')
parser.add_argument('pdb_file', help='pdb file name')
parser.add_argument(
    'residue', help='Chain name(optional) and Residue number in format [chain]res_num')
parser.add_argument('--model_num', dest='model_num', type=int, default=0,
                    help='model number , by default its 0')

args = parser.parse_args()

print('PDB:', args.pdb_file, 'Residue:', args.residue, '\n')
if args.residue[0].isalpha():
    chain = args.residue[0]
    res_num = int(args.residue[1:])
else:
    chain = '*'  # all chains
    res_num = int(args.residue)


pdb_parser = PDBParser()

st = pdb_parser.get_structure('', args.pdb_file)
model = args.model_num
# if provided optional model which does not exist:
if model >= len(st):
    print("Id's of the models available in given file are :", end=' ')
    for i in range(len(st)):
        print(i, end=' ')
    print('\n')
    sys.exit()

for ch in st[model]:
    if ch.id == chain or chain == '*':
        for res in ch.get_residues():
            if res.id[1] == res_num:
                print("Atoms with their coordinates of Residue", res.get_resname(),
                      res_num, 'of Chain', ch.id, 'are:\n')
                for atom in res.get_atoms():
                    print(atom_id(atom), atom.get_serial_number(),
                          atom.get_coord())
    print('\n')
