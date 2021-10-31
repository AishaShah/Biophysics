#!/usr/bin/env python
"""
Determine all possible hydrogen bonds (Polar atoms at less than 3.5 Å).  
Parameters: PDB file name. Optional: cut-off distance (defaults to 3.5) 
"""

from Bio.PDB import *
import argparse
from Bio.PDB.NeighborSearch import NeighborSearch


def residue_id(res):
    return '{} {} {}'.format(res.get_resname(), res.get_parent().id, res.id[1])


def atom_id(at):
    return '{}.{}'.format(residue_id(at.get_parent()), str(at.id))


# make instance of obj doing something lieke parser
parser = argparse.ArgumentParser(description='get possible H-bonds')

parser.add_argument('pdb_file', help='file name')
parser.add_argument('--cutoff', dest='cutoff', type=float, default=3.5,
                    help='cut-off distance optional, by default it will be 3.5')

args = parser.parse_args()  # execute parser


print('PDB:', args.pdb_file, 'Cut-off:', args.cutoff, '\n')

pdb_parser = PDBParser(PERMISSIVE=True)

st = pdb_parser.get_structure('', args.pdb_file)
MAXDIST = args.cutoff
select = []
polar = ['O', 'N', 'S']
for atom in st.get_atoms():
    if atom.id[0] in polar:
        select.append(atom)


nbsearch = NeighborSearch(select)
npair = 0
for at1, at2 in nbsearch.search_all(MAXDIST):
    npair += 1
    print("H_bond ", npair, ':\n', atom_id(at1), '\n',
          atom_id(at2), '\nDistance : ', at2-at1, '\n')
