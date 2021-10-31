#!/usr/bin/env python
"""
Id 4, but for disulphide bonds. 


Disulphide bonds are formed between S atoms of Cys Residues when they are at the 
appropriate distance (around 1.9 Å). Using the same approach as 1, 3, or 5, find S-S contacts. 
Allow some more distance to access structure variability. 
"""

from Bio.PDB.NeighborSearch import NeighborSearch
import argparse
from Bio.PDB.PDBParser import PDBParser


def residue_id(res):
    return '{} {} {}'.format(res.get_resname(), res.get_parent().id, res.id[1])


def atom_id(at):
    return '{}.{}'.format(residue_id(at.get_parent()), str(at.id))


parser = argparse.ArgumentParser(
    description='list of all disulphide bonds')

parser.add_argument('pdb_file', help='name of pdb file')
parser.add_argument('--cutoff', dest='cutoff', type=float, default=2.5,
                    help='cut-off distance optional, by default it will be 2.5 i.e around 1.9')

args = parser.parse_args()
print('PDB:', args.pdb_file, 'cut off distance:', args.cutoff, '\n')

parser = PDBParser(PERMISSIVE=True)
MAXDIST = args.cutoff
st = parser.get_structure('', args.pdb_file)
res = "CYS"

select = []

for atom in st.get_atoms():
    if 'S' in atom.id:

        if atom.get_parent().get_resname() == res:
            select.append(atom)

if select != []:
    print('list of S atoms of all CYS residues:')
    for atom in select:
        print(atom_id(atom))
    print('\n')

# finding possible S-S bonds
npair = 0
if select != []:
    nbsearch = NeighborSearch(select)

    for at1, at2 in nbsearch.search_all(MAXDIST):
        npair += 1
        print('Contact : ', npair)
        print('Atom 1 : ', atom_id(at1))
        print('Atom 2 :', atom_id(at2))
        print('distance : ', at2-at1)
        print('\n')

if npair == 0:
    print('There are no S atoms around the distance {} in any two Cystine residues so no disulphide bonds exist'.format(MAXDIST))
else:
    print('Number of disulphide bonds are: {}'.format(npair))
