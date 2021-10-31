#!/usr/bin/env python
"""
5. Generate a list of backbone connectivity (i.e. which residues are linked by ordinary peptide bonds). 
Parameters: PDB file name. Optional: Cut-off distance for peptide bonds (defaults to 2.5)
"""
from Bio.PDB import *
import argparse


def residue_id(res):
    # get parent will giv chain
    # return res.get_resname()+' '+str(res.get_parent().id) + str(res.id[1])
    return '{} {} {}'.format(res.get_resname(), res.get_parent().id, res.id[1])


def atom_id(at):
    return '{}.{}'.format(residue_id(at.get_parent()), str(at.id))


parser = argparse.ArgumentParser(
    description='Generates a list of backbone connectivity')

parser.add_argument('pdb_file', help='file name')
parser.add_argument('--cutoff', dest='cutoff', type=float, default=2.5,
                    help='cut-off distance optional, by default it will be 2.5')

args = parser.parse_args()
print('PDB:', args.pdb_file, 'Cut-off:', args.cutoff, '\n')
pdb_parser = PDBParser()

st = pdb_parser.get_structure('', args.pdb_file)
MAXDIST = args.cutoff
select = []
req_atoms = ['N', 'C']

# Collecting N and C atoms

for atom in st.get_atoms():
    if atom.id in req_atoms:
        select.append(atom)

# searching neighbours
# if distance less than given cutoff (by default 2.5) then two atoms are involved in backbone connectivity

nbsearch = NeighborSearch(select)
npair = 0
print('list of backbone connectivity')
for at1, at2 in nbsearch.search_all(MAXDIST):

    if at1.get_parent().id[1] != at2.get_parent().id[1]:
        npair += 1
        print("Peptide Bond", npair, ':\n', atom_id(at1), '\n', atom_id(
            at2), '\n', 'Distance : ', at2-at1, '\n')
