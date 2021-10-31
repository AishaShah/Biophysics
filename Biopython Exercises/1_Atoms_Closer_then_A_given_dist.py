#!/usr/bin/env python
'''
1. Determine the list of pairs of residues whose CA atoms are closer than a given distance

Parameters: PDB file name, distance
'''
import argparse
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch


parser = argparse.ArgumentParser(
    description='list of pairs of residues whose CA atoms are closer than a given distance')
parser.add_argument('pdb_file', help='file name ')
parser.add_argument('cutoff', type=float, help='cut-off distance')
args = parser.parse_args()


print('PDB:', args.pdb_file, 'Cut-off:', args.cutoff)

pdb_parser = PDBParser(PERMISSIVE=True)
st = pdb_parser.get_structure('', args.pdb_file)


select = []
# residues having CA atoms

for atom in st.get_atoms():
    if atom.id == 'CA':
        select.append(atom)

# searching neighbours
nbsearch = NeighborSearch(select)

npair = 0
for at1, at2 in nbsearch.search_all(args.cutoff):
    npair += 1
    print("contact ", npair, ':', '\n', at1.get_parent().get_parent().id, at1.get_parent().get_resname()+'.' + str(at1.get_parent().id[1]), str(at1.get_id())+'.'+str(at1.get_serial_number(
    )), '\n', at2.get_parent().get_parent().id, at2.get_parent().get_resname()+'.'+str(at2.get_parent().id[1]), str(at2.get_id())+'.'+str(at2.get_serial_number()), '\ndistance:', at2-at1, '\n')
