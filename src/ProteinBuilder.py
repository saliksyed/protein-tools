from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from operator import itemgetter


def load_geometry():
	names = open("data/names.txt", "r").readlines()
	residues = {}
	for name in names:
		row = name.split()
		fname = "_".join(row[2:]).lower() + ".pdb"
		contents = filter(lambda x : x[0]=='ATOM', map(lambda x : x.split(), open("data/v3PDB/" + fname, "r").readlines()[:-1]))
		atoms = map(itemgetter(2), contents) 
		residue = map(itemgetter(3), contents)
		x = map(float, map(itemgetter(5), contents))
		y = map(float, map(itemgetter(6), contents))
		z = map(float, map(itemgetter(7), contents))
		all_data = zip(atoms, residue, x, y, z)
		atom_geometry = {}
		for atom in all_data:
			atom_geometry[atom[0]] = list(atom[-3:])
		residues[residue[0]] = atom_geometry
	return residues

def transform_geometry(atoms, transform):
	return [atoms[i] + transform[i] for i in  xrange(0, 3)]

# TODO: lookup peptide bond length from forcefield instead of hard coding
def build_topology(seq, forcefield, peptide_bond_length=1.43):
	geometry = load_geometry()
	names_map = dict(map(lambda x : (x.split()[1].upper(), x.split()[0].upper()) , open("data/names.txt").readlines()))
	t = Topology()
	c = t.addChain('c1')
	idx_offset = 0
	atoms = []
	positions = []
	transformStack = [[0, 0, 0]]
	stackOffsets = []
	last_carbon_idx = None
	# keey track of backbone atoms
	for i, symbol in enumerate(seq):
		residueSymbol = names_map[symbol]
		geometry_for_residue = geometry[residueSymbol]
		
		if residueSymbol == "HIS":
			residueSymbol = "HIE"
		r = t.addResidue(residueSymbol, c)
		
		if i == 0:
			residueSymbol = "N" + residueSymbol
		if i == len(seq) - 1:
			residueSymbol = "C" + residueSymbol

		for atom in forcefield._templates[residueSymbol].atoms:
			atoms.append(t.addAtom(atom.name, forcefield._atomTypes[atom.type].element, r))
			if atom.name in geometry_for_residue:
				positions.append(geometry_for_residue[atom.name])
			else:
				#print "no geometry for:"
				#print atom.name
				positions.append([0, 0, 0])
		
		for bond in forcefield._templates[residueSymbol].bonds:
			t.addBond(atoms[bond[0] + idx_offset], atoms[bond[1] + idx_offset])


		curr_nitrogen_idx = None
		curr_carbon_idx = None

		for bond in forcefield._templates[residueSymbol].externalBonds:
			a = atoms[bond + idx_offset]
			if a.element.name == "carbon":
				curr_carbon_idx = bond + idx_offset
			elif a.element.name =="nitrogen":
				curr_nitrogen_idx = bond + idx_offset

		if curr_nitrogen_idx:
			loc_to = positions[curr_nitrogen_idx]
		else:
			loc_to = None

		if last_carbon_idx:
			loc_from = positions[last_carbon_idx]
		else:
			loc_from = [0, 0, 0]

		if loc_to != None and loc_from != None:
			curr_t = [loc_from[i] - loc_to[i] for i in  xrange(0, 3)]
			curr_t[0] -= peptide_bond_length
			transformStack.append(curr_t)

		# add bond for backchain
		if last_carbon_idx != None and curr_nitrogen_idx !=None:
			t.addBond(atoms[last_carbon_idx], atoms[curr_nitrogen_idx])
		
		if curr_carbon_idx:
			last_carbon_idx = curr_carbon_idx
		idx_offset += len(forcefield._templates[residueSymbol].atoms)
		stackOffsets.append(len(forcefield._templates[residueSymbol].atoms))
	
	transformed_positions = []
	idx = 0
	currStackIdx = 0
	transformStack.append([0, 0, 0]) # No transform for end of chain!
	curr_t = transformStack[0]
	for offset in stackOffsets:
		for i in xrange(0, offset):
			transformed_positions.append(transform_geometry(positions[idx], curr_t))
			idx += 1
		currStackIdx += 1
		curr_t = [curr_t[i] + transformStack[currStackIdx][i] for i in  xrange(0, 3)]
	return t, transformed_positions


