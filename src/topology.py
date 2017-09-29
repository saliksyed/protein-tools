
from operator import itemgetter
import xml.etree.ElementTree as etree

names = open("data/names.txt", "r").readlines()

class Bond:
	def __init__(self, atom1, atom2, forcefield):
		self.atom1 = atom1
		self.atom2 = atom2
		self.forcefield = forcefield

	def get_atom_1(self):
		return self.atom1

	def get_atom_2(self):
		return self.atom2

class Atom:
	def __init__(self, atom_name, pos, forcefield, bonds=[]): # TODO look up data from amber for atom based on name
		self.name = atom_name
		self.pos = pos
		self.bonds = bonds
		self.forcefield = forcefield



class Topology:
	def __init__(self):
		self.atoms = []
		self.bonds = []

	def add_atom(self, atom):
		self.atoms.append(atom)
		return len(self.atoms) - 1

	def add_bond(self, atom_idx1, atom_idx2):
		b = Bond(self.atoms[atom_idx1], self.atoms[atom_idx2])
		self.bonds.append(b)
		return len(self.bonds) - 1

class Residue:
	def __init__(self, symbol, forcefield):
		self.symbol = symbol
		# TODO: build topology for residue from amber + residue data
		self.forcefield = forcefield
		self.init()

	def init(self):
		self.topology = Topology()
		geometry = self.forcefield.get_geometry_for_residue(self.symbol)
		atoms = self.forcefield.get_atoms_for_residue(self.symbol)
		bonds = self.forcefield.get_bonds_for_residue(self.symbol)
		for atom in geometry:
			symbol = atom[0]
			#a = self.forcefield.get_atom(symbol, pos=atom[-3:])
			#self.topology.add_atom(a)

		# generate the bonds
		print atoms
		print bonds

	def add_bond_to_nitrogen(self, next_res):
		pass
	def add_bond_to_carbon(self, last_res):
		pass

class Chain:
	def __init__(self, forcefield):
		self.residues = []
		self.forcefield = ForceField

	def add_residue(self, residue):
		self.residues.append(residue)
		# todo: link the residues one to the next


class ForceField:
	def __init__(self, data_file="data/amber99sb.xml"):
		pdb_map = {}
		symbol_map = {}
		for name in names:
			row = name.split()
			fname = "_".join(row[2:]).lower() + ".pdb"
			pdb_map[row[0].upper()] = fname
			symbol_map[row[1]] = fname

		residues = {}
		for key in pdb_map:
			fname = pdb_map[key]
			# TODO: read positional data
			contents = map(lambda x : x.split(), open("data/" + fname, "r").readlines()[:-1])
			atoms = map(itemgetter(2), contents) 
			residue = map(itemgetter(3), contents)
			x = map(itemgetter(5), contents)
			y = map(itemgetter(6), contents)
			z = map(itemgetter(7), contents)
			all_data = zip(atoms, residue, x, y, z)
			residues[residue[0]] = all_data

		self.residues = residues
		tree = etree.parse(data_file)
		root = tree.getroot()

		self.field_data = {}
		for child in root:
			for sub_child in child:
				if not sub_child.tag in self.field_data:
					self.field_data[sub_child.tag] = []
				#print sub_child.attrib
				self.field_data[sub_child.tag].append(sub_child)


	def get_atoms_for_residue(self, residue_name):
		for residue in self.field_data['Residue']:
			if residue_name == residue.attrib['name']:
				return filter(lambda x : x.tag == 'Atom', residue)

	def get_bonds_for_residue(self, residue_name):
		for residue in self.field_data['Residue']:
			if residue_name == residue.attrib['name']:
				return filter(lambda x : x.tag == 'Bond', residue)

	def get_geometry_for_residue(self, residue_name):
		if residue_name in self.residues:
			return self.residues[residue_name]
		else:
			return None

	def get_bond_length_params(self, bond):
		return None

	def get_torsion_params(self, bond_1, bond_2, bond_3):
		return None

	def get_angle_params(self, bond):
		return None

	def get_charge_params(self, atom):
		return None
f = ForceField()
a = Residue('CYS', f)