import xml.etree.ElementTree as etree
import numpy as np
from operator import itemgetter
from pyrr import Matrix44, Vector3, Vector4
import pyrr
import random

class SYMBOL:
	NITROGEN = "N"
	CARBON = "C"
	ALA = "ALA"
	ARG = "ARG"
	ASN = "ASN"
	ASP = "ASP"
	CYS = "CYS"
	GLN = "GLN"
	GLU = "GLU"
	GLY = "GLY"
	HIS = "HIS"
	HIS_DEFAULT_ISOMER = "HIE"
	ILE = "ILE"
	LEU = "LEU"
	LYS = "LYS"
	MET = "MET"
	PRO = "PRO"
	PHE = "PHE"
	SER = "SER"
	THR = "THR"
	TRP = "TRP"
	TYR = "TYR"
	VAL = "VAL"

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
	def __init__(self, atom_name, pos, forcefield, bonds=[]):
		self.name = atom_name
		self.pos = pos
		self.transformed_pos = Vector3(pos)
		self.bonds = bonds
		self.forcefield = forcefield

	def get_position(self):
		return self.pos

	def get_transformed_position(self):
		return self.transformed_pos

	def apply_transform(self, mat):
		self.transformed_pos = pyrr.matrix44.apply_to_vector(mat, Vector4(self.get_position() + [1.0]))

class Residue:
	def __init__(self, symbol, forcefield):
		self.symbol = symbol
		self.atoms = []
		self.bonds = []
		self.forcefield = forcefield
		self.child_atom_idx = None
		self.parent_atom_idx = None
		self.child_peptide = None
		self.conformation = None
		self.color = [random.random(), random.random(), random.random()]

	def get_default_conformation(self):
		bond_length = 1.47
		torsion_angle = 0.0
		ret = [[bond_length, torsion_angle]]
		if self.child_peptide:
			return ret + self.child_peptide.get_default_conformation()
		else:
			return ret 

	def get_random_conformation(self):
		bond_length = 1.47
		torsion_angle = random.random() * 360.0
		ret = [[bond_length, torsion_angle]]
		if self.child_peptide:
			return ret + self.child_peptide.get_random_conformation()
		else:
			return ret 

	def add_atom(self, atom):
		self.atoms.append(atom)
		return len(self.atoms) - 1

	def add_bond(self, atom_idx1, atom_idx2):
		b = Bond(self.atoms[atom_idx1], self.atoms[atom_idx2], self.forcefield)
		self.bonds.append(b)
		return len(self.bonds) - 1

	def set_child_idx(self, idx):
		self.child_atom_idx = idx

	def set_parent_idx(self, idx):
		self.parent_atom_idx = idx
	

	def get_child_atom(self):
		if self.child_atom_idx == None:
			return None
		return self.atoms[self.child_atom_idx]

	def get_parent_atom(self):
		if self.parent_atom_idx == None:
			return None
		return self.atoms[self.parent_atom_idx]

	def add_child(self, child):
		# TODO: Is this right? seems like all backbone bonds should be 180 deg bond angle. Looks correct visually
		self.bond_axis = -1.0 * Vector3([1,0,0])

		loc_from = self.atoms[self.parent_atom_idx].get_position()
		loc_to = child.atoms[child.child_atom_idx].get_position()
		t = [loc_from[i] - loc_to[i] for i in  xrange(0, 3)]

		self.bond_transform = Vector3(t)
		self.child_peptide = child

	def get_conformation(self):
		return self.conformation

	def get_bond_length(self):
		return self.conformation[self.conformation_idx][0]

	def get_torsion_angle(self):
		return self.conformation[self.conformation_idx][1]


	def get_pdb_str(self, atom_idx, peptide_idx):
		ret = ""
		count = 0
		for i, atom in enumerate(filter(lambda x : x!=None, self.atoms)):
			line = []
			line.append(str(atom_idx + i + 1))
			line.append(atom.name)
			line.append(self.symbol)
			line.append(str(peptide_idx + 1))
			pt = atom.get_transformed_position()
			line.append(str(pt[0]))
			line.append(str(pt[1]))
			line.append(str(pt[2]))
			line.append("1.00")
			line.append("0.00")
			ret += "ATOM" + "".join(["%10s"]*len(line)) % tuple(line) + "\n"

			count += 1
		return ret, count

	def write_to_pdb(self, file):
		with open(file, "w") as f:
			f.write("REMARK Generated by Protein-tools https://github.com/saliksyed/protein-tools\n")
			c = self
			i = 0
			atoms_seen = 0
			while c != None:
				lines, newly_seen = c.get_pdb_str(atoms_seen, i)
				f.write(lines)
				atoms_seen += newly_seen
				c = c.child_peptide
				i += 1


	def set_conformation(self, conformation, mat=None, idx=0):
		self.conformation = conformation
		self.conformation_idx = idx
		bond_length, torsion_angle = conformation[idx]
		if (type(mat) != type(None)):
			for atom in self.atoms:
				if atom:
					atom.apply_transform(mat)
			
		if self.child_peptide:
			if (type(mat) == type(None)):
				mat = Matrix44.identity()

			transformed_axis = Vector3(pyrr.vector.normalise(pyrr.matrix44.apply_to_vector(mat, self.bond_axis)))
			rmat = pyrr.matrix44.create_from_axis_rotation(transformed_axis, torsion_angle)
			tmat = pyrr.matrix44.create_from_translation(self.bond_transform + bond_length * self.bond_axis)
			curr_t = pyrr.matrix44.multiply(tmat, rmat)
			curr_t = pyrr.matrix44.multiply(curr_t, mat)
			self.child_peptide.set_conformation(conformation, curr_t, idx + 1)


class ForceField:
	def __init__(self, data_file="data/amber99sb.xml"):
		pdb_map = {}
		names = open("data/names.txt", "r").readlines()
		self.symbol_to_resname = {}
		for name in names:
			row = name.split()
			fname = "_".join(row[2:]).lower() + ".pdb"
			pdb_map[row[0].upper()] = fname
			self.symbol_to_resname[row[1]] = row[0].upper()
		# Read in geometric information for each residue from PDB:
		residues = {}
		for key in pdb_map:
			fname = pdb_map[key]
			contents = filter(lambda x : x[0]=='ATOM', map(lambda x : x.split(), open("data/v3PDB/" + fname, "r").readlines()[:-1]))
			atoms = map(itemgetter(2), contents) 
			residue = map(itemgetter(3), contents)
			x = map(float, map(itemgetter(5), contents))
			y = map(float, map(itemgetter(6), contents))
			z = map(float, map(itemgetter(7), contents))
			all_data = zip(atoms, residue, x, y, z)
			residues[residue[0]] = all_data
		self.residues = residues
		tree = etree.parse(data_file)
		root = tree.getroot()

		# read in forcefield data:
		self.field_data = {}
		for child in root:
			for sub_child in child:
				if not sub_child.tag in self.field_data:
					self.field_data[sub_child.tag] = []
				self.field_data[sub_child.tag].append(sub_child)

	def get_residue_symbols(self):
		return self.symbol_to_resname.keys()

	def create_chain(self, sequence):
		root = None
		curr = None
		color = 1.0 / len(sequence)
		for i, symbol in enumerate(sequence):
			if not symbol.rstrip():
				continue
			if root == None:
				curr = root = self.get_residue(symbol, i==0, i==len(sequence) - 1)
			else:
				tmp = curr
				curr = self.get_residue(symbol, i==0, i==len(sequence) - 1)
				tmp.add_child(curr)
		root.set_conformation(root.get_default_conformation())
		return root

	def _symbol_to_residue_name(self, symbol, is_start_of_chain=False, is_end_of_chain=False):
		res_name = self.symbol_to_resname[symbol]
		if res_name == SYMBOL.HIS:
			res_name = SYMBOL.HIS_DEFAULT_ISOMER
		if is_start_of_chain:
			return SYMBOL.CARBON + res_name
		elif is_end_of_chain:
			return SYMBOL.NITROGEN + res_name
		else:
			return res_name

	def get_residue(self, symbol, is_start_of_chain=False, is_end_of_chain=False):
		raw_residue_name = self.symbol_to_resname[symbol]
		forcefield_residue_name = self._symbol_to_residue_name(symbol, is_start_of_chain, is_end_of_chain)
		geometry = self._get_geometry_for_symbol(symbol)
		atoms = self._get_atoms_for_residue(forcefield_residue_name)
		bonds = self._get_bonds_for_residue(forcefield_residue_name)
		external_bonds = self._get_external_bonds_for_residue(forcefield_residue_name)

		atom_names = map(lambda x : x.attrib["name"], atoms)
		residue_atoms = []
		atom_geometry = {}

		for atom in geometry:
			atom_geometry[atom[0]] = list(atom[-3:])
			
		ret = Residue(raw_residue_name, self)

		for atom in atoms:
			name = atom.attrib["name"]
			if name in atom_geometry:
				residue_atoms.append(Atom(name, atom_geometry[name], self))
			else:
				residue_atoms.append(None)

		for atom in residue_atoms:
			ret.add_atom(atom)

		for bond in bonds:
			ret.add_bond(int(bond.attrib["from"]), int(bond.attrib["to"]))

		for bond in external_bonds:
			# TODO: make bonding order configurable?
			bond_from_name = residue_atoms[int(bond.attrib["from"])].name
			if bond_from_name == SYMBOL.NITROGEN:
				ret.set_parent_idx(int(bond.attrib["from"]))
			if bond_from_name == SYMBOL.CARBON:
				ret.set_child_idx(int(bond.attrib["from"]))
		return ret

	def _get_atoms_for_residue(self, residue_name):
		atoms = []
		for residue in self.field_data['Residue']:
			if residue_name == residue.attrib['name']:
				for elem in filter(lambda x : x.tag == 'Atom', residue):
					atoms.append(elem)
		return atoms

	def _get_bonds_for_residue(self, residue_name):
		for residue in self.field_data['Residue']:
			if residue_name == residue.attrib['name']:
				return filter(lambda x : x.tag == 'Bond', residue)

	def _get_external_bonds_for_residue(self, residue_name):
		for residue in self.field_data['Residue']:
			if residue_name == residue.attrib['name']:
				return filter(lambda x : x.tag == 'ExternalBond', residue)

	def _get_geometry_for_symbol(self, symbol):
		residue_name = self.symbol_to_resname[symbol]
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
