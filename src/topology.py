import xml.etree.ElementTree as etree
import numpy as np
from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
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
		self.child_atom_idx = 0
		self.parent_atom_idx = 0
		self.child_peptide = None
		self.bond_axis = None
		self.bond_length = 1.0
		self.color = [random.random(), random.random(), random.random()]

	def get_default_conformation(self):
		bond_length = 1.0
		torsion_angle = 0.0
		ret = [(bond_length, torsion_angle)]
		if self.child_peptide:
			return ret + self.child_peptide.get_default_conformation()
		else:
			return ret 

	def get_random_conformation(self):
		bond_length = 4.5
		torsion_angle = 0.0 #random.random() * 360.0
		ret = [(bond_length, torsion_angle)]
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

	def set_side_chain_geometry(self, geometry):
		pass

	def add_child(self, child):
		loc_from = self.atoms[self.child_atom_idx].get_position()
		loc_to = child.atoms[child.parent_atom_idx].get_position()
		t = [loc_from[i] - loc_to[i] for i in  xrange(0, 3)]
		self.bond_axis = Vector3(pyrr.vector.normalise(Vector3(t)))
		self.child_peptide = child

	def set_conformation(self, conformation, mat=None):
		bond_length, torsion_angle = conformation.pop(0)
		if (type(mat) != type(None)):
			for atom in self.atoms:
				if atom:
					atom.apply_transform(mat)
			
		if self.child_peptide:
			rmat = pyrr.matrix44.create_from_axis_rotation(self.bond_axis, torsion_angle)
			tmat = pyrr.matrix44.create_from_translation(bond_length * self.bond_axis)
			curr_t = pyrr.matrix44.multiply(tmat, rmat)
			if (type(mat) != type(None)):
				curr_t = pyrr.matrix44.multiply(mat, curr_t)
			self.child_peptide.set_conformation(conformation, curr_t)

	def render(self):
		# render each of the atoms
		for i, atom in enumerate(self.atoms):
			if atom != None:
				r = 0.1
				pos = atom.get_transformed_position()
				glPushMatrix()
				glTranslatef(pos[0], pos[1], pos[2])
				if i == self.child_atom_idx:
					glColor3f(self.color[0], self.color[1], self.color[2])
					r = 0.2
				elif i == self.parent_atom_idx:
					glColor3f(self.color[0], self.color[1], self.color[2])
					r = 0.2
				else:
					glColor3f(0.5, 0.5, 0.5)
				glutSolidSphere(r, 20, 20)
				glPopMatrix()

		glLineWidth(3.0)
		glBegin(GL_LINES)
		for bond in self.bonds:
			a1 = bond.get_atom_1()
			a2 = bond.get_atom_2()
			if a1 == None or a2 == None:
				continue
			pos1 = a1.get_transformed_position()
			pos2= a2.get_transformed_position()
			glColor3f(self.color[0], self.color[1], self.color[2])
			glVertex3f(pos1[0], pos1[1], pos1[2])
			glVertex3f(pos2[0], pos2[1], pos2[2])
		glEnd()
		# render the child
		if self.child_peptide:
			glLineWidth(8.0)
			glBegin(GL_LINES)
			glColor3f(1.0, 0.0, 0.0)
			loc_from = self.atoms[self.child_atom_idx].get_transformed_position()
			loc_to = self.child_peptide.atoms[self.child_peptide.parent_atom_idx].get_transformed_position()
			glVertex3f(loc_from[0], loc_from[1], loc_from[2])
			glColor3f(0.0, 0.0, 1.0)
			glVertex3f(loc_to[0], loc_to[1], loc_to[2])
			glEnd()
			self.child_peptide.render()

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
			contents = filter(lambda x : x[0]=='ATOM', map(lambda x : x.split(), open("data/" + fname, "r").readlines()[:-1]))
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

	def create_chain(self, sequence):
		root = None
		curr = None
		for i, symbol in enumerate(sequence):
			if not symbol.rstrip():
				continue
			res_name = self.symbol_to_resname[symbol]
			if i == 0:
				res_name = SYMBOL.NITROGEN + res_name
				root = curr = self._get_residue(res_name)
			else:
				if i == len(sequence) - 1 :
					res_name = SYMBOL.CARBON + res_name
				tmp = curr
				curr = self._get_residue(res_name)
				tmp.add_child(curr)
		return root

	def _get_residue(self, residue_name):
		geometry = self._get_geometry_for_residue(residue_name)

		# TODO: Replace HIS more intelligently?
		if residue_name == SYMBOL.HIS:
			residue_name = SYMBOL.HIS_DEFAULT_ISOMER

		atoms = self._get_atoms_for_residue(residue_name)
		bonds = self._get_bonds_for_residue(residue_name)
		external_bonds = self._get_external_bonds_for_residue(residue_name)

		atom_names = map(lambda x : x.attrib["name"], atoms)
		residue_atoms = []
		atom_geomtery = {}
		for atom in geometry:
			symbol = atom[0]
			letters = filter(lambda x : not x.isdigit(), symbol)
			numbers = filter(lambda x : x.isdigit(), symbol)
			symbol_normalized = "".join(letters + numbers)
			if symbol_normalized in atom_names:
				atom_geomtery[symbol_normalized] = list(atom[-3:])

		ret = Residue(residue_name, self)
		for atom in atoms:
			name = atom.attrib["name"]
			if name in atom_geomtery:
				residue_atoms.append(Atom(name, atom_geomtery[name], self))
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

	def _get_geometry_for_residue(self, residue_name):
		if len(residue_name) == 4:
			residue_name = residue_name[1:]
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
