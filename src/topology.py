
from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
from operator import itemgetter
import xml.etree.ElementTree as etree


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
		self.bonds = bonds
		self.forcefield = forcefield

	def get_position(self):
		return self.pos



class Residue:
	def __init__(self, symbol, forcefield):
		self.symbol = symbol
		self.atoms = []
		self.bonds = []
		self.forcefield = forcefield
		self.child_atom_idx = 0
		self.parent_atom_idx = 0
		self.child_peptide = None

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
		bond_length = 0.1
		self.child_peptide = child

	def render(self):
		# render each of the atoms
		for i, atom in enumerate(self.atoms):
			if atom != None:
				pos = atom.get_position()
				print pos
				glPushMatrix()
				glTranslatef(pos[0], pos[1], pos[2])
				if i == self.child_atom_idx:
					glColor3f(1.0, 0.0, 0.0)
				elif i == self.parent_atom_idx:
					glColor3f(0.0, 0.0, 1.0)
				else:
					glColor3f(0.5, 0.5, 0.5)
				glutSolidSphere(0.1, 20, 20)
				glPopMatrix()

		glLineWidth(3.0)
		glBegin(GL_LINES)
		for bond in self.bonds:
			a1 = bond.get_atom_1()
			a2 = bond.get_atom_2()
			if a1 == None or a2 == None:
				continue
			pos1 = a1.get_position()
			pos2= a2.get_position()
			glColor3f(1.0, 1.0, 1.0)
			glVertex3f(pos1[0], pos1[1], pos1[2])
			glVertex3f(pos2[0], pos2[1], pos2[2])
		glEnd()
		# render the child
		if self.child_peptide:
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
			res_name = self.symbol_to_resname[symbol]
			if i == 0:
				res_name = 'N' + res_name
				root = curr = self._get_residue(res_name)
			else:
				if i == len(sequence) - 1 :
					res_name = 'C' + res_name
				tmp = curr
				curr = self._get_residue(res_name)
				tmp.add_child(curr)
		return root

	def _get_residue(self, residue_name):
		geometry = self._get_geometry_for_residue(residue_name)
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
				atom_geomtery[symbol_normalized] = atom[-3:]

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
			if "N" == residue_atoms[int(bond.attrib["from"])].name:
				ret.set_parent_idx(int(bond.attrib["from"]))
			if "C" == residue_atoms[int(bond.attrib["from"])].name:
				ret.set_child_idx(int(bond.attrib["from"]))

		return ret

	def _get_atoms_for_residue(self, residue_name, include_carbon=True, include_nitrogen=True):
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
