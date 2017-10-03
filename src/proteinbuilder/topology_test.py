import unittest
from topology import *
from pyrr import Matrix44, Vector3, Vector4

class TestProteinBuilder(unittest.TestCase):

	def test_all_internal_residue_link_atoms(self):
		"""	
			The child atom is the atom of this amino acid that connects
			to the parent amino acid (the predecessor in the chain).
			In the example below, the child atoms are in [] while parent atoms are in {}:
			[C]==C=={N}==>[C]==C=={N}==>[C]==C=={N}==> 
		"""

		f = ForceField()
		symbols = f.get_residue_symbols()
		for symbol in symbols:
			c = f.get_residue(symbol)
			self.assertEqual(c.get_parent_atom().name, SYMBOL.NITROGEN)
			self.assertEqual(c.get_child_atom().name, SYMBOL.CARBON)


			# Check that if it is the beginning of the chain no child exists
			c = f.get_residue(symbol, True)

			self.assertEqual(c.get_parent_atom().name, SYMBOL.NITROGEN)
			self.assertEqual(c.get_child_atom(), None)

			# Check that if it is the end of the chain no child exists
			c = f.get_residue(symbol, False, True)
			self.assertEqual(c.get_child_atom().name, SYMBOL.CARBON)
			self.assertEqual(c.get_parent_atom(), None)

	def no_overlapping_atoms(self, c):
		for i in xrange(0, len(c.atoms)):
			for j in xrange(i + 1, len(c.atoms)):
				d = pyrr.vector.length(Vector3(c.atoms[i].get_position()) - Vector3(c.atoms[j].get_position()))
				if d < 0.0001:
					print "Overlap found"
					print c.symbol
					print c.atoms[i].name
					print c.atoms[j].name
					print "=========="
					return i, j
		return True


	def test_geometry_is_unique(self):
		"""
			No two atoms should be assigned the same geometry
		"""
		f = ForceField()
		symbols = f.get_residue_symbols()
		for symbol in symbols:
			c = f.get_residue(symbol)
			self.assertEquals(self.no_overlapping_atoms(c), True)
			c = f.get_residue(symbol, True)
			self.assertEquals(self.no_overlapping_atoms(c), True)
			c = f.get_residue(symbol, False, True)
			self.assertEquals(self.no_overlapping_atoms(c), True)
		return True

if __name__ == '__main__':
    unittest.main()





