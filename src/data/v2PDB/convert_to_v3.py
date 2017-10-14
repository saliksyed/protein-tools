import glob
import xml.etree.ElementTree as etree
"""
	This script converts all the PDBv2 files in the directory
	in to v3 and copies them in to ../v3PDB
"""


conversion_map = {
	'1HZ': 'HZ1',
	'2HZ': 'HZ2',
	'3HZ': 'HZ3',
	'1HD': 'HD2',
	'2HD': 'HD3',
	'1HB': 'HB2',
	'2HB': 'HB3',
	'1H' : 'H1',
	'2H' : 'H2'
}

def convert_all():
	not_found = {}
	normalized = {}
	for file in glob.glob('*.pdb'):
		with open('../v3PDB/' + file, "w") as output:
			with open(file) as f:
				for line in f.readlines():
					columns = line.split()
					if columns[0] != 'ATOM':
						output.write(line)
					else:
						if columns[2] in conversion_map:
							columns[2] = conversion_map[columns[2]]
						else:
							letter_found = False
							# normalize the names between the PDB files and amber naming convention
							suffix = symbol_normalized = ""
							symbol = columns[2]
							for char in symbol:
								if not char.isdigit():
									letter_found = True
								if letter_found:
									symbol_normalized += char
								else:
									suffix += char
							symbol_normalized += suffix
							columns[2] = symbol_normalized
							normalized[symbol_normalized] = True
							not_found[columns[2]] = True
						output.write("\t".join(columns)+"\n")
	print "SYMBOLS NOT FOUND"
	print not_found.keys()

	print "SYMBOLS REPLACED"
	print normalized.keys()

	print "UNKNOWN SYMBOLS"
	set(normalized.keys()) - set(not_found.keys())


	tree = etree.parse('../amber99sb.xml')
	root = tree.getroot()

	# read in forcefield data:
	field_data = {}
	for child in root:
		for sub_child in child:
			if not sub_child.tag in field_data:
				field_data[sub_child.tag] = []
			field_data[sub_child.tag].append(sub_child)

	atoms = []
	for residue in field_data['Residue']:
		for elem in filter(lambda x : x.tag == 'Atom', residue):
			atoms.append(elem.attrib['name'])

	print "REPLACED SYMBOLS NOT FOUND IN AMBER"
	set(normalized.keys()) - set(atoms).intersection(set(normalized.keys()))


if __name__ == '__main__':
	convert_all()