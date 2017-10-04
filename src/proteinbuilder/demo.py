import sys
import math
import traceback
from topology import *
from visualizer import *
import threading
import time
import random

f = ForceField()
villinN68H = """LSDEDFKAVFGMTRSAFANLPLWKQQHLKKEKGLF"""
c = f.create_chain(villinN68H)
c.write_to_pdb("data/output.pdb")

# Sample Conformations:
def sample_conformation():
	global count, start
	conformation = c.get_conformation()
	while True:
		conformation[random.randint(0, len(conformation)-1)][1] += 0.5
		c.set_conformation(conformation)
		print f.get_energy(c)
		time.sleep(0.01)

t = threading.Thread(target=sample_conformation)
t.start()
r = Visualizer(1920, 1080, c)
