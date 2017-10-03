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
count = 0
start = time.time()
def sample_conformation():
	global count, start
	conformation = c.get_conformation()
	while True:
		try:
			conformation[random.randint(0, len(conformation)-1)][1] += 0.5
			c.set_conformation(conformation)
			count += 1
			if count%10 == 0:
			    print "%d samples evaluated in %.2f seconds" % (count, (time.time() - start))
			time.sleep(0.1)
		except:
			traceback.print_exc()

t = threading.Thread(target=sample_conformation)
t.start()
r = Visualizer(1920, 1080, c)
