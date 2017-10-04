protein-tools
============
Tools to generate unfolded protein models from amino-acid sequences

### Features
* Convert amino-acid sequence into unfolded topology
* Set backbone conformation using torsion angles / bond-lengths
* Export topology to PDB v3.0 format

### Dependencies
* [Numpy](http://www.numpy.org/)
* [Pyrr matrix library](https://github.com/adamlwgriffiths/Pyrr)
* [PyOpenGL](https://github.com/mcfletch/pyopengl) (for visualizer)
* [lxml](http://lxml.de/)

### Usage

NOTE: This is an personal, experimental project. You should check the results prior to using this for academic work!

```
import sys
import math
import traceback
from topology import *
from visualizer import *
import threading
import time

f = ForceField()
villinN68H = """LSDEDFKAVFGMTRSAFANLPLWKQQHLKKEKGLF"""
c = f.create_chain(villinN68H)
c.write_to_pdb("data/output.pdb")

# Sample Conformations:
count = 0
start = time.time()
def sample_conformation():
	global count, start
	while True:
	    c.set_conformation(c.get_random_conformation())
	    count += 1
	    if count%10 == 0:
	        print "%d samples evaluated in %.2f seconds" % (count, (time.time() - start))
	    time.sleep(0.1)

t = threading.Thread(target=sample_conformation)
t.start()
r = Visualizer(1920, 1080, c)

```


