protein-tools
============
Tools to generate unfolded protein models from amino-acid sequences

### Features
* Convert amino-acid sequence into unfolded topology
* (TODO) Set backbone conformation using torsion angles

![Screenshot](https://github.com/saliksyed/protein-tools/blob/master/docs/screenshot.png)

### Usage

NOTE: This is an unfinished, personal, experimental project. Master may be
completely broken at any given time! This should *NOT* be used for
any academic/research work. 

Visualizer is currently limited to 15 amino acids in the sequence

```
import sys
import math
import traceback
from topology import Forcefield

f = ForceField()
villin = """MTKLSAQVKGSLNITTPGLQIWRIE"""
c = f.create_chain(villin)
r = Visualizer(Visualizer.RES1080P, c)
r.run()
```


