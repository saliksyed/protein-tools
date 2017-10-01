import sys
import math
import traceback
from topology import *
from visualizer import *

f = ForceField()
villin = """MTKLSAQVKGSLNITTPGLQIWRIE"""
c = f.create_chain(villin)
r = Visualizer(Visualizer.RES1080P, c)
r.run()
