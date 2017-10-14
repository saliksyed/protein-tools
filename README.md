protein-tools
============
Tools to generate unfolded protein models from amino-acid sequences for OpenMM

### Features
* Convert amino-acid sequence into unfolded topology
* Export topology to PDB v3.0 format

### Dependencies
* [OpenMM](https://github.com/pandegroup/openmm)

### Usage

NOTE: This is an personal, experimental project. You should check the results prior to using this for academic work!

```
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from operator import itemgetter
from ProteinBuilder import build_topology
from sys import stdout

forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
villinN68H = """LSDEDFKAVFGMTRSAFANLPLWKQQHLKKEKGLF"""
topology, positions = build_topology(villinN68H, forcefield)

topology.setUnitCellDimensions([6.0, 6.0, 6.0])
print "Creating System"
system = forcefield.createSystem(topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds, ignoreExternalBonds=True)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
print "Stepping system"
simulation.step(10000)
```


