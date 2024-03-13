from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# Load the PDB file and forcefield
pdb = PDBFile('input.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Create the system
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

# Setup the simulation
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Serialize the system to system.xml
with open('system.xml', 'w') as system_file:
    xml = XmlSerializer.serialize(system)
    system_file.write(xml)

# No simulation steps are executed because simulation.step(n) is never called. 
# Therefore state.xml reflects the initial positions as defined in the PDB file, and has the same step as system.xml
# Serialize the state (including initial positions) to state.xml
state = simulation.context.getState(getPositions=True)
with open('state.xml', 'w') as state_file:
    xml = XmlSerializer.serialize(state)
    state_file.write(xml)

print("System and initial state serialized to system.xml and state.xml respectively.")