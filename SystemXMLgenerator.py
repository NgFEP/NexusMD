from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout


# Load the PDB file
# Load the forcefield
# Create the system

pdb = PDBFile('input.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)


# Serialize the system to an XML file
with open('system.xml', 'w') as outfile:
    xml = XmlSerializer.serialize(system)
    outfile.write(xml)

print("System serialized to system.xml")
