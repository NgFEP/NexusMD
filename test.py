# test_openmm.py
from openmm import version
print(version.full_version)
from openmm.app import PDBFile
pdb = PDBFile('input.pdb')
print("PDB file loaded successfully.")
