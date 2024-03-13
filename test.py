# test_openmm.py
from openmm import version
print(version.full_version)
from openmm.app import PDBFile
pdb = PDBFile('input.pdb')
print("PDB file loaded successfully.")


import sys, math

# conversion from radians to degrees and vice versa
rad2deg = 180.0 / math.pi
deg2rad = 1.0 / rad2deg


# calculate distance between two 3-d cartesian coordinates
def get_r12(coords1, coords2):
    r2 = 0.0
    for p in range(3):
        r2 += (coords2[p] - coords1[p])**2
    r = math.sqrt(r2)
    return r

# calculate dot product between two unit vectors
def get_udp(uvec1, uvec2):
    udp = 0.0
    for p in range(3):
        udp += uvec1[p] * uvec2[p]
    udp = max(min(udp, 1.0), -1.0)
    return udp

# calculate unit vector between to 3-d cartesian coordinates
def get_u12(coords1, coords2):
    r12 = get_r12(coords1, coords2)
    u12 = [0.0 for p in range(3)]
    for p in range(3):
        u12[p] = (coords2[p] - coords1[p]) / r12
    return u12



# calculate unit cross product between two unit vectors
def get_ucp(uvec1, uvec2):
    ucp = [0.0 for p in range(3)]
    cos_12 = get_udp(uvec1, uvec2)
    sin_12 = math.sqrt(1 - cos_12**2)
    ucp[0] = (uvec1[1]*uvec2[2] - uvec1[2]*uvec2[1]) / sin_12
    ucp[1] = (uvec1[2]*uvec2[0] - uvec1[0]*uvec2[2]) / sin_12
    ucp[2] = (uvec1[0]*uvec2[1] - uvec1[1]*uvec2[0]) / sin_12
    return ucp

# calculate torsion angle between four 3-d cartesian coordinates
def get_t1234(coords1, coords2, coords3, coords4):
    u21 = get_u12(coords2, coords1)
    u23 = get_u12(coords2, coords3)
    u32 = get_u12(coords3, coords2)
    u34 = get_u12(coords3, coords4)
    u21c23 = get_ucp(u21, u23)
    u32c34 = get_ucp(u32, u34)
    dp = get_udp(u21c23, u32c34)
    sign = 2 * float(get_udp(u21c23, u34) < 0) - 1
    t1234 = rad2deg * sign * math.acos(dp)
    return t1234