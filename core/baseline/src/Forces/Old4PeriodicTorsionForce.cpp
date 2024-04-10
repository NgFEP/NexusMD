#include "PeriodicTorsionForce.h"
#include "Coords3D.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "SystemXMLParser.h" //this contains definitions for TorsionParameters

#ifndef PI
#define PI 3.14159265358979323846
#endif

const double tenm3 = std::numeric_limits<double>::epsilon(); // Smallest difference between two doubles


using namespace std;
using namespace BaseLine; 

PeriodicTorsionForce::PeriodicTorsionForce() {
}

PeriodicTorsionForce::~PeriodicTorsionForce() {
}

// Calculate the distance between two points
double PeriodicTorsionForce::atom_distance(const Coords3D& coords1, const Coords3D& coords2) {
    Coords3D diff = coords2 - coords1;
    return sqrt(diff.dot(diff));
}


// Calculate torsion angle (in radians) between four points
double PeriodicTorsionForce::torsion_angle(const Coords3D& pi, const Coords3D& pj, const Coords3D& pk, const Coords3D& pl) {
    
    //3 vectors between four atoms of i, j, k, l 
    Coords3D rij = pi - pj;
    Coords3D rkj = pk - pj;
    Coords3D rkl = pk - pl;

    Coords3D n1 = rij.cross(rkj);//normal vector of the ijk plane
    Coords3D n2 = rkj.cross(rkl);//normal vector of the jkl plane

    //start of previous method
    Coords3D m1 = n1.cross(rkj).normalize();

    double x = n1.dot(n2);
    double y = m1.dot(n2);
    double angle2 = atan2(y, x);
    //end of previous method

    double fxi = n1.dot(n1);
    double fyi = n2.dot(n2);
    double ct = n1.dot(n2);

    double z1, z2, z12, fzi;
    // Branch if linear dihedral:
    if (tenm3 <= fxi) {
        z1 = 1.0 / fxi;
    }
    else {
        z1 = 0.0;
    }

    if (tenm3 <= fyi) {
        z2 = 1.0 / fyi;
    }
    else {
        z2 = 0.0;
    }

    z12 = z1 * z2;

    if (z12 != 0.0) {
        fzi = 1.0;
    }
    else {
        fzi = 0.0;
    }

    double s = rkj.dot(n1.cross(n2));//   xkj* (dz.cross(gy) - dy.cross(gz)) + ykj * (dx.cross(gz) - dz.cross(gx)) + zkj * (dy.cross(gx) - dx.cross(gy));

    double angle = PI - copysign(acos(max(-1.0, min(1.0, ct * z12))), s); // Assuming s.dot(Coords3D(0,0,1)) simulates the sign function.

    return angle;
}

// Calculate forces for a given torsion
vector<Coords3D> PeriodicTorsionForce::calculateForces(const vector<Coords3D>& atomPositions, const TorsionParameters& params) {
    // Assuming atomPositions contains positions for the four atoms involved in the torsion
    double torsionAngle = torsion_angle(atomPositions[0], atomPositions[1], atomPositions[2], atomPositions[3]);

    // Calculate energy derivative with respect to the torsion angle
    double dEdPhi = params.k * params.periodicity * sin(params.periodicity * torsionAngle - params.phase);

    Coords3D r1 = atomPositions[0], r2 = atomPositions[1], r3 = atomPositions[2], r4 = atomPositions[3];

    Coords3D b1 = r2 - r1;
    Coords3D b2 = r3 - r2;
    Coords3D b3 = r4 - r3;

    Coords3D n1 = b1.cross(b2).normalize();
    Coords3D n2 = b2.cross(b3).normalize();
    Coords3D m1 = n1.cross(b2).normalize();

    // Derivatives of the torsion angle with respect to the positions of atoms
    double x = n1.dot(n2);
    double y = m1.dot(n2);
    double angle = atan2(y, x);

    // Convert dEdPhi to forces
    Coords3D f1, f2, f3, f4;
    // Compute vectors perpendicular to bonds
    Coords3D p1 = b2.cross(b1).normalize();
    Coords3D p2 = b1.cross(b2).normalize() + b3.cross(b2).normalize();
    Coords3D p3 = b2.cross(b3).normalize();

    double factor = -dEdPhi / sin(angle); // Adjusted for correct units and direction

    f1 = p1 * factor;
    f2 = p2 * factor;
    f3 = -p2 * factor; // p3 should be equal in magnitude but opposite in direction to p2 for internal consistency
    f4 = -p3 * factor;

    vector<Coords3D> forces = { f1, f2, f3, f4 };

    return forces;
}