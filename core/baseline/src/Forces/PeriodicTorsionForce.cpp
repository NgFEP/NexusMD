#include "PeriodicTorsionForce.h"
#include "Coords3D.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "SystemXMLParser.h" //this contains definitions for TorsionParameters

#ifndef PI
#define PI 3.14159265358979323846
#endif

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
double PeriodicTorsionForce::torsion_angle(const Coords3D& r1, const Coords3D& r2, const Coords3D& r3, const Coords3D& r4) {
    Coords3D b1 = r2 - r1;
    Coords3D b2 = r3 - r2;
    Coords3D b3 = r4 - r3;

    Coords3D n1 = b1.cross(b2).normalize();
    Coords3D n2 = b2.cross(b3).normalize();
    Coords3D m1 = n1.cross(b2).normalize();

    double x = n1.dot(n2);
    double y = m1.dot(n2);
    double angle = atan2(y, x);

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