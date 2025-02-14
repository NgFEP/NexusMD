//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>

#include "PeriodicTorsionForce.h"
#include "Coords3D.h"
#include <cmath>
#include <iomanip>

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

// To calculate the distance between two points
double PeriodicTorsionForce::atom_distance(const Coords3D& coords1, const Coords3D& coords2) {
    Coords3D diff = coords2 - coords1;
    return sqrt(diff.dot(diff));
}


// To calculate torsion angle (in radians) between four points
double PeriodicTorsionForce::torsion_angle(const Coords3D& p1, const Coords3D& p2, const Coords3D& p3, const Coords3D& p4, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    // Calculate vectors between atoms

    Coords3D v0 = PeriodicBoundaryCondition::minimumImageVector(p1, p2, boxInfo);
    Coords3D v1 = PeriodicBoundaryCondition::minimumImageVector(p3, p2, boxInfo);
    Coords3D v2 = PeriodicBoundaryCondition::minimumImageVector(p3, p4, boxInfo);


    // Cross products to get normals to planes formed by atoms
    Coords3D cp0 = v0.cross(v1);
    Coords3D cp1 = v1.cross(v2);

    Coords3D m1 = cp0.cross(cp1).normalize();


    // Normalizing cross products to use for angle calculation
    double cosangle = cp0.normalize().dot(cp1.normalize());
    double angle;

    // To check if cosangle is near the singularity of acos()
    if (cosangle > 0.99f || cosangle < -0.99f) {// no f to make sure double precision is applied
        Coords3D cross_prod = cp0.cross(cp1);
        double scale = cp0.dot(cp0) * cp1.dot(cp1);
        angle = asin(sqrt(cross_prod.dot(cross_prod) / scale));
        if (cosangle < 0) {
            angle = PI - angle;
        }
    }
    else {
        angle = acos(cosangle);
    }

    // To adjust angle based on the direction relative to the v0 vector
    angle = (v0.dot(cp1) >= 0 ? angle : -angle);

    return angle;
}



vector<Coords3D> PeriodicTorsionForce::calculateForces(const vector<Coords3D>& atomPositions, const PTorsionParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    // Assuming atomPositions contains positions for the four atoms involved in the torsion
    double torsionAngle = torsion_angle(atomPositions[0], atomPositions[1], atomPositions[2], atomPositions[3], boxInfo);

    // To calculate energy derivative with respect to the torsion angle
    double deltaAngle = params.periodicity * torsionAngle - params.phase;
    // New line for energy calculation
    double energy = params.k * (1.0 + cos(deltaAngle));
    totalPEnergy += energy;
    double dEdPhi = -params.k * params.periodicity * sin(deltaAngle);

    Coords3D r1 = atomPositions[0], r2 = atomPositions[1], r3 = atomPositions[2], r4 = atomPositions[3];

    // To calculate vectors between atoms


    Coords3D v0 = PeriodicBoundaryCondition::minimumImageVector(r1, r2, boxInfo);
    Coords3D v1 = PeriodicBoundaryCondition::minimumImageVector(r3, r2, boxInfo);
    Coords3D v2 = PeriodicBoundaryCondition::minimumImageVector(r3, r4, boxInfo);

    // To calculate cross products
    Coords3D cp0 = v0.cross(v1);
    double uuu1= v0[1] * v1[2] - v0[2] * v1[1];
    double uuu2 = v0[1] * v1[2] ;

    double uuu3 = - v0[2] * v1[1];

    Coords3D cp1 = v1.cross(v2);

    // To calculate norms of cross products and vectors
    double normCross1 = cp0.dot(cp0);
    double normSqrBC = v1.dot(v1);
    double normBC = sqrt(normSqrBC);
    double normCross2 = cp1.dot(cp1);
    double dp = 1.0 / normSqrBC;


    // Force factors calculation
    double factor1 = (-dEdPhi * normBC) / normCross1;
    double factor2 = v0.dot(v1) * dp;
    double factor3 = v2.dot(v1) * dp;
    double factor4 = (dEdPhi * normBC) / normCross2;


    // To calculate forces based on cross products and factors
    Coords3D force1 = cp0 * factor1;
    Coords3D force4 = cp1 * factor4;
    Coords3D s = force1 * factor2 - force4 * factor3;

    Coords3D force2 = s - force1;
    Coords3D force3 = -s - force4;


    if (abs(force1[1]) > 200) {
        cout << "";
    }
    if (abs(force2[1]) > 200) {
        cout << "";
    }
    if (abs(force3[1]) > 200) {
        cout << "";
    }
    if (abs(force4[1]) > 200) {
        cout << "";
    }
     vector<Coords3D> forces = { force1, force2, force3, force4 };

    return forces;
}
